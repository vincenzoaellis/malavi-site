/* =============================================================================
 * sequence-check.mjs
 * -----------------------------------------------------------------------------
 * The submit page's browser-side sequence checker.
 *
 * WHAT THIS IS FOR
 *   A submitter pastes one cytochrome b sequence and wants to know, before they
 *   write to a curator: is this already a named MalAvi lineage, or is it new?
 *
 * THE ONE ERROR IT MUST NOT MAKE
 *   Telling someone their sequence is new when it is already in MalAvi. That
 *   sends a duplicate name into a paper and into GenBank, where it is very hard
 *   to undo. Everything below is arranged around not doing that:
 *
 *   - It compares against EVERY lineage in the release, never a sample.
 *   - It is exact and deterministic: seed-and-extend over a k-mer index, then a
 *     full ungapped comparison. No heuristic scoring, no cutoffs that could
 *     silently drop a true match.
 *   - IUPAC ambiguity codes compare by base set, so a query with an N or an R
 *     at a position still matches a reference that has a concrete base there.
 *     An ambiguous site is not a difference, which is also what the submission
 *     guidance tells people.
 *   - When it cannot find a match it says so plainly and points at BLAST,
 *     rather than guessing at a nearest neighbour.
 *
 * It is a helper, not a verdict. A curator confirms every new lineage. The page
 * copy says so, and this module should never be made to sound more certain.
 *
 * Pure functions with no DOM access, so the whole thing runs under Node against
 * the real release. See tests/test_sequence_check.mjs.
 * ============================================================================= */

/* The MalAvi barcode window. Reference sequences are all stored at exactly this
   length; `min_length` is the shortest unambiguous run we will name from. Both
   are carried in the index payload (from config/project.yml) rather than being
   written here, so the stated rule and the applied rule cannot drift apart. */

/* ---- IUPAC nucleotide codes ------------------------------------------------
   Each code maps to the set of concrete bases it can stand for. Two positions
   are "compatible" when their base sets overlap -- that is what makes an N
   match everything and an R (A or G) match an A but not a C. */
const IUPAC = {
  A: "A", C: "C", G: "G", T: "T", U: "T",
  R: "AG", Y: "CT", S: "CG", W: "AT", K: "GT", M: "AC",
  B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT"
};

/* Precomputed compatibility table: COMPATIBLE[a][b] is true when codes a and b
   share at least one concrete base. Building it once keeps the inner comparison
   loop to a pair of object lookups. */
const COMPATIBLE = (() => {
  const codes = Object.keys(IUPAC);
  const table = {};
  for (const a of codes) {
    table[a] = {};
    for (const b of codes) {
      table[a][b] = [...IUPAC[a]].some((base) => IUPAC[b].includes(base));
    }
  }
  return table;
})();

const COMPLEMENT = {
  A: "T", T: "A", U: "A", G: "C", C: "G",
  R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K",
  B: "V", V: "B", D: "H", H: "D", N: "N", "-": "-"
};

/** Reverse-complement a sequence, preserving ambiguity codes. */
export function reverseComplement(seq) {
  let out = "";
  for (let i = seq.length - 1; i >= 0; i--) out += COMPLEMENT[seq[i]] || "N";
  return out;
}

/**
 * Pull a single nucleotide sequence out of whatever the user pasted.
 *
 * Accepts a bare sequence or FASTA. Drops header lines, whitespace, digits and
 * punctuation (so sequence copied out of a paper with line numbers still
 * works), uppercases, and keeps only letters and gap characters. Anything else
 * that survives is reported to the user as ignored rather than silently
 * dropped -- see countContent.
 */
export function cleanSequence(raw) {
  if (raw == null) return "";
  return String(raw)
    .split(/\r?\n/)
    .filter((line) => !line.trim().startsWith(">"))
    .join("")
    .toUpperCase()
    .replace(/[\s0-9.*]/g, "")
    .replace(/[^A-Z-]/g, "");
}

/** Break a cleaned sequence down into the categories the checker reports on. */
export function countContent(sequence) {
  const unambiguous = (sequence.match(/[ACGTU]/g) || []).length;
  const ambiguous = (sequence.match(/[RYSWKMBDHVN]/g) || []).length;
  const gaps = (sequence.match(/-/g) || []).length;
  return {
    unambiguous,
    ambiguous,
    gaps,
    // Letters that are not nucleotide codes at all (E, F, P, ... ) -- usually a
    // sign that protein sequence or stray prose was pasted in.
    invalid: sequence.length - unambiguous - ambiguous - gaps,
    length: sequence.length
  };
}

/* ---- the reference index ---------------------------------------------------
   Built once from docs/assets/data/lineage_sequences.json, then reused for
   every query. The k-mer index maps each exact k-length word in every
   reference to the (entry, position) pairs where it occurs, which is what lets
   a query find its candidate references without comparing against all 5,000+.

   k = 16 is chosen so that the index stays small while guaranteeing candidate
   recovery for anything close enough to matter: a 470 bp query differing from a
   reference at d positions must contain an exact run of at least
   floor((470 - d) / (d + 1)) matching bases, which stays above 16 for d up to
   about 26 differences. Real "is this a new lineage?" cases differ at a handful
   of positions, far inside that bound. Anything more divergent is not a
   near-match question and is sent to BLAST. */
const K = 16;

/* Seeds are taken every SEED_STRIDE bases rather than at every position. The
   whole point of a seed is to nominate candidates for exact scoring, and
   overlapping seeds nominate the same candidates over and over: at stride 1 a
   single query costs millions of index lookups, because conserved 16-mers occur
   in thousands of references at once.

   Stride 4 leaves a 470 bp query with ~114 seeds. A reference differing from it
   at d positions loses at most 4 seeds per difference (a difference falls
   inside 16/4 = 4 sampled windows), so d would have to exceed about 28 before
   every seed is destroyed -- far beyond the handful of differences that a
   "is this a new lineage?" question ever involves, and anything more divergent
   than that is a BLAST question, which is where the checker sends it. */
const SEED_STRIDE = 4;

export function buildIndex(payload) {
  const entries = payload.entries.map((e) => ({
    ...e,
    // Gaps are alignment padding, not biology. Comparison happens on the
    // ungapped string, with the gapped original kept for display.
    ungapped: e.seq.replace(/-/g, "")
  }));

  /* References are indexed at every position; only the QUERY is strided. That
     way a seed taken at any query offset still finds its reference, whatever
     the alignment offset between them happens to be. */
  const kmers = new Map();
  entries.forEach((entry, entryIndex) => {
    const s = entry.ungapped;
    for (let pos = 0; pos + K <= s.length; pos++) {
      const word = s.substr(pos, K);
      // A word containing an ambiguity code cannot be used as an exact seed,
      // because the query may legitimately carry a different code there.
      if (/[^ACGT]/.test(word)) continue;
      let bucket = kmers.get(word);
      if (!bucket) kmers.set(word, (bucket = []));
      bucket.push({ entryIndex, pos });
    }
  });

  return {
    entries,
    kmers,
    windowLength: payload.window_length,
    minLength: payload.min_length,
    release: payload.release,
    nLineages: payload.n_lineages
  };
}

/**
 * Count differences between a query and a reference at a fixed alignment
 * offset, over the region where they overlap.
 *
 * `offset` is where the query's first base sits within the reference; it may be
 * negative when the query extends past the reference's 5' end. Positions where
 * either side is ambiguous but compatible are NOT counted as differences.
 */
export function compareAtOffset(query, reference, offset) {
  const start = Math.max(0, -offset);
  const end = Math.min(query.length, reference.length - offset);
  let differences = 0;
  let compared = 0;
  let ambiguousSites = 0;

  for (let i = start; i < end; i++) {
    const q = query[i];
    const r = reference[i + offset];
    const qSet = IUPAC[q];
    const rSet = IUPAC[r];
    // A gap on either side is alignment padding; skip it entirely rather than
    // scoring it, and do not count it toward the compared length.
    if (!qSet || !rSet) continue;
    compared++;
    if (COMPATIBLE[q][r]) {
      if (qSet.length > 1 || rSet.length > 1) ambiguousSites++;
    } else {
      differences++;
    }
  }
  return { differences, compared, ambiguousSites, offset };
}

/** Keep one candidate per reference entry, preserving order. A candidate can be
    produced twice — once per orientation, or by both the exact scan and the
    seeded search. */
function dedupeByEntry(candidates) {
  const seen = new Set();
  const out = [];
  for (const candidate of candidates) {
    if (seen.has(candidate.entry)) continue;
    seen.add(candidate.entry);
    out.push(candidate);
  }
  return out;
}

/**
 * Find the best-matching reference lineages for a query.
 *
 * Tries the query in both orientations. Returns every reference tied at the
 * lowest difference count, because a sequence shared by several lineage names
 * must not be reported as just one of them.
 */
export function findMatches(index, querySequence, { maxCandidates = 400 } = {}) {
  const cleanedQuery = querySequence.replace(/-/g, "");
  if (!cleanedQuery) return null;

  /* Exact path. When the query carries no ambiguity codes, an exact match can
     be found by straight string containment, and -- importantly -- EVERY such
     match can be found, by scanning the whole reference set.

     That exhaustiveness is the point, not the speed. A short sequence can sit
     inside many different lineages: in this release the 210 bp RBQ16 is
     contained in 13 of them. Reporting only the first one found would tell a
     submitter "this is RBQ16" when the sequence cannot distinguish RBQ16 from
     twelve others. The scan costs a couple of milliseconds and removes that
     whole class of over-claim. */
  /* Candidates are ranked by how many bases they actually explain, not by
     difference count alone. Ranking on differences alone is actively dangerous:
     a 479 bp sequence one base away from LBPIP1 happens to contain the 210 bp
     RBQ16 exactly, and "0 differences over 210 bp" would beat "1 difference
     over 479 bp" -- reporting a genuinely new lineage as already named. Matched
     bases gets this right: 478 beats 210. */
  const matchedBases = (candidate) => candidate.compared - candidate.differences;
  const better = (a, b) => {
    if (!b) return true;
    if (matchedBases(a) !== matchedBases(b)) return matchedBases(a) > matchedBases(b);
    return a.differences < b.differences;
  };

  const exactCandidates = [];
  if (!/[^ACGT]/.test(cleanedQuery)) {
    for (const [orientation, query] of [
      ["forward", cleanedQuery],
      ["reverse", reverseComplement(cleanedQuery)]
    ]) {
      for (const entry of index.entries) {
        const reference = entry.ungapped;
        // Either the query sits inside the reference (a truncated submission)
        // or the reference sits inside the query (the submitter sequenced more
        // than the barcode window). Both are exact matches.
        if (reference.includes(query) || query.includes(reference)) {
          exactCandidates.push({
            entry,
            orientation,
            differences: 0,
            compared: Math.min(query.length, reference.length),
            ambiguousSites: 0,
            offset: 0
          });
        }
      }
    }
  }

  /* Short-circuit only when a match explains the WHOLE query with no
     differences. No other candidate can beat that, because matched bases can
     never exceed the query length -- so skipping the search is provably safe
     rather than merely quick. This is the common case: someone pasting a
     sequence straight out of MalAvi or GenBank. */
  const perfect = exactCandidates.filter((c) => c.compared === cleanedQuery.length);
  if (perfect.length) {
    const ties = dedupeByEntry(perfect);
    return {
      best: ties[0],
      ties,
      names: [...new Set(ties.flatMap((t) => t.entry.names))].sort()
    };
  }

  let best = null;
  for (const candidate of exactCandidates) if (better(candidate, best)) best = candidate;
  // Every candidate we actually scored, kept so ties can be gathered from this
  // set rather than by re-scanning the whole reference (which would cost
  // millions of character comparisons per query and stall the page).
  const scored = [];

  for (const [orientation, query] of [
    ["forward", cleanedQuery],
    ["reverse", reverseComplement(cleanedQuery)]
  ]) {
    /* Seed: collect votes for (reference, offset) pairs from shared k-mers.
       A vote means "this reference lines up with the query at this offset". */
    const votes = new Map();
    for (let pos = 0; pos + K <= query.length; pos += SEED_STRIDE) {
      const word = query.substr(pos, K);
      if (/[^ACGT]/.test(word)) continue;
      const bucket = index.kmers.get(word);
      if (!bucket) continue;
      for (const hit of bucket) {
        // offset = where query position 0 sits inside the reference
        const offset = hit.pos - pos;
        const key = hit.entryIndex + ":" + offset;
        votes.set(key, (votes.get(key) || 0) + 1);
      }
    }
    if (votes.size === 0) continue;

    /* Extend: score the best-supported candidates exactly. Sorting by vote
       count first means the true match is scored early; the cap only ever
       discards candidates that shared fewer exact words than 400 others. */
    const ranked = [...votes.entries()]
      .sort((a, b) => b[1] - a[1])
      .slice(0, maxCandidates);

    for (const [key] of ranked) {
      const sep = key.lastIndexOf(":");
      const entryIndex = Number(key.slice(0, sep));
      const offset = Number(key.slice(sep + 1));
      const entry = index.entries[entryIndex];
      const result = compareAtOffset(query, entry.ungapped, offset);
      if (result.compared === 0) continue;

      const candidate = { entry, orientation, ...result };
      scored.push(candidate);
      if (better(candidate, best)) best = candidate;
    }
  }

  if (!best) return null;

  /* Gather every candidate tied with the winner, so a sequence carried by more
     than one lineage name is reported under all of them. Exact matches came
     from an exhaustive scan; seeded matches are drawn from `scored`, and
     identical sequences always seed identically, so the relevant entries are
     present whenever the winner is. */
  const ties = dedupeByEntry(
    [...exactCandidates, ...scored].filter(
      (c) => c.differences === best.differences && c.compared === best.compared
    )
  );

  const names = [...new Set(ties.flatMap((t) => t.entry.names))].sort();
  return { best, ties, names };
}

/**
 * The full check: content validation plus identity, as a plain data structure.
 * Rendering lives in the page; this returns only findings.
 *
 * `verdict` is one of:
 *   "stop"    something is wrong with the input, or it is too short to judge
 *   "known"   an exact match to one or more named lineages
 *   "new"     no exact match, but a close relative was found
 *   "unknown" nothing close enough to say anything useful -- go to BLAST
 */
export function checkSequence(index, raw) {
  const sequence = cleanSequence(raw);
  const content = countContent(sequence);
  const checks = [];

  if (!sequence) {
    return {
      verdict: "stop",
      title: "Nothing to check",
      message: "Paste one nucleotide sequence, with or without a FASTA header.",
      content,
      checks
    };
  }

  if (content.invalid > 0) {
    checks.push({
      state: "fail",
      label: "Content",
      text:
        `${content.invalid} character${content.invalid === 1 ? "" : "s"} that are not ` +
        "nucleotide codes were found and ignored. Check that this is a DNA sequence."
    });
  } else {
    checks.push({ state: "pass", label: "Content", text: "Nucleotide codes only." });
  }

  if (content.ambiguous > 0) {
    checks.push({
      state: "warn",
      label: "Ambiguity",
      text:
        `${content.ambiguous} ambiguous position${content.ambiguous === 1 ? "" : "s"} ` +
        "(N, R, Y…). These are not counted as differences: a difference at an " +
        "ambiguous site does not make a new lineage."
    });
  }

  if (content.gaps > 0) {
    checks.push({
      state: "warn",
      label: "Gaps",
      text:
        `${content.gaps} gap character${content.gaps === 1 ? "" : "s"}. Cytochrome b ` +
        "should not carry indels — please check for a frameshift."
    });
  }

  const minLength = index.minLength;
  if (content.unambiguous < minLength) {
    checks.push({
      state: "fail",
      label: "Length",
      text:
        `${content.unambiguous} bp of unambiguous sequence — below the ${minLength} bp ` +
        "minimum. Too short to show it is genuinely distinct."
    });
  } else {
    checks.push({
      state: "pass",
      label: "Length",
      text: `${content.unambiguous} bp of unambiguous sequence.`
    });
  }

  /* Identity is assessed even for short sequences: telling someone their 300 bp
     read already matches a named lineage is useful and safe. What a short
     sequence cannot support is naming something NEW, which the verdict below
     refuses to do. */
  const match = findMatches(index, sequence);

  if (!match) {
    checks.push({
      state: "warn",
      label: "Identity",
      text:
        "No MalAvi lineage shares a substantial exact stretch with this sequence. " +
        "That can mean it is very divergent, is not avian haemosporidian cytochrome b, " +
        "or is not in the standard barcode region."
    });
    return {
      verdict: "unknown",
      title: "No close match found",
      message:
        "Nothing in the release lines up with this sequence closely enough to judge. " +
        "Run it through BLAST for a full search, and send it to us either way.",
      content,
      match: null,
      checks
    };
  }

  const { best, names } = match;
  const nameList = names.join(", ");

  checks.push({
    state: best.orientation === "reverse" ? "warn" : "pass",
    label: "Orientation",
    text:
      best.orientation === "reverse"
        ? "Reverse-complemented. Flipped for you — no need to resubmit."
        : "Forward, as expected."
  });

  if (best.differences === 0) {
    /* An exact match means the sequence is certainly already in MalAvi, which
       is the safe and important half of the answer. WHICH lineage it is can
       still be undecidable: a short sequence may sit inside many lineages, and
       claiming one of them would be a guess dressed as a result. */
    const ambiguousIdentity = names.length > 1;

    checks.push({
      state: "pass",
      label: "Genus",
      text: [...new Set(match.ties.flatMap((t) => t.entry.genus))].join(" / ")
    });
    checks.push({
      state: ambiguousIdentity ? "warn" : "pass",
      label: "Identity",
      text: ambiguousIdentity
        ? `Exact match over ${best.compared} bp, but to ${names.length} lineages: ${nameList}. ` +
          "This sequence is too short to tell them apart."
        : `Exact match to ${nameList} over ${best.compared} bp.`
    });
    checks.push({
      state: "pass",
      label: "Naming",
      text: "No new name required — this sequence is already in the database."
    });

    if (ambiguousIdentity) {
      return {
        verdict: "known",
        title: `Already in MalAvi — matches ${names.length} lineages`,
        message:
          `This sequence matches ${names.length} MalAvi lineages exactly over ` +
          `${best.compared} bp: ${nameList}. It is already in the database, so it does not ` +
          "need a new name, but it is too short to say which of these it is. A longer " +
          "sequence, or a curator, can settle that.",
        content,
        match,
        checks
      };
    }

    return {
      verdict: "known",
      title: `Already in MalAvi — this is ${nameList}`,
      message:
        `An exact match over ${best.compared} bp to ${nameList}` +
        (best.entry.acc.length ? ` (GenBank ${best.entry.acc.join(", ")})` : "") +
        ". Use the existing name; no new lineage is needed.",
      content,
      match,
      checks
    };
  }

  const identity = (100 * (1 - best.differences / best.compared)).toFixed(1);

  if (content.unambiguous < minLength) {
    checks.push({
      state: "fail",
      label: "Identity",
      text:
        `Closest is ${nameList} at ${best.differences} difference` +
        `${best.differences === 1 ? "" : "s"} over ${best.compared} bp, but the sequence ` +
        "is too short to name a new lineage from."
    });
    return {
      verdict: "stop",
      title: "Too short to name",
      message:
        `Only ${content.unambiguous} bp of unambiguous sequence. Its closest match is ` +
        `${nameList} (${identity}% identical), but we would ask you to re-sequence before ` +
        "naming this as new. Do still get in touch if that is not possible.",
      content,
      match,
      checks
    };
  }

  checks.push({
    state: "warn",
    label: "Genus",
    text:
      `Probably ${best.entry.genus.join(" / ")}, inferred from the closest match at ` +
      `${identity}% identity. A curator confirms this.`
  });
  checks.push({
    state: "pass",
    label: "Identity",
    text:
      `Closest is ${nameList} at ${best.differences} difference` +
      `${best.differences === 1 ? "" : "s"} over ${best.compared} bp.`
  });
  checks.push({
    state: "warn",
    label: "Naming",
    text: "Tell us the host species and we will give you the next free number for its acronym."
  });

  return {
    verdict: "new",
    title: "Looks like a new lineage",
    message:
      `Closest match is ${nameList} — ${best.differences} difference` +
      `${best.differences === 1 ? "" : "s"} over ${best.compared} bp (${identity}% identical). ` +
      "One difference is enough to make it new, so this needs a name. " +
      "A curator confirms every new lineage before it is added.",
    content,
    match,
    checks
  };
}
