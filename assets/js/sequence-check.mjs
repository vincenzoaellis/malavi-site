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
 *     at a position still matches a reference that has a concrete base there --
 *     and so does the reverse, a reference with an N against a submitter's clean
 *     read. An ambiguous site is not a difference, which is also what the
 *     submission guidance tells people. Ambiguous positions are handled the same
 *     way at every stage: in the comparison, in the exact scan, and in seeding.
 *   - When it cannot find a match it says so plainly and points at BLAST,
 *     rather than guessing at a nearest neighbour. And when it could not
 *     properly look -- a read so ambiguous that no stretch of it is definite
 *     enough to search with -- it says THAT, instead of reporting no match.
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

/* An ambiguous position is compared by base set (above), but it cannot be looked
   up in an exact k-mer index by letter -- the index holds concrete words only.
   Rather than throw such a window away, we expand it into the concrete words it
   could stand for and look each one up. This is what stops a mediocre read from
   blinding the checker entirely; see EXPANDING AMBIGUOUS SEEDS below.

   The budget is on the number of realizations, not the number of codes, because
   that is what the work is proportional to: one N is 4 words, two Ns are 16, one
   R is 2, an R and an N together are 8. Sixteen keeps the worst case at ~114
   windows x 16 lookups per orientation, negligible beside the comparison work
   that follows. */
const MAX_SEED_EXPANSION = 16;

/**
 * Expand a k-mer into the concrete ACGT words it is compatible with.
 *
 * Returns the single word unchanged when it is already concrete, an array of
 * realizations when it is ambiguous but within MAX_SEED_EXPANSION, and null when
 * it is too ambiguous to expand or contains something that is not a nucleotide
 * code at all (a gap, a stray letter). Null means "this window cannot be used as
 * a seed", which is the same answer the module gave for every ambiguous window
 * before.
 */
export function expandSeedWord(word, maxVariants = MAX_SEED_EXPANSION) {
  // Fast path: the overwhelming majority of windows are already concrete.
  if (!/[^ACGT]/.test(word)) return [word];

  /* Count the realizations before building any of them, so an unusable window
     costs a short loop rather than a large intermediate array. */
  let variants = 1;
  for (const code of word) {
    const bases = IUPAC[code];
    if (!bases) return null; // not a nucleotide code -- never seedable
    variants *= bases.length;
    if (variants > maxVariants) return null;
  }

  /* Build the cross product one position at a time. */
  let words = [""];
  for (const code of word) {
    const bases = IUPAC[code];
    const next = [];
    for (const prefix of words) {
      for (const base of bases) next.push(prefix + base);
    }
    words = next;
  }
  return words;
}

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

/* ---- CORROBORATION REQUIRED OF AN EXACT MATCH -------------------------------
   containmentOffset() calls an alignment exact when every overlapping position
   is COMPATIBLE. Compatible is not the same as observed: N is compatible with
   every base, so a query shorter than a reference's run of ambiguity codes
   slides wholesale into that run and comes back "0 differences over its whole
   length" having been told nothing at all. BUTJAM13 carries 175 N of 479 bp,
   and that one reference tied with the genuine match for 149 of 150 measured
   150 bp queries (malavi_sanger, METHODS_draft.md 5H.8).

   So an exact alignment must additionally be CORROBORATED: it must contain at
   least one unbroken stretch of positions where both sides are a single
   concrete base and they agree.

   K is reused as that length rather than a new number being invented. It is
   already this module's standard for "a stretch long enough to nominate a
   candidate" -- the seeded path accepts nothing shorter -- so the exact path is
   now held to the same evidentiary bar, by construction rather than by
   coincidence.

   Measured separation in the 2026-03-23 release: every genuine member of a
   reported equivalence class scores at least 25, because 25 is the shortest
   concrete run any of the 191 ambiguity-carrying references has anywhere; every
   wildcard match scores 0 to 2. No reference in the release can be excluded by
   this test at its own true offset. */
const EXACT_MATCH_CORROBORATION = K;

/* ---- EXPANDING AMBIGUOUS SEEDS ---------------------------------------------
   The reasoning above is about seeds destroyed by DIFFERENCES, and it is correct
   about those. Ambiguity codes destroy seeds by a different mechanism: they are
   excluded from the index and from seeding outright, because they cannot be
   matched as exact letters. An ambiguity code every s positions leaves clean
   runs of exactly s - 1 bases, so at s <= 16 a K = 16 seed has nowhere to sit
   and NOT ONE seed survives anywhere in the query.

   That was measured (malavi_sanger, METHODS_draft.md 5H.2): with an N every
   17 bp all 60 test queries were named correctly; with an N every 16 bp all 60
   returned "no close match" while a position-wise comparison named the source
   lineage from a median of 478 informative positions. An N every 16 bp is 30
   ambiguous calls in 479 bp -- a mediocre but ordinary Sanger trace. One base of
   spacing separated full function from total failure, with no warning.

   Two changes remove that cliff:
     1. Ambiguous query windows are expanded into their concrete realizations and
        each is looked up (expandSeedWord). One N per window costs 4 lookups, so
        the s = 16 case seeds normally again, and spacings below it degrade
        gradually instead of falling off an edge.
     2. When a query really is too ambiguous to seed at all, checkSequence says
        exactly that rather than asserting that no lineage matches -- because
        that assertion was false.
   Only query-side seeding changed; the reference index is built as before. */

/* ---- THE SITE PROFILE -------------------------------------------------------
   Every reference is stored padded to exactly the 479 bp alignment window, so
   column c of one sequence is homologous to column c of every other. That makes
   a per-column tally of the bases MalAvi has ever seen possible, and the tally
   is what lets the checker say something a plain nearest-neighbour comparison
   cannot: "the base you have here has never been seen at this position in any
   lineage."

   Counts are weighted by how many LINEAGE NAMES carry a sequence, not by how
   many distinct sequences do, because the question a submitter is asking is
   about lineages.

   It is built here, from the index itself, rather than shipped as its own data
   file, so it can never describe a different release from the sequences beside
   it. Cost is one pass over 5,359 x 479 characters -- far less than the k-mer
   index built in the same function.

   Columns with poor coverage are excluded from the novelty judgement. Column 1
   is covered by only 2,337 of the 5,367 lineages, because both standard inner
   forward primers end one base into the window (see the primer frame reference
   in reference/cytb_primer_frame_reference/), and "never seen before" means
   very little where most lineages have nothing recorded at all. The 50%
   threshold is the same one malaviR's lineage_screen() uses, and in this
   release it excludes column 1 and nothing else. */
const MIN_SITE_COVERAGE = 0.5;

const BASE_INDEX = { A: 0, C: 1, G: 2, T: 3 };

/* ---- READING FRAME AND STOP CODONS -----------------------------------------
   The MalAvi window is in frame: column 1 is a first codon position, so codons
   are columns 1-3, 4-6, and so on. malaviR relies on the same convention
   (R/internal.R, .qc_codon_position), and it holds -- of the 5,359 sequences in
   the 2026-03-23 release only 18 carry a stop codon read this way.

   The genetic code is NCBI table 4, mold/protozoan/coelenterate mitochondrial,
   which is the right one for Plasmodium, Haemoproteus and Leucocytozoon. It
   differs from the standard code in exactly the way that matters here: TGA
   codes tryptophan rather than stop, so TAA and TAG are the only stops. Reading
   these sequences under the vertebrate mitochondrial code (table 2) invents
   stop codons that are not there; malaviR carries a regression test for that
   same mistake. */
const STOP_CODONS = new Set(["TAA", "TAG"]);

export function buildIndex(payload) {
  const entries = payload.entries.map((e) => {
    // Gaps are alignment padding, not biology. Comparison happens on the
    // ungapped string, with the gapped original kept for display.
    const ungapped = e.seq.replace(/-/g, "");
    return {
      ...e,
      ungapped,
      /* Whether this reference carries an ambiguity code, computed once here so
         the exact scan can route each reference to the right comparison: 191 of
         the 5,359 entries in the 2026-03-23 release do, and only those need the
         slower compatibility-aware containment test. */
      hasAmbiguity: /[^ACGT]/.test(ungapped)
    };
  });

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

  /* Per-column base counts over the whole release. See THE SITE PROFILE. */
  const windowLength = payload.window_length;
  const siteCounts = new Int32Array(windowLength * 4); // A, C, G, T per column
  const siteCoverage = new Int32Array(windowLength);   // lineages with a base there
  for (const entry of entries) {
    const weight = entry.names.length;
    const gapped = entry.seq;
    const columns = Math.min(gapped.length, windowLength);
    for (let column = 0; column < columns; column++) {
      // Gaps and ambiguity codes are not evidence that a base was observed, so
      // they contribute to neither the counts nor the coverage.
      const base = BASE_INDEX[gapped[column]];
      if (base === undefined) continue;
      siteCounts[column * 4 + base] += weight;
      siteCoverage[column] += weight;
    }
  }

  /* Every lineage name in the release, sorted. Used by suggestLineageName to
     work out which numbers a host acronym has already used. */
  const allNames = [...new Set(entries.flatMap((e) => e.names))].sort();

  return {
    entries,
    kmers,
    windowLength,
    minLength: payload.min_length,
    release: payload.release,
    nLineages: payload.n_lineages,
    siteCounts,
    siteCoverage,
    allNames
  };
}

/**
 * How many lineages carry `base` at 1-based alignment column `column`?
 *
 * Returns -1 when the column is outside the window, or is too poorly covered
 * across the release for its tally to mean anything (see MIN_SITE_COVERAGE).
 * Callers treat -1 as "no opinion" rather than as zero, which is the whole
 * point of separating it from a genuine count of nought.
 */
export function siteBaseCount(index, column, base) {
  if (column < 1 || column > index.windowLength) return -1;
  const slot = BASE_INDEX[base];
  if (slot === undefined) return -1;
  if (index.siteCoverage[column - 1] < MIN_SITE_COVERAGE * index.nLineages) return -1;
  return index.siteCounts[(column - 1) * 4 + slot];
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
  /* Positions where BOTH sides are a single concrete base. These are the only
     positions that carry POSITIVE evidence about which lineage this is:
     agreement at a position where either side is unresolved is a
     non-contradiction, not an observation, and conflating the two is what let a
     mostly-N reference pass for a match (see EXACT_MATCH_CORROBORATION). */
  let jointlyInformative = 0;
  /* The longest unbroken stretch of those positions on which the two sides also
     AGREE. A bare count cannot tell "corroborated by one solid stretch" apart
     from "one base here, one base there", and it is the stretch that decides
     whether an alignment was found on evidence. */
  let longestConcreteAgreement = 0;
  let currentAgreementRun = 0;

  for (let i = start; i < end; i++) {
    const q = query[i];
    const r = reference[i + offset];
    const qSet = IUPAC[q];
    const rSet = IUPAC[r];
    // A gap on either side is alignment padding; skip it entirely rather than
    // scoring it, and do not count it toward the compared length. It also
    // breaks the agreement run, because a run has to be contiguous to mean
    // anything.
    if (!qSet || !rSet) {
      currentAgreementRun = 0;
      continue;
    }
    compared++;

    const bothConcrete = qSet.length === 1 && rSet.length === 1;
    if (bothConcrete) jointlyInformative++;

    if (!COMPATIBLE[q][r]) {
      differences++;
      currentAgreementRun = 0;
      continue;
    }
    if (!bothConcrete) {
      // Compatible, but only because at least one side is unresolved.
      ambiguousSites++;
      currentAgreementRun = 0;
      continue;
    }
    // Two concrete bases, compatible -- therefore identical. Real evidence.
    currentAgreementRun++;
    if (currentAgreementRun > longestConcreteAgreement) {
      longestConcreteAgreement = currentAgreementRun;
    }
  }
  return {
    differences,
    compared,
    ambiguousSites,
    jointlyInformative,
    longestConcreteAgreement,
    offset
  };
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
 * Where does `inner` sit inside `outer`, comparing by IUPAC base set rather than
 * by letter? Returns the first such offset, or -1 if there is none.
 *
 * `ambiguous` says whether either string can carry a code outside A/C/G/T. When
 * neither does -- true for 5,168 of the 5,359 references, against a concrete
 * query -- set containment and literal containment are the same question, and
 * the engine's own substring search answers it far faster than a JS loop.
 *
 * This matters because the letter-wise test alone was wrong: a reference with an
 * N where the submitter has a real base is the submitter's OWN lineage, and
 * `includes()` cannot see it. See findMatches's exact path.
 *
 * `minAgreementRun` is the corroboration requirement (EXACT_MATCH_CORROBORATION):
 * an accepted offset must also carry an unbroken run of at least this many
 * positions where both sides are a single concrete base and agree. It is tested
 * INSIDE the scan rather than applied afterwards to the returned offset, because
 * a reference can contain the query at more than one offset -- and if an
 * uninformative one came first, filtering afterwards would throw away a genuine
 * match sitting further along. Default 0 preserves the plain containment
 * question for callers that only want that.
 */
export function containmentOffset(outer, inner, ambiguous, minAgreementRun = 0) {
  if (inner.length > outer.length) return -1;
  if (!ambiguous) {
    /* Neither string carries a code outside A/C/G/T, so every position of a
       literal containment is two concrete bases in agreement: the run is the
       whole of `inner`, and the requirement reduces to a length test. */
    if (inner.length < minAgreementRun) return -1;
    return outer.indexOf(inner);
  }

  const lastOffset = outer.length - inner.length;
  for (let offset = 0; offset <= lastOffset; offset++) {
    let matches = true;
    let currentRun = 0;
    let longestRun = 0;
    for (let i = 0; i < inner.length; i++) {
      const a = inner[i];
      const b = outer[i + offset];
      const aSet = IUPAC[a];
      const bSet = IUPAC[b];
      /* A character that is not a nucleotide code is not "compatible with
         everything" -- it cannot be part of an exact match at all. */
      if (!aSet || !bSet || !COMPATIBLE[a][b]) {
        matches = false;
        break;
      }
      if (aSet.length === 1 && bSet.length === 1) {
        currentRun++;
        if (currentRun > longestRun) longestRun = currentRun;
      } else {
        // Compatible only because one side is unresolved: no evidence here.
        currentRun = 0;
      }
    }
    if (matches && longestRun >= minAgreementRun) return offset;
  }
  return -1;
}

/**
 * How many of the seed windows a query would sample can actually be used.
 *
 * Reported for the better of the two orientations, since findMatches tries both.
 * This exists so that "we found nothing" can be told apart from "we could not
 * look" -- the difference between a sequence that is genuinely unlike anything in
 * MalAvi and one that is simply too ambiguous for an exact-word index. Only
 * called when a query fails to match, so it costs nothing in the normal case.
 */
export function countUsableSeeds(sequence) {
  const cleaned = sequence.replace(/-/g, "");

  const countOneOrientation = (query) => {
    let sampled = 0;
    let usable = 0;
    for (let pos = 0; pos + K <= query.length; pos += SEED_STRIDE) {
      sampled++;
      if (expandSeedWord(query.substr(pos, K))) usable++;
    }
    return { sampled, usable };
  };

  const forward = countOneOrientation(cleaned);
  const reverse = countOneOrientation(reverseComplement(cleaned));
  return forward.usable >= reverse.usable ? forward : reverse;
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
        //
        // Containment is judged by base set, not by letter. Using `includes()`
        // alone made the short-circuit below unsound: it missed every reference
        // carrying an ambiguity code where the query has a concrete base, and
        // then returned some OTHER lineage that happened to contain the query
        // literally. That produced 19 confidently mis-named queries in the
        // 2026-08-02 sweep -- "Already in MalAvi -- this is X", with the
        // submitter's actual lineage not among the names shown
        // (METHODS_draft.md 5H.5). The set-wise test makes `perfect` the
        // complete equivalence class, which is what the short-circuit assumes.
        //
        // Set-wise containment on its own over-shoots in the other direction --
        // a mostly-N reference absorbs any short query -- so the offset must
        // also be corroborated by a real stretch of agreeing concrete bases.
        // See EXACT_MATCH_CORROBORATION.
        const queryInside = containmentOffset(
          reference, query, entry.hasAmbiguity, EXACT_MATCH_CORROBORATION
        );
        const referenceInside =
          queryInside < 0
            ? containmentOffset(query, reference, entry.hasAmbiguity, EXACT_MATCH_CORROBORATION)
            : -1;
        if (queryInside >= 0 || referenceInside >= 0) {
          /* Score the hit with the same function every other candidate goes
             through, at the offset where it was found (negative when the
             reference sits inside the query). That keeps `compared`,
             `differences` and `ambiguousSites` consistent across both paths --
             `better()` compares candidates from each -- and reports the
             ambiguous sites honestly instead of assuming an exact match has
             none, which stopped being true once containment became set-wise. */
          const offset = queryInside >= 0 ? queryInside : -referenceInside;
          exactCandidates.push({
            entry,
            orientation,
            ...compareAtOffset(query, reference, offset)
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
  /* `differences === 0` is guaranteed by how these candidates were found -- they
     are scored at an offset already shown to be compatible at every position --
     but it is asserted rather than assumed, because the whole soundness of
     returning early rests on it. */
  const perfect = exactCandidates.filter(
    (c) => c.compared === cleanedQuery.length && c.differences === 0
  );
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
    /* Each candidate carries two vote counts, because expanded seeds and
       concrete seeds have to be ranked separately -- see THE TWO RANKINGS. */
    const votes = new Map();
    for (let pos = 0; pos + K <= query.length; pos += SEED_STRIDE) {
      /* An ambiguous window is expanded into the concrete words it could be and
         all of them are looked up, instead of the window being discarded. Two
         realizations of one window can never match the same reference position,
         so a window still contributes at most one vote per (reference, offset)
         and the vote counts keep meaning "how many sampled windows agree". */
      const word = query.substr(pos, K);
      const words = expandSeedWord(word);
      if (!words) continue;
      const isConcrete = words.length === 1 && words[0] === word;
      for (const realization of words) {
        const bucket = index.kmers.get(realization);
        if (!bucket) continue;
        for (const hit of bucket) {
          // offset = where query position 0 sits inside the reference
          const offset = hit.pos - pos;
          const key = hit.entryIndex + ":" + offset;
          let tally = votes.get(key);
          if (!tally) votes.set(key, (tally = { total: 0, concrete: 0 }));
          tally.total++;
          if (isConcrete) tally.concrete++;
        }
      }
    }
    if (votes.size === 0) continue;

    /* ---- THE TWO RANKINGS ---------------------------------------------------
       Extend: score the best-supported candidates exactly. Sorting by vote count
       first means the true match is scored early; the cap only ever discards
       candidates that shared fewer exact words than 400 others.

       Two rankings are taken and their union scored, because expanded seeds
       systematically DISADVANTAGE one particular candidate: the reference the
       query actually came from, when the query carries ambiguity codes of its
       own. The index holds no word containing an ambiguity code, so at exactly
       the windows where the query is ambiguous its own source cannot be voted
       for -- while every other reference, which has concrete bases there, can be.
       Ranking on expanded votes alone therefore pushed real sources (RBQ15,
       RBQ18, RBQ19 and others) out of the top 400 and reported near-misses
       instead.

       Ranking by concrete votes alone is exactly what this module did before
       expansion existed, and that ranking recognized every sequence in the
       release. So it is kept as a floor: whatever the old seeding would have
       scored is still scored, and expansion only ever ADDS candidates. */
    const rankedKeys = new Set();
    for (const [property, eligible] of [
      ["total", () => true],
      ["concrete", (tally) => tally.concrete > 0]
    ]) {
      const ranking = [...votes.entries()]
        .filter(([, tally]) => eligible(tally))
        .sort((a, b) => b[1][property] - a[1][property])
        .slice(0, maxCandidates);
      for (const [key] of ranking) rankedKeys.add(key);
    }

    for (const key of rankedKeys) {
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

/* ---- PLACING THE QUERY IN THE BARCODE WINDOW --------------------------------
   Everything below the identity question -- how much of the barcode region was
   actually read, whether the sequence reads through in frame, whether a base is
   one MalAvi has never seen at that position -- needs to know which alignment
   COLUMN each base of the submitted sequence sits in. The match gives us that
   for free: findMatches already worked out the offset at which the query lines
   up with its closest reference, and every reference is padded to the window,
   so the reference's own gap pattern converts an offset into columns.

   Bases that fall outside columns 1..479 are the case the checker used to say
   nothing about at all: someone who sequenced past the barcode region on either
   side was told the length of everything they pasted and left to assume all of
   it had been compared, when only the part inside the window ever is. */

/**
 * Work out the alignment column of every base of the submitted sequence.
 *
 * `query` is the cleaned, ungapped sequence as pasted (forward orientation);
 * `best` is the winning candidate from findMatches. Returns the oriented
 * sequence that was actually compared, the 1-based column of each of its
 * positions (which may be below 1 or above the window length, meaning the base
 * lies outside the barcode region), and a tally of how much fell where.
 */
export function placeInWindow(index, query, best) {
  const oriented = best.orientation === "reverse" ? reverseComplement(query) : query;

  /* Column of each ungapped position of the reference. A short lineage is
     stored padded -- RBQ18 is 133 bp of sequence inside 479 columns -- so this
     is not the identity mapping and cannot be assumed to be. */
  const referenceColumns = [];
  const gapped = best.entry.seq;
  for (let column = 0; column < gapped.length; column++) {
    if (gapped[column] !== "-") referenceColumns.push(column + 1);
  }

  const columns = new Array(oriented.length);
  let inside = 0;
  let before = 0;
  let after = 0;
  let firstColumn = 0;
  let lastColumn = 0;

  for (let i = 0; i < oriented.length; i++) {
    // Where this base sits along the reference's own ungapped sequence.
    const referencePosition = i + best.offset;
    let column;
    if (referencePosition >= 0 && referencePosition < referenceColumns.length) {
      column = referenceColumns[referencePosition];
    } else if (referencePosition < 0) {
      /* The query starts before the reference does. Step back through whatever
         padding the reference carries at its 5' end: those columns are still
         inside the barcode region, and only once the padding runs out is the
         base genuinely outside it. */
      column = referenceColumns[0] + referencePosition;
    } else {
      column = referenceColumns[referenceColumns.length - 1] +
        (referencePosition - referenceColumns.length + 1);
    }
    columns[i] = column;

    if (column < 1) {
      before++;
    } else if (column > index.windowLength) {
      after++;
    } else {
      inside++;
      if (!firstColumn || column < firstColumn) firstColumn = column;
      if (column > lastColumn) lastColumn = column;
    }
  }

  return { oriented, columns, inside, before, after, firstColumn, lastColumn };
}

/**
 * Report positions the way the submitter can act on them: as a 1-based index
 * into the sequence they pasted, not into the reverse complement we may have
 * compared, and not into a copy with the gap characters taken out.
 *
 * `gapMap[i]` is the position, in the cleaned sequence including gaps, of the
 * i-th ungapped base.
 */
function submissionPosition(orientedIndex, orientedLength, reversed, gapMap) {
  const ungappedIndex = reversed ? orientedLength - 1 - orientedIndex : orientedIndex;
  return gapMap[ungappedIndex] + 1;
}

/** Positions of the non-gap characters of a cleaned sequence, in order. */
function ungappedToCleaned(sequence) {
  const map = [];
  for (let i = 0; i < sequence.length; i++) if (sequence[i] !== "-") map.push(i);
  return map;
}

/** "positions 231 and 402" / "positions 12, 231 and 402" -- no Oxford comma. */
function joinList(items) {
  if (items.length <= 1) return items.join("");
  return items.slice(0, -1).join(", ") + " and " + items[items.length - 1];
}

/**
 * How much of the barcode region the submission covers, and how much of it sits
 * outside.
 */
function coverageCheck(index, placement) {
  const { inside, before, after, firstColumn, lastColumn } = placement;
  const outside = before + after;
  const window = index.windowLength;

  const span =
    inside === window
      ? `the whole ${window} bp barcode region`
      : `${inside} bp of the ${window} bp barcode region (positions ` +
        `${firstColumn}–${lastColumn})`;

  if (outside === 0) {
    return { state: "pass", label: "Coverage", text: `Covers ${span}.` };
  }

  /* Sequenced past the region. The comparison is still sound -- only the
     overlap is ever scored -- but the submitter has no way of knowing that from
     a length alone, and the flanking bases are neither checked nor named from. */
  const where = [];
  if (before) where.push(`${before} bp before it`);
  if (after) where.push(`${after} bp after it`);
  return {
    state: "warn",
    label: "Coverage",
    text:
      `Covers ${span}, plus ${joinList(where)}. Those ${outside} bp lie outside the ` +
      "barcode region: they were not compared, and they play no part in whether this " +
      "is a new lineage. Trim to the barcode region before depositing."
  };
}

/**
 * Read the submission in frame and report stop codons.
 *
 * Only codons that lie wholly inside the window and are wholly concrete are
 * read: an ambiguity code could stand for a stop or for something else, and
 * guessing either way would be inventing a finding.
 */
function frameCheck(placement, position) {
  const { oriented, columns } = placement;

  // Column -> index in the oriented sequence, so a codon can be assembled from
  // its three columns rather than from three adjacent bases (which is only the
  // same thing when nothing is missing).
  const byColumn = new Map();
  for (let i = 0; i < columns.length; i++) byColumn.set(columns[i], i);

  let codonsRead = 0;
  const stops = [];
  for (const [column, index] of byColumn) {
    if ((column - 1) % 3 !== 0) continue; // not a first codon position
    const second = byColumn.get(column + 1);
    const third = byColumn.get(column + 2);
    if (second === undefined || third === undefined) continue;
    const codon = oriented[index] + oriented[second] + oriented[third];
    if (/[^ACGT]/.test(codon)) continue;
    codonsRead++;
    if (!STOP_CODONS.has(codon)) continue;
    // Report the span in the submitter's own coordinates, low end first --
    // a reverse-complemented read runs the other way.
    const ends = [position(index), position(third)].sort((a, b) => a - b);
    stops.push({ codon, column, from: ends[0], to: ends[1] });
  }

  if (codonsRead === 0) return null; // nothing readable in frame; say nothing

  if (stops.length === 0) {
    return {
      state: "pass",
      label: "Frame",
      text: `Reads through all ${codonsRead} codons without a stop.`
    };
  }

  stops.sort((a, b) => a.from - b.from);
  const listed = stops
    .slice(0, 4)
    .map((s) => `${s.codon} at positions ${s.from}–${s.to}`);
  const rest = stops.length - listed.length;

  return {
    state: "fail",
    label: "Frame",
    text:
      `Stop codon${stops.length === 1 ? "" : "s"} in the reading frame: ` +
      joinList(listed) + (rest > 0 ? `, and ${rest} more` : "") +
      ". Cytochrome b should read through the whole barcode region, so this " +
      "usually means a mis-called base, an insertion or deletion that has shifted " +
      "the frame, or sequence pasted in from outside the region. " +
      "(TAA and TAG are the only stops in these parasites; TGA codes tryptophan.)"
  };
}

/**
 * Compare each base against everything MalAvi has ever seen at that column.
 *
 * A base at a position where no lineage carries it is the single most useful
 * thing this profile can say: 97% of the lineages in the release carry no such
 * base at all when held out of their own reference, so it is rare enough to be
 * worth a second look at the trace, and common enough in bad reads to catch
 * real errors.
 */
function sitesCheck(index, placement, position) {
  const { oriented, columns } = placement;
  const unseen = [];
  let scarce = 0;
  let judged = 0;

  for (let i = 0; i < oriented.length; i++) {
    const base = oriented[i];
    const count = siteBaseCount(index, columns[i], base);
    if (count < 0) continue; // outside the window, ambiguous, or a thin column
    judged++;
    if (count === 0) unseen.push({ base, position: position(i), column: columns[i] });
    else if (count <= 2) scarce++;
  }

  if (judged === 0) return null;

  const scarceNote = scarce
    ? ` ${scarce} other base${scarce === 1 ? " sits" : "s sit"} at a position where ` +
      `${scarce === 1 ? "it is" : "they are"} carried by no more than two lineages.`
    : "";

  if (unseen.length === 0) {
    return {
      state: "pass",
      label: "Sites",
      text: `Every base is one MalAvi has seen at that position before.${scarceNote}`
    };
  }

  unseen.sort((a, b) => a.position - b.position);
  const listed = unseen.slice(0, 6).map((u) => `${u.base} at position ${u.position}`);
  const rest = unseen.length - listed.length;

  return {
    state: "warn",
    label: "Sites",
    text:
      `${joinList(listed)}${rest > 0 ? `, and ${rest} more,` : ""} ` +
      `${unseen.length === 1 && rest === 0 ? "is a base" : "are bases"} found at ` +
      "that position in no other MalAvi lineage. That can be perfectly real, but it is " +
      "also what a mis-called base looks like — please check the trace at " +
      `${unseen.length === 1 ? "that position" : "those positions"} before naming it.` +
      scarceNote
  };
}

/* ---- NAMING A NEW LINEAGE ---------------------------------------------------
   From the submission guide, in Staffan Bensch's words: a new lineage gets "a
   5-6 letter acronym based on the scientific name of the first encountered host
   species followed by a two-digit number. The name must be unique for the
   database."

   The acronym is genus letters plus epithet letters, and BOTH widths are in
   real use -- Turdus migratorius has 24 lineages named TUMIG and 17 named
   TURMIG. So the checker does not pick one. It works out both candidates, says
   which numbers each has already used, and leaves the choice where it belongs:
   with the curator, who is confirming the name anyway. */

/** Two-digit zero padding, widening past 99 rather than truncating. */
function padNumber(n) {
  return String(n).padStart(2, "0");
}

/**
 * Which lineage names already use this acronym, and what the next number is.
 *
 * `claims` is the list of names claimed by submissions that have been received
 * but are not yet in a release, each with the date it was claimed. A claimed
 * number is not free: the submission guide's reason for sending names before
 * publication is that they are held for you, and priority runs by the date the
 * submission arrived. So the next proposal steps over the claims as well as over
 * the release, and the claims are reported rather than silently skipped -- being
 * told "TUMIG25 is spoken for, you would be TUMIG26" is the useful answer.
 */
function acronymUsage(index, acronym, claims) {
  const pattern = new RegExp("^" + acronym + "(\\d+)$");

  const numbers = [];
  for (const name of index.allNames) {
    const hit = pattern.exec(name);
    if (hit) numbers.push(Number(hit[1]));
  }

  const claimed = [];
  for (const claim of claims) {
    const hit = pattern.exec(claim.name);
    if (hit) claimed.push({ name: claim.name, claimed: claim.claimed, number: Number(hit[1]) });
  }
  claimed.sort((a, b) => a.number - b.number);

  const highestInRelease = numbers.length ? Math.max(...numbers) : 0;
  const highestClaimed = claimed.length ? claimed[claimed.length - 1].number : 0;
  const next = Math.max(highestInRelease, highestClaimed) + 1;

  return {
    acronym,
    taken: numbers.length,
    highest: numbers.length ? acronym + padNumber(highestInRelease) : null,
    claims: claimed.map((c) => ({ name: c.name, claimed: c.claimed })),
    proposal: acronym + padNumber(next)
  };
}

/**
 * Turn a host species name into the lineage names that are still free.
 *
 * Pure and offline, like everything else here: it reads the lineage names in
 * the loaded release, plus whatever names pending submissions have claimed, and
 * nothing else. It proposes; it does not assign.
 *
 * `reservations` is the reserved_names.json payload, or its `names` array, or
 * nothing at all. Nothing at all is a supported answer, not a failure: the feed
 * is optional, and a page that could not load it still gives the right answer
 * about the release -- it just cannot speak for submissions in the queue.
 *
 * `taxonomy` is the bird_names.json payload, and is optional in the same way.
 * When supplied, the result carries a `taxonomy` verdict on the host name (see
 * checkHostName). It never changes what is proposed: an unrecognized name still
 * gets its acronym, because taxonomy moves and a dated checklist is not the
 * authority on what a submitter found.
 */
export function suggestLineageName(index, hostName, reservations, taxonomy) {
  const claims = Array.isArray(reservations)
    ? reservations
    : (reservations && reservations.names) || [];

  /* Keep only things that look like words. This drops "sp.", "cf.", authorities
     and any punctuation, so "Turdus sp." is correctly treated as a genus with no
     epithet rather than as a binomial. */
  const words = String(hostName == null ? "" : hostName)
    .split(/[^A-Za-z]+/)
    .filter((word) => word.length >= 3);

  if (words.length < 2) {
    return {
      ok: false,
      message:
        "Give the host's scientific name — genus and species, as in " +
        "Turdus migratorius. The acronym is built from both."
    };
  }

  const genus = words[0].toUpperCase();
  const epithet = words[1].toUpperCase();

  /* Both widths that MalAvi uses, longer first: 3,336 of the 5,367 names in this
     release use a six-letter acronym and 1,343 a five-letter one. A duplicate is
     dropped, which is what happens when the genus is only two letters long. */
  const acronyms = [...new Set([
    genus.slice(0, 3) + epithet.slice(0, 3),
    genus.slice(0, 2) + epithet.slice(0, 3)
  ])];

  const options = acronyms.map((acronym) => acronymUsage(index, acronym, claims));
  // An acronym already in use for this host is the one to follow; show it first.
  options.sort((a, b) => b.taken - a.taken);

  return {
    ok: true,
    host: words[0].charAt(0).toUpperCase() + words[0].slice(1).toLowerCase() +
      " " + epithet.toLowerCase(),
    options,
    inUse: options.some((o) => o.taken > 0),
    // Whether any pending claim bore on this host at all, so the page can say
    // when it has checked the queue and found nothing rather than staying silent.
    claimsChecked: claims.length > 0,
    claimed: options.some((o) => o.claims.length > 0),
    // What the avian checklist makes of the name, when one was supplied.
    taxonomy: checkHostName(taxonomy, words[0], words[1])
  };
}

/**
 * Is this a real bird, according to the eBird/Clements checklist?
 *
 * Returns a verdict, never a veto. The checker's job here is to catch a typo or
 * an invented name before it becomes a lineage name in a paper and in GenBank --
 * not to decide what exists. Taxonomy moves, this checklist is a dated snapshot,
 * and a submitter naming a bird it has not caught up with is right more often
 * than the file is. So an unrecognized name is reported and the acronym is
 * offered anyway.
 *
 * Deliberately checked against the WHOLE checklist rather than MalAvi's own host
 * list. The interesting case for a new lineage is a parasite sequenced from a
 * host nobody has screened before: 11,167 species are in the checklist and only
 * about 2,300 have ever appeared in MalAvi, so validating against MalAvi's hosts
 * would reject the great majority of legitimate new host species.
 *
 * `checklist` is the bird_names.json payload, or nothing. Nothing is a supported
 * answer -- the page still works, it just cannot speak about the name.
 *
 * The returned `status` is one of:
 *   "unchecked"  no checklist was loaded
 *   "accepted"   a current eBird species
 *   "synonym"    a valid older name; `current` gives the name it resolves to
 *   "genus-only" the genus is real but the epithet is not one of its species
 *   "unknown"    neither the genus nor the name is in the checklist
 */
export function checkHostName(checklist, rawGenus, rawEpithet) {
  if (!checklist || !checklist.accepted) return { status: "unchecked" };

  /* The checklist stores names as they are written: capitalised genus, lower
     case epithet. Submitters type all sorts of things, so compare on a
     normalised form rather than demanding they get the casing right. */
  const genus = rawGenus.charAt(0).toUpperCase() + rawGenus.slice(1).toLowerCase();
  const epithet = rawEpithet.toLowerCase();
  const binomial = genus + " " + epithet;

  const species = checklist.accepted[genus];
  if (species && species.indexOf(epithet) !== -1) {
    return {
      status: "accepted",
      name: binomial,
      family: (checklist.families && checklist.families[genus]) || null
    };
  }

  const synonyms = checklist.synonyms && checklist.synonyms[genus];
  const current = synonyms && synonyms[epithet];
  if (current) {
    /* An older name. Both forms are reported and neither is chosen: the acronym
       should follow the name the submitter will publish, and which that is
       depends on what their paper says, not on what this file prefers. */
    return { status: "synonym", name: binomial, current };
  }

  if (species) {
    return {
      status: "genus-only",
      name: binomial,
      genus,
      family: (checklist.families && checklist.families[genus]) || null
    };
  }
  return { status: "unknown", name: binomial };
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
    /* Distinguish "we found nothing" from "we could not look". A query with an
       ambiguity code in every 16 bp window leaves the exact-word index with
       nothing to seed from, and the old copy below then told the submitter their
       sequence was probably very divergent or not cytochrome b at all -- a
       confident claim the checker had no evidence for, and one that was
       measurably false: these sequences were named lineages
       (METHODS_draft.md 5H.2). Expanded seeds handle the ordinary degraded read;
       when even that is not enough, say so plainly.

       Only for input that is actually nucleotide sequence: something with
       non-nucleotide characters in it (protein, pasted prose) also has no usable
       seeds, but "too ambiguous" is the wrong explanation for it -- the Content
       check above has already said the right one. */
    const seeds = countUsableSeeds(sequence);
    if (content.invalid === 0 && seeds.sampled > 0 && seeds.usable === 0) {
      checks.push({
        state: "fail",
        label: "Identity",
        text:
          "Not checked — every part of this sequence carries too many ambiguous " +
          "positions for the quick comparison. This is not evidence either way."
      });
      return {
        verdict: "unknown",
        title: "Too ambiguous to check here",
        message:
          "The ambiguous positions are spread so evenly that no stretch of this " +
          "sequence is definite enough for the quick check to use. That says nothing " +
          "about whether it is new — it may well be a named lineage. Re-check the " +
          "trace if you can, run it through BLAST, and send it to us either way.",
        content,
        match: null,
        checks
      };
    }

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

  /* With a match in hand the submission can be placed in the alignment window,
     which is what the remaining checks all read from. */
  const placement = placeInWindow(index, sequence.replace(/-/g, ""), best);
  const gapMap = ungappedToCleaned(sequence);
  const position = (i) =>
    submissionPosition(i, placement.oriented.length, best.orientation === "reverse", gapMap);

  checks.push(coverageCheck(index, placement));

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

  /* Reading frame and per-position novelty, the two checks malaviR's
     lineage_qc() runs on a candidate lineage. They are worth running here for
     the same reason they exist there: a sequence that differs from everything
     named is either a new lineage or a bad read, and these are what tell those
     two apart.

     They are run only when the sequence DOES differ from its closest match. An
     exact match to a named lineage is already in the database with whatever
     properties it has -- 18 lineages in this release carry a stop codon -- and
     reporting the database's own quirks back to someone who has simply
     re-sequenced a known lineage would be noise, not a finding. */
  const sequenceQuality = [
    frameCheck(placement, position),
    sitesCheck(index, placement, position)
  ].filter(Boolean);

  if (content.unambiguous < minLength) {
    checks.push({
      state: "fail",
      label: "Identity",
      text:
        `Closest is ${nameList} at ${best.differences} difference` +
        `${best.differences === 1 ? "" : "s"} over ${best.compared} bp, but the sequence ` +
        "is too short to name a new lineage from."
    });
    for (const quality of sequenceQuality) checks.push(quality);
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
  for (const quality of sequenceQuality) checks.push(quality);
  checks.push({
    state: "warn",
    label: "Naming",
    text:
      "Needs a name. Enter the host species below and we will show you which numbers " +
      "its acronym has already used."
  });

  return {
    verdict: "new",
    /* Tells the page to offer the naming box. The suggestion itself is made by
       suggestLineageName once the submitter has told us the host, which is the
       one thing a sequence cannot tell us. */
    naming: true,
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
