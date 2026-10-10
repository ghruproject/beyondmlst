# Independent assembly tree

`chronoclade tree <selection.json> --out tree/ --threads 4` consumes a validated
cgMLST or ESM2 selection. It checks every source checksum and builds precisely the
ordered selected IDs. It never retrieves other context or chooses replacements.
Linked local assemblies are validated and copied. Missing assemblies are acquired
by exact Pathogenwatch `source_genome_id` using the configured API key; credentials
are never saved. Acquisition failure leaves an explicit `assembly_audit.json` in
the failed job and preserves any previous complete result.

SKA reference mapping, IQ-TREE GTR+G and ClonalFrameML reuse the existing scientific
engines. The reference is the selected assembly with greatest N50, with a stable
ID tie-break. Alignment/coherence failures fail the exact selection rather than
silently dropping genomes. The saved result records raw and corrected trees,
filtered alignment and callable-site counts, executable hashes/package versions,
reference, importation/filtering evidence, SNP pair distances and the original
selection/dataset hashes. Undated samples remain in the phylogeny.

`tree.json` is published last; artifacts live in fingerprinted job directories.
An unchanged complete run validates and resumes without rebuilding. Parameter,
local assembly or executable changes invalidate that job. `--force` creates a
fresh job without deleting the previous complete result. A failed replacement
never publishes a replacement `tree.json`.

The offline report retains the full-report style, corrected biological tree,
country/host metadata including missing categories, exact location-history network
and temporal results. SNP neighbours have SNP/site evidence and do not inherit
cgMLST allele-distance labels. The corrected tree uses substitutions/site.

## Temporal gate

The default runs 20 date permutations with the least-squares root refitted for the
observed clock and every permutation. Maximal monophyletic clades sharing the
same exact collection date string are the exchangeability units, following
[Murray et al. 2016](https://doi.org/10.1111/2041-210X.12466). Tree topology is inferred
without dates; same-date clades are then identified on that frozen topology. One
date is permuted per cluster and expanded to its members. Year/month/day strings
are retained, permitting TreeTime's partial-date interpretation. Cluster date
multiplicities are preserved; individual-tip date multiplicities can change with
unequal cluster sizes. Tip weights remain unchanged in clock fitting.

The saved assessment includes every eligible ID, excluded undated IDs, actual
cluster memberships, seed, requested/successful refits, root policy and a
descriptive pairwise genetic-distance/date-gap correlation. That correlation is
not an independent-pair significance test. Positive rate and the existing
root-to-tip permutation threshold gate are required; incomplete permutations,
fewer than three distinct cluster dates or disabled assessment remain unsupported.
A passing clustered screen does not establish model adequacy or remove every
sampling bias. This is not the full-model CR2 uncertainty-overlap test.

`--no-assess-temporal` skips clock/permutation tools. Fewer than three usable
distinct dates writes a complete tree product with an unsupported assessment.
Tree never creates a time-scaled phylogeny; it saves `assessment.json` for time.
