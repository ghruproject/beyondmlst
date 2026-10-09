# Native query-typing validation

On 9 October 2026, ChronoClade installed the pinned Pathogenwatch cgMLST caller,
plincer and hclink through `setup-typing` on an Apple Silicon laptop. The Python
tools used isolated uv environments, and the cgMLST caller used the installed
Node 22 runtime and Pixi BLAST. Both organisms were then typed through the actual
`type-queries` command, using controlled three-locus assemblies and reference
fixtures. Both returned `1_1_1`; plincer and hclink returned the expected fixture
assignments. Tool commits, input/reference hashes and outputs are in
[results.json](results.json).

These small references exercise the real upstream programs and the integration.
They are not production Klebsiella or E. coli schemes and do not validate biological
assignments. Production reference databases were not built because the required
Pasteur and EnteroBase credentials were unavailable. The default setup configuration
correctly reports that it is not ready for typing.

The same validation also applied refinement to the real frozen public ST147
catalogue. Of 5,372 eligible same-ST samples, 1,724 had compatible priority evidence
for at least one of the three public focal stand-ins. An eight-genome budget selected
six prioritised samples and two balanced background samples. All three queries had
comparable evidence. This was a selection replay; it did not rerun the old pilot tree
with those new selections. Unversioned cgLIN comparisons were scoped to the exact
frozen export, rather than treated as a known nomenclature release.

Automated end-to-end tests additionally cover the native-provider path, imported
query/public assignments, E. coli context discovery, partial/ambiguous calls,
reference version mismatches, fair per-query prioritisation and the background
budget. See [the usage guide](../../docs/query-typing.md) for installation and database
preparation commands.
