# Profile-first workflow validation

Validated locally on 9 October 2026.

- Live Pathogenwatch collection `jX5cwsoyJ1KDquMqssUzAD`: 23 query genomes, all profiles retrieved through the collection's advertised cgMLST job, five locus bootstrap replicates, 253 comparable pairs, no assembly downloads. This is the supplied *S. aureus* collection, not the Klebsiella ST147 dataset.
- Exact source-ID accession input for Klebsiella ST147: one query plus eight bounded public context profiles, all nine profiles available, 36 comparable pairs. The frozen discovery pool contains 7,807 records and 5,721 deduplicated sample units; 5,720 remain eligible after excluding the query. No context assemblies were downloaded by fast mode. The observed export locus union is not claimed to establish canonical scheme coverage.
- The seven-genome ST147 assembly pilot exercised the public full-stage CLI, producing a recombination-adjusted tree and corrected report without temporal testing. An isolated native finish test retained one undated genome and ran three successful date permutations. Temporal signal was unsupported and no dated tree was reported. Three permutations verify wiring, not biological evidence.
- Reader-facing profile reports were inspected in the browser. Earlier stages preserve independent figures, tables and bundles. Future collection dates are preserved as source metadata and excluded from temporal calculations.

Production typing of previously untyped assemblies still requires the appropriate licensed reference datasets. The Pasteur API key request remains pending; upstream native tool fixtures and already typed public genomes do not substitute for that production validation. SLURM execution requires access to the user's cluster and has not been exercised here.

Reproduce the query-only collection report:

```sh
pixi run chronoclade run --collection jX5cwsoyJ1KDquMqssUzAD \
  --mode fast --profile-limit 0 --profile-bootstraps 5 --output collection-results
```

Use `--mode full` to add selected assemblies and a corrected tree, or `--mode finish` to add the temporal gate. The defaults retain every query and select up to 50 context assemblies per species/ST analysis lineage.
