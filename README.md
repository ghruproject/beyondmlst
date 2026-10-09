# ChronoClade

ChronoClade tests whether sampling dates contain enough information to estimate a
molecular clock in a bacterial lineage. It builds a whole-genome phylogeny,
accounts for recombination with ClonalFrameML, examines the root-to-tip
relationship, and runs a date-randomisation test. TreeTime estimates a dated
phylogeny only when the temporal-signal test passes.

Two analysis modes are available. `full` uses ClonalFrameML and can use either
the quick root-to-tip permutation screen or complete TreeTime refits. `fast`
uses PhiPack Profile to screen for recombination, then builds an uncorrected fast
tree for root-to-tip permutations. It is triage: PhiPack does not localise or
remove recombinant tracts here, and fast mode does not produce a dated tree or
an epidemiological interpretation.

The workflow is intended for longitudinal surveillance within a species and
lineage, such as an MLST sequence type. Public context genomes can be added to
help distinguish a sampled local lineage from separate introductions. The final
HTML report starts with countries, closest analysed relatives and a genetic tree,
then explains interpretation and dating. Detailed methods and the full public
country atlas are available separately.

ChronoClade does not infer direct transmission. A dated tree estimates ancestral
times under a molecular-clock model; it does not identify who infected whom.

## Installation

[Pixi](https://pixi.sh/) installs ChronoClade and its native dependencies on
macOS or Linux.

```bash
git clone https://github.com/ghruproject/chronoclade.git
cd chronoclade
pixi install
pixi run chronoclade preflight
```

Windows users can run the Pixi environment under WSL2.

## Basic use

ChronoClade reads a CSV containing one row per assembly. Samples are analysed
separately by `species` and `lineage`.

```csv
sample_id,assembly,collection_date,location,species,lineage,origin,is_reference,patient_id
KPN001,assemblies/KPN001.fasta,2023-03,KIMS,Klebsiella_pneumoniae,ST15,local,true,P001
KPN002,assemblies/KPN002.fasta,2024-01-17,KIMS,Klebsiella_pneumoniae,ST15,local,false,P002
KPN_REF,context/KPN_REF.fasta,2021,Philippines,Klebsiella_pneumoniae,ST15,context,false,
```

Check the metadata and planned analyses before running the workflow:

```bash
pixi run chronoclade validate metadata.csv
pixi run chronoclade run metadata.csv --dry-run
```

Run the analysis:

```bash
pixi run chronoclade run metadata.csv \
  --output chronoclade_results \
  --threads 8 \
  --lineage-jobs 2 \
  --randomisation-jobs 4 \
  --date-randomisations 100
```

For a quick screen:

```bash
pixi run chronoclade run metadata.csv --mode fast --date-randomisations 100
```

For the more demanding date-randomisation test, refitting the full TreeTime
model after every permutation:

```bash
pixi run chronoclade run metadata.csv \
  --mode full \
  --date-randomisation-method full-tree \
  --date-randomisations 100
```

After SKA mapping, ChronoClade checks whether any genome is an extreme
raw-distance outlier within its declared lineage. A failed check stops before
IQ-TREE and ClonalFrameML and records the per-genome evidence in
`lineage_coherence.tsv`; verify the accession and lineage assignment before
rerunning. A second check on the recombination-filtered alignment writes
`clonal_lineage_coherence.tsv` before temporal analysis.

Open `chronoclade_results/index.html` when the run completes. Each lineage has
its own report and a ZIP archive containing the figures, tables, trees and
machine-readable results used in that report.

## Public context genomes

Pathogenwatch is the primary public-context provider, initially for
*Klebsiella pneumoniae*. ChronoClade freezes the full public same-ST catalogue,
imports comparable cgLIN assignments, then downloads a bounded country/year
balanced candidate pool for SKA screening. Country figures use QC-passing,
accession-deduplicated public records, including undated records, before
selection. cgLIN depths are prefix levels, not SNP cutoffs or transmission groups.

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme klebsiella --st 147 --output context/ST147 \
  --candidate-pool 40 --max-context 12 --threads 4 --dry-run

# Inspect the frozen catalogue and audit, then rerun without --dry-run.
pixi run chronoclade prepare-context focal_metadata.csv \
  --scheme klebsiella --st 147 --output context/ST147 \
  --candidate-pool 40 --max-context 12 --threads 4

pixi run chronoclade run context/ST147/combined_metadata.csv \
  --context-manifest context/ST147/context_manifest.tsv \
  --output chronoclade_results --threads 4
```

Downloads use `PATHOGENWATCH_API_KEY` or a user-only configuration file; see the
[provider documentation](docs/pathogenwatch.md). A validated `--cglin-export`
and `--focal-crosswalk` can supply assignments without uploading focal genomes.
Missing assignments remain explicit. See the [country figures](docs/context-geography.md)
and [laptop/SLURM pilot](docs/pathogenwatch-pilot.md).

AllTheBacteria is retained as an explicit legacy/deferred route:

```bash
pixi run chronoclade prepare-context focal_metadata.csv \
  --context-source atb --scheme ecoli_achtman_4 --st 131 \
  --source auto --output context/ST131
```

`--source` continues to select the ATB transport (`auto`, `aws`, `osf`).
Unsupported Pathogenwatch organisms do not silently fall back to ATB.

## Worked example

The repository contains a real-data analysis of 72 *Staphylococcus aureus*
ST239 genomes sampled between 1980 and 2012. The estimated clock rate was
1.693 x 10^-6 substitutions per site per year, the inferred root date was 1946,
and the date-randomisation p-value was 0.0099. These estimates agree closely
with the published analysis of the collection.

The [worked example](https://ghruproject.github.io/chronoclade/worked-example/)
shows the commands, root-to-tip plot, randomisation result, dated tree with node
intervals, clonal distances and the resulting public-health interpretation.

The repository also includes an
[accession-defined ST131-H30 validation](validation/komori2024_st131/README.md)
against the 96 C0/C1 chromosomes from Komori et al. (2024). It records the
identity audit, both lineage-coherence checks and the distinction between
detectable temporal signal and an imprecise TreeTime root date.

## Documentation

The [ChronoClade documentation](https://ghruproject.github.io/chronoclade/)
covers input preparation, temporal-signal testing, report interpretation and
the command-line interface.

To build it locally:

```bash
pixi run docs
```

## Methods and limitations

ChronoClade uses SKA2 for reference-ordered whole-genome alignment, IQ-TREE for
the starting maximum-likelihood phylogeny, ClonalFrameML to account for
recombination, and TreeTime for clock analysis and time scaling. The observed
root-to-tip fit is compared with fits obtained after permuting collection dates.
The default screen requires a positive clock rate and an empirical R²
permutation p-value of 0.05 or less. The optional full-tree test reruns the
complete TreeTime fit and applies the stricter CR2 rate-interval rule.

Root-to-tip regression is a diagnostic rather than a formal test. Population
structure, biased sampling, date uncertainty and residual recombination can all
affect molecular-clock estimates. ChronoClade retains the source tables and
figures so these assumptions can be reviewed.

See the [methods page](https://ghruproject.github.io/chronoclade/methods/) for
the full workflow and references.

## Development

```bash
pixi run lint
pixi run test
pixi run docs
```

ChronoClade is under active development. Validate results against epidemiology
and the underlying phylogenetic outputs before using them to guide an
investigation.
