# ESM2 temporal validation fixture

Five simulated samples, four loci and six unique proteins. Sample dates, countries,
roles and allele identifiers are synthetic and do not describe a real outbreak.
Reference proteins are the first four frozen public Pasteur allele-1 proteins in
`../fixture-manifest.json`. Two labelled amino-acid substitutions were introduced.
SIM_02 has a different DNA allele identifier but the same protein as SIM_01;
SIM_05 has no collection date.

Vectors were generated with real ESM2 8M inference on Apple MPS. The saved manifest
records the checkpoint hash, pooling, device and vector checksum. Tests reuse these
small frozen vectors without importing Torch or downloading weights. This verifies
the data/report pipeline, not biological temporal signal or a molecular clock.
