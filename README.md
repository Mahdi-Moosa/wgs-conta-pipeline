# WGS contamination pipeline

The WGS contamination pipeline detects and quantifies unintended human-reference-aligned
sequence in bacterial isolate short-read whole-genome sequencing libraries. It
combines public-archive metadata and read acquisition with conservative human-fragment
calling, chromosome summaries and exact-assembly adjudication.

## Scope

Software version **0.1.0rc1** implements conservative human-sequence detection
and quantification using reference alignment, mapping-quality, alignment-coverage
and sequence-identity filters, with layout-specific fragment counting and
coordinate deduplication. See [method fidelity](docs/methods_consistency.md) for
detailed scientific definitions and implementation conventions.

Supported combinations:

| Archive platform | Layout | Counted unit |
|---|---|---|
| ILLUMINA | PAIRED | coordinate-deduplicated fragment |
| BGISEQ / DNBSEQ | PAIRED | coordinate-deduplicated fragment |
| ION_TORRENT | SINGLE | coordinate-deduplicated read |

The pipeline was developed and validated primarily for bacterial isolate
short-read WGS datasets. Its core methods for human-reference alignment and
contamination quantification may also apply to other sample types, including
metagenomic datasets, but these applications have not been systematically
validated. Archive eligibility rules and bacterial assembly-based adjudication
are specific to the bacterial isolate workflow and may require adaptation for
other applications. The current CLI supports only the platform/layout
combinations listed above; potential applicability does not imply end-to-end
support for arbitrary sample types. Long-read sequencing and additional
platform/layout combinations have not been validated and are unsupported.
Archive labels alone do not establish biological eligibility.

## Capabilities

The pipeline supports the following analyses and workflows:

- **Archive exploration:** Retrieve sequencing metadata, evaluate bacterial-isolate eligibility, and select independent cohorts using metadata-based sampling.
- **Read acquisition:** Download and validate sequencing reads from public archives.
- **Human-sequence detection:** Identify high-confidence human-aligned reads and fragments using stringent alignment-quality, coverage and sequence-identity criteria.
- **Contamination quantification:** Measure human-sequence burden and chromosome-level representation, with Kraken2-based corroboration.
- **Assembly-based assessment:** Evaluate human-aligned reads against bacterial assemblies linked to the same BioSample.
- **Batch processing:** Analyze multiple sequencing libraries with configurable parallel execution, checkpoints, failure tracking and resume support.
- **Results export:** Generate run-level results, chromosome summaries and aggregate outputs for downstream analysis.

## Installation

Linux is required for the production runner (file locking and subprocess tools).
The core Python package uses only the standard library and requires Python >=3.10.
The documented Conda environment was tested with Python 3.14.10. See
[validation](docs/validation.md) for installation tests, integration coverage and
known limits.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install .
wgs-conta --help
wgs-conta run --help
python -m unittest discover -s tests -v
python examples/synthetic_smoke.py
```

Tests use artificial records and include real samtools filtering and a
minimap2/samtools assembly challenge when those tools are available. They do not
reconstruct historical manuscript cohorts.

For external scientific software, [environment.yml](environment.yml) records the
validated tool versions:

```bash
conda env create -f environment.yml
conda activate wgs-conta-release
python -m pip install .
```

Full runs require minimap2, samtools, seqkit, Kraken2, curl, SRA Toolkit, EDirect
and NCBI Datasets. AWS CLI is needed for AWS acquisition; HTTPS ENA acquisition
is the example default. Aspera is optional and not bundled. Tool/reference licenses
are separate. See [dependencies](docs/dependencies.md).

## Reference preparation

The pipeline uses the **T2T-CHM13v2.0 masked-Y + rCRS** human reference genome:
the Y pseudoautosomal regions are masked, and the mitochondrial sequence is rCRS.
Reference sequences and minimap2 indexes are **not distributed with the software**.

Obtain `chm13v2.0_maskedY_rCRS.fa.gz` separately from the
[official T2T resource page](https://genomeinformatics.github.io/CHM13v2/), or use
the [reference FASTA download](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz).
Download and decompress the reference, then prepare a minimap2 short-read index
(`-x sr`). The pipeline checks the FASTA against its validated SHA-256 and the
index against the configured SHA-256.

A compatible Kraken2 database must also be provided separately, with a verified
version/checksum manifest. See [reference preparation](docs/references.md) for
download, checksum, indexing and configuration commands.

Copy `configs/pipeline.example.ini` to `configs/local.ini` and set paths and
reference identity. Paths are relative to the INI file, not to a research tree.
Keep reference resources outside all workspace/cleanup roots. By default all
runtime files live in `workspace/`, which is excluded from Git. Keep source,
reference and working directories separate. Do not put private data in this repo.

```bash
wgs-conta validate --config configs/local.ini --paths-only
# After preparing the external resources:
wgs-conta validate --config configs/local.ini
```

## Metadata and independent cohort selection

`configs/archive.example.json` is an illustrative, user-editable time window,
not a historical manuscript sampling frame. Declare the window before outcomes
are measured. Metadata capture queries bacterial WGS/GENOMIC records for the
supported platform/layout combinations, uses half-open timestamp partitions,
checks counts, verifies partition integrity and resumes completed partitions.
Capture contacts ENA and may take substantial time; it is never run on import.

```bash
wgs-conta-archive --config configs/archive.example.json --output workspace/archive
# On a selected raw partition, decompress TSV if needed:
wgs-conta classify --input workspace/ena_records.tsv --output workspace/classified.tsv
```

Classification emits explicit exclusions and direct-isolate evidence. Taxonomy
comes from the archive bacterial query; imported arbitrary TSVs require the same
upstream taxonomic restriction. Clinical specimen words in isolation_source alone
are not sufficient to establish direct-specimen sequencing. Review contradictory
or incomplete records before freezing the selection. Optional `--taxon-deny-list`
accepts an explicit JSON list of names/IDs for prospectively declared exclusions.
These later deterministic metadata rules support independent studies; they are
not represented as the exact historical E selection producer.

The included selector preserves the reviewed HMAC/Hamilton implementation. It
requires a **reviewed canonical frame**, including `country_compare_status=MATCH`
and external ISO3/M49 mappings. ENA-to-NCBI country agreement is not inferred by
this release: obtain and validate it separately, and do not fabricate MATCH.
`classify` can populate country fields with `--m49-map` plus `--iso-json`, but does
not itself perform that independent cross-check. Schemas are in
[inputs and outputs](docs/schemas.md). The selection API may also be used directly
when designing an independent sampling protocol.

```bash
wgs-conta new-seed --output workspace/private/selection.key
wgs-conta select --input workspace/reviewed_frame.tsv --m49-map references/m49.tsv \
  --seed-file workspace/private/selection.key --projects 20 \
  --biosamples-per-project 3 --runs-per-biosample 1 --cohort-label INDEPENDENT \
  --output workspace/private/registry.tsv
```

This creates a fresh user-controlled seed; no original manuscript seed is
included. Roster/seed files are private local artifacts. Independent samples can
overlap unless an explicit, prospectively declared exclusion is applied. The
software does not reconstruct or promise the exact historical cohorts or results.

## Execute, inspect, resume

Alternatively supply your independently reviewed registry directly; column
`Accession` is required. Never paste the historical private manuscript roster.

```bash
wgs-conta run --config configs/local.ini --accessions workspace/private/registry.tsv --dry-run
wgs-conta run --config configs/local.ini --accessions workspace/private/registry.tsv --source all
# Same command resumes scientifically incomplete work; complete compatible results skip.
wgs-conta summarize --input workspace/master/master_results.tsv --output workspace/aggregate.json
```

Dry-run makes no network/tool calls but may create the configured master directory.
SRA cache lookup is accession-driven, bounded and nonrecursive. Acquired SRA is
validated; biological PASS/CRITERIA reads are retained and biological REJECT
reads/spots excluded with accounting checks. Technical/index rejection alone does
not discard a valid biological pair. ENA files are MD5-verified with bounded
retries. Equal paired FASTQ record counts are required. The original input
normalization and QNAME conventions are retained; there is no trimming or read
correction step. User-supplied cached SRA is not deleted. The example retains new
SRA too; work intermediates are cleaned only after a successful permanent commit.

The pipeline keeps a private runtime SQLite queue and permanent TSV master.
Completion requires a result plus all 25 chromosome rows with matching summed
strict units. Interrupted commits can be reconciled. Failures remain explicitly
failed/incomplete. A release identity guard rejects silent reuse after science
code, reference identity or tracked tool-version changes; use a new workspace.
Changes in an external database must be reflected in its verified identity file.
Database content is not recursively scanned or rehashed on every run.

## Endpoints and interpretation

Primary alignments exclude unmapped, secondary and supplementary records
(`0x904`). Strict reads require MAPQ >=30, query coverage >=80% and identity
>=98%. Paired libraries require **both** mates passing, proper-pair flags and
same reference sequence. Deduplicate by chromosome and outer reference-coordinate
span. Single-end Ion counts use chromosome/start/end, ignoring strand. This is
coordinate deduplication, not proof of independent original molecules. See
[method fidelity](docs/methods_consistency.md) for exact CIGAR/NM semantics.

Burden is strict units per million normalized **input pairs** or **input reads**,
according to layout. Chromosome representation is not base-by-base genome coverage.
Kraken human calls are corroboration; they do not define the strict count.
Kraken2 is a required validated pipeline component, not an optional development
dependency: the runner requires its executable/database, and strict-positive
runs invoke it. The real database wrapper is tested separately with artificial
read pairs; this does not establish performance for every positive dataset.
Exact assembly challenge requires >=80% query coverage and >=95% identity in at
least one mate/read. Same-species substitutions are not accepted. `PASS`,
`BACTERIAL_MATCH_PRESENT`, `NO_EXACT_ASSEMBLY` and `FAILED` are distinct statuses.
The aggregate exporter does not convert failed/missing challenges into confirmed
human positives or measured negatives. Historical exploratory adjudication used
a separate rule, documented in the method fidelity table.

## Resources and privacy

The human FASTA is about 3.16 GB uncompressed, with a multi-GB minimap2 index.
Kraken memory depends strongly on database choice; reads/SRA and temporary BAMs
can require tens to hundreds of GB per active library. Start with one global
analysis slot and monitor disk/RAM. The example has four AWS fetch workers, two
AWS analysis workers constrained by one global science slot. This is configurable,
not a performance guarantee or a universal resource estimate.

Generated run-level outcomes, metadata, strict reads and genetic evidence can be
sensitive even when inputs are public. Outputs default to user-only permissions.
Do not commit the workspace. Aggregate output still needs small-cell/linkage
review. Pseudonyms are not anonymity guarantees. See [privacy](docs/privacy.md).

## Limitations

Validation includes synthetic tests and one bounded paired-end public-archive
integration. Full cross-platform integration and large positive/negative cohort
regressions were not performed. See [validation](docs/validation.md) for coverage.
Database, tool and archive revisions can change outputs. Exact historical selections
and some historical analysis producers remain restricted or unrecovered. Reproducing
the algorithms on independently selected data does not reproduce the exact historical
manuscript populations or results.

## License and citation

The original pipeline code is licensed under the [MIT License](LICENSE),
Copyright (c) 2026 Mahdi Moosa. Separately installed executables and reference data
retain their own licenses and terms; see [dependencies](docs/dependencies.md).

Please cite the software title, author and version using [CITATION.cff](CITATION.cff).
