# Release-candidate validation

Validated locally on 3 November 2026; no publication, GitHub repository, commit,
remote or release tag was created. Mahdi Moosa confirmed the MIT copyright and
software authorship. A manuscript DOI is not required for software use/release.

## Fresh installation

A new Conda prefix was created with `conda env create --prefix NEW_PREFIX --file
environment.yml --solver libmamba`. It was not cloned from or installed into the
research environment. Conda may reuse downloaded package-cache objects; that is
not environment reuse. Python 3.14.10, minimap2 2.31-r1302, samtools/htslib1.24,
SeqKit2.13.0, Kraken2 2.17.1 and SRA Toolkit3.4.1 were present. Exact Conda package
builds/URLs/checksums and commands are retained in the separate review records.
The wrapper installed from its own source with no scientific algorithm changes.
The build requires setuptools>=77 for SPDX license metadata.

From the repository root, use the documented Conda procedure, then:

```bash
python -m pip install .
python -B -m unittest discover -s tests -v
python examples/synthetic_smoke.py
```

All 20 packaged tests passed, including actual samtools primary filtering and an
actual minimap2/samtools challenge of artificial read pairs against an artificial
bacterial assembly. Test records contain no reference-derived human sequence.
Configuration parsing, missing-reference handling, layout/QC boundaries, output
accounting, completion/resume and deterministic selection tests passed.

## One bounded public-archive integration

A single E. coli isolate paired-end Illumina WGS library was selected using
ordinary public metadata and smallest advertised compressed size, after the
historical sampling-frame cutoff. No private cohort list or contamination outcome
was consulted. The download was about 54 MB compressed, under a declared 250 MB
budget. This budget is a validation resource bound, not scientific eligibility.
It is not an attempt to find a positive sample or reproduce a historical cohort.

The installed public package completed metadata-only eligibility classification,
ENA metadata retrieval and FASTQ acquisition, MD5 checking, paired record
validation, mapping to the verified CHM13 masked-Y/rCRS reference, strict counting,
25-row chromosome accounting, result-schema generation, checkpoint completion and
aggregate export. Repeating the same run left the scientific result/chromosome
files byte-identical. An exact-BioSample assembly was retrieved by the existing
resolver, but a completed candidate-read assembly challenge was not exercised by
this library. Its algorithmic alignment path passed the real-tool artificial
assembly test instead. No positive sample or successful clinical/genetic recovery
was required or sought. Individual human evidence and run-linked outcomes are not
public validation fixtures.

The externally provisioned human FASTA/index and Kraken database were checksum
verified. No genomes, Kraken database or public read dataset is redistributed.
No full cohort, model, Hostile workflow, variant extraction or figure was run.

## Kraken scope

Kraken2 is a required validated core component. It corroborates strict candidate
reads and is invoked by the positive branch; it does not define strict positivity.
A standalone call through the packaged paired-read wrapper processed one wholly
artificial pair against the real PlusPFP-16 database successfully. This establishes
wrapper/database execution, not sensitivity/specificity or a golden regression
for every real positive library. The no-completed-challenge integration case does
not replace that separate check.

## Limits

Only one paired-end ENA library was integrated. Full real Ion Torrent, SRA/AWS
acquisition, Aspera and positive Kraken/exact-assembly golden suites were not run.
Their original implementations are preserved and the supported scientific
predicates have synthetic tests; no broader live-integration coverage is claimed.
Original historical cohort rosters/seeds/private genetic observations are not
included. Country agreement for stratified selection and external reference assets
remain explicit prerequisites. Companion manuscript analyses remain separate.

Full checksums, artifact inventories, isolated installation and disclosure audit
results are supplied in the local release-review directory, outside this repo.
Passing scans and synthetic tests does not guarantee anonymity or universal
scientific equivalence across changing tools, databases and archive records.

## 5 November technical-failure correction

The original 20 tests and nine new failure-propagation controls pass (29 total). Empty successful samtools streams remain valid negative challenges. Nonzero/partial failed streams raise technical errors and are reported FAILED/NA. Worker/master/aggregate controls preserve explicit failed status; resume completion is strict-processing completion, not certification of a failed assembly challenge. No scientific QC predicate changed. Both rebuilt distributions are installed and tested in a fresh isolated Python environment, jointly with SRAFEGUARD. Original dated integration evidence remains historical; no archive download or cohort reprocessing was repeated.
