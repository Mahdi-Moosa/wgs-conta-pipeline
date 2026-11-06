# Privacy and release boundary

Public source code does not make linked human evidence safe to publish. A
public bacterial accession can become sensitive when linked to human alleles,
splice evidence, inferred affinity or clinical annotations. No such observations
or historical selected lists are included here.

The included code takes user-selected records and intentionally writes local
run-level/evidence outputs for audit. They are private research artifacts, not
automatically releasable tables. Entry-point outputs use umask 077; keep a
user-only workspace on encrypted/controlled storage if appropriate. Do not commit
outputs, seeds, mapping keys, evidence BAM/FASTQ, logs or configuration secrets.
A .gitignore does not protect files manually added or uploaded elsewhere.

Release algorithms, explicit configuration schemas and independently configurable
selection procedures. Withhold original sampling seeds and exact historical
privacy-sensitive rosters. Public methods can be replicated on newly selected,
outcome-blind samples but need not reproduce the historical sample or estimates.
No seed value from the original project is included. The artificial unit-test key
is clearly synthetic and never used by a production default.

The source allowlist excludes original Git metadata, data freezes, notebooks,
credentials, machine logs, private paths, pseudonym crosswalks, country mappings,
real genetic profiles, raw reference genomes and read data. Dedicated checks cover
source, configs, fixtures, wheel contents, hidden files, symlinks and unexpected
large files. Pattern scanning alone does not prove safety. Numerical fingerprints,
small cells and repeated aliases need human cross-table review before any data
release. Do not infer individual donor identity, ancestry, genotype or diagnosis.

Before a first public commit, review the author-confirmed MIT/citation metadata,
review the dependency boundary and final file manifest, inspect Git's staged
files (if initialized), and re-run disclosure checks on the exact release archive.
This candidate intentionally has no Git repository or remote.
