# Tools, attribution and licensing

The author selected MIT and confirmed Mahdi Moosa as copyright holder and
software author. LICENSE contains the standard MIT text. No institutional
ownership or approval is asserted. Incorporated local source modules have no
identified third-party source headers; external tools/data are not vendored.

External executables are called as subprocesses, not bundled:

| Tool | Recorded workflow version | Purpose |
|---|---|---|
| minimap2 | 2.31 | short-read human and bacterial mapping |
| samtools | 1.24 | BAM sorting and primary-alignment stream |
| seqkit | 2.13.0 | FASTQ statistics and exact QNAME extraction |
| Kraken2 | 2.17.1 | candidate-read taxonomic corroboration |
| SRA Toolkit | 3.4.1 | SRA validation, biological read extraction/filter accounting |
| pigz | 2.8 | optional compression/reference preparation |
| EDirect | version not pinned by retained environment | exact-BioSample assembly/metadata queries |
| NCBI Datasets CLI | version not pinned by retained environment | assembly retrieval |
| curl | version not pinned by retained environment | API/file transport |
| AWS CLI | optional; version must be recorded by user | anonymous SRA mirror acquisition |

The package also supports optional user-installed Aspera; no binary, bypass key
or license is supplied. No credential or private API key is needed by defaults.
ENA, NCBI and AWS archive terms/rate limits apply independently. Do not download
protected data without authorization.

Tool software licenses must be respected separately: inspect the exact upstream
releases and their notices before redistribution. These descriptions are not a
license clearance for a bundled environment. Review transitive dependencies and
reference-provider terms; this candidate bundles none of the binaries/databases.

Citations should include the minimap2, SAMtools, Kraken2, SeqKit and SRA/ENA
resource publications relevant to actual use. Software authorship is author-confirmed in CITATION.cff. No DOI or repository
URL is invented; the manuscript DOI may be added when available.

## Licenses reported by the installed Conda distributions

The fresh environment package metadata reports the following identities. These
are separate from the wrapper MIT license; exact upstream and transitive notices
remain authoritative. No executable or reference data is bundled.

| Dependency | Installed version | Package-reported license |
|---|---|---|
| minimap2 | 2.31 | MIT |
| samtools | 1.24 | MIT |
| seqkit | 2.13.0 | MIT |
| kraken2 | 2.17.1 | GPL-3.0-or-later |
| sra-tools | 3.4.1 | Public Domain |
| pigz | 2.8 | Zlib |
| entrez-direct | 26.0 | Public Domain |
| ncbi-datasets-cli | 18.40.0 | Public Domain |
| curl | 8.23.0 | curl |
| awscli | 2.38.3 | Apache-2.0 |
| python | 3.14.10 | Python-2.0 |
| setuptools | 84.2.1 | MIT |
