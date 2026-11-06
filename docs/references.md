# External references

Reference genomes, minimap2 indexes, Kraken2 databases, population VCFs and
clinical resources must be provisioned separately; they are not bundled with
software version 0.1.0rc1.

## Human mapping analysis set

Use the **T2T-CHM13v2.0 masked-Y + rCRS** analysis set from the
[official T2T resource page](https://genomeinformatics.github.io/CHM13v2/):
[`chm13v2.0_maskedY_rCRS.fa.gz`](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz).

This analysis set contains the completed Y chromosome with its pseudoautosomal
regions hard masked with Ns, and uses rCRS for the mitochondrial sequence. It
neither excludes Y nor masks the entire chromosome. Do not substitute the unmasked-Y,
no-Y or other mitochondrial-reference analysis sets.

The retained reference preparation and validation records identify:

| Resource | Recorded identity |
|---|---|
| Compressed FASTA | MD5 `b17cd2bdd8eb5544e195d6a8882e417b`; 980,387,160 bytes |
| Decompressed FASTA | SHA-256 `a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3` |
| Historical short-read minimap2 index | SHA-256 `2158acb5aef0f7b9e8c09fcab1ac7e3c9a7281817f0a9a081380074580a72555` |

These are retained validation values, not newly calculated hashes. The official
resource page and a header-only check of the FASTA download URL were verified on
6 November 2026; the genome was not downloaded for that check. The S3 ETag is not
the documented MD5 checksum.

### Download and index

Run from the repository root. Keep references outside all workspace and cleanup
roots. These commands download the external reference; installation does not do so.
Stop if a checksum check fails.

```bash
mkdir -p references
curl --fail --location --continue-at - \
  --output references/chm13v2.0_maskedY_rCRS.fa.gz \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
printf '%s  %s\n' b17cd2bdd8eb5544e195d6a8882e417b references/chm13v2.0_maskedY_rCRS.fa.gz | md5sum --check
gzip -dc references/chm13v2.0_maskedY_rCRS.fa.gz > references/chm13v2.0_maskedY_rCRS.fa
printf '%s  %s\n' a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3 references/chm13v2.0_maskedY_rCRS.fa | sha256sum --check
minimap2 -x sr -d references/chm13v2.0_maskedY_rCRS.sr.mmi references/chm13v2.0_maskedY_rCRS.fa
samtools faidx references/chm13v2.0_maskedY_rCRS.fa
sha256sum references/chm13v2.0_maskedY_rCRS.sr.mmi
```

Use the documented minimap2 version in [environment.yml](../environment.yml).
Record your verified local index hash: an index rebuilt with a different executable
may differ from the historical hash. Do not use a `map-ont` or `asm` index. Retain
the index-building command, tool version and source FASTA identity externally;
checking an index hash alone does not prove which FASTA or preset produced it.

### Configure and validate

Copy [configs/pipeline.example.ini](../configs/pipeline.example.ini) to
`configs/local.ini`. Set these entries, retaining the other example settings:

```ini
[references]
human_index = ../references/chm13v2.0_maskedY_rCRS.sr.mmi
kraken_db = ../references/kraken2

[reference_identity]
human_fasta = ../references/chm13v2.0_maskedY_rCRS.fa
human_index_sha256 = SET_VERIFIED_64_CHARACTER_SHA256
kraken_identity_file = ../references/kraken2_identity.json
```

Replace `SET_VERIFIED_64_CHARACTER_SHA256` with the verified SHA-256 of your local
short-read index. Paths are relative to the INI file. The runtime validator checks
the decompressed FASTA against its fixed validated hash above and the index against
`human_index_sha256`; setting an arbitrary FASTA hash in the configuration does not
permit another reference build. Provide the separately verified Kraken2 database
and identity manifest described below before full validation.

```bash
wgs-conta validate --config configs/local.ini --paths-only
wgs-conta validate --config configs/local.ini
```

## Kraken2 and exact bacterial assemblies

The manuscript used frozen PlusPFP-16, release 9 December 2022, with Kraken2
2.17.1. This package does not retrieve or redistribute that database. Obtain it
from the database provider under its terms, retain an exact version and SHA-256
manifest for hash/opts/taxo files, and set `kraken_identity_file` to that JSON or
text manifest. Its fingerprint is used to reject silent resume incompatibility;
verify it against the actual files before running. A different database is a
new analysis and is not exact manuscript reproduction.

Exact bacterial assemblies are resolved against the same BioSample through
NCBI, with ENA fallback. Optional user-curated tables must establish exact
isolate/BioSample matching. Never substitute a species representative. API and
assembly release drift can change availability. Save downloaded assembly identity
and integrity records with each independent analysis.

## Population-analysis dependency (not implemented here)

Use complete five-superpopulation ascertainment and scoring: AFR, AMR, EAS, EUR,
SAS. The final complete autosomal resource is named
`chm13_1kgp_five_population_delta010.tsv.gz`, SHA-256
`6e5655a095f1b071915e5c5f527a4a698e96dcf08a6cd58cba141f7a7ecdc112`.
Its Delta-AF thresholds .10/.20/.30 retain 6,947,747 / 4,013,865 / 2,156,198 SNPs.
The upstream release is the T2T CHM13v2.0 1KGP unrelated_samples_2504 allele_freq
sites resource, autosomes 1--22, biallelic A/C/G/T SNVs and preserved ambiguity
filters. The upstream base is:
`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/unrelated_samples_2504/allele_freq/`
Chromosome filename: `1KGP.CHM13v2.0.chrN.recalibrated.snp_indel.pass.withafinfo.vcf.gz`.

Construction, full upstream object checksums and scoring need the separately
reviewed manuscript-analysis package. Merely recalculating Delta AF in an old
four-population subset is incomplete. True AF=0 can be recovered from explicit
AC=0 with AN>0; AN=0 remains missing. Do not silently replace this dependency with
a differently defined reference or bundle allele observations.

## Metadata resources

ISO-3166 and UN M49 tables are optional external inputs to geography-stratified
selection, not patient/source-country mappings. Record their version and checksum
and verify redistribution terms. Do not bundle institution-specific mappings.
