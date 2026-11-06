# External references

No genomes, VCFs, human sequences, databases or clinical resources are bundled.

## Human mapping analysis set

Source recorded in the original preparation log:
`https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz`

This is CHM13v2.0 with the completed Y chromosome, masked Y pseudoautosomal
regions and rCRS mitochondrion. It is not a reference with the entire Y masked.
The retained upstream download records give compressed MD5
`b17cd2bdd8eb5544e195d6a8882e417b` and compressed size 980,387,160 bytes.
The retained frozen reference configuration identifies decompressed FASTA SHA-256:
`a5128da5a5729d1ef92bf2040494092d784a5ff07e4d08f7b9f385b5498d84d3`.
These are preserved provenance values; remote availability was not re-probed in
this packaging audit.

Example external preparation (not automatically executed by installation):

```bash
mkdir -p references
curl --fail --location --continue-at - \
  --output references/chm13v2.0_maskedY_rCRS.fa.gz \
  https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_maskedY_rCRS.fa.gz
printf '%s  %s\n' b17cd2bdd8eb5544e195d6a8882e417b references/chm13v2.0_maskedY_rCRS.fa.gz | md5sum --check
gzip -dc references/chm13v2.0_maskedY_rCRS.fa.gz > references/chm13v2.0_maskedY_rCRS.fa
sha256sum references/chm13v2.0_maskedY_rCRS.fa
minimap2 -x sr -d references/chm13v2.0_maskedY_rCRS.sr.mmi references/chm13v2.0_maskedY_rCRS.fa
samtools faidx references/chm13v2.0_maskedY_rCRS.fa
sha256sum references/chm13v2.0_maskedY_rCRS.sr.mmi
```

Compare the FASTA hash with the pinned value above. Configure your verified local
index hash. The retained historical short-read index hash is
`2158acb5aef0f7b9e8c09fcab1ac7e3c9a7281817f0a9a081380074580a72555`;
an index rebuilt with a different executable may differ. Do not substitute a
map-ont/asm index. The provenance of the index-to-FASTA construction must be
retained externally; checking a supplied index hash alone does not prove its source.

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
