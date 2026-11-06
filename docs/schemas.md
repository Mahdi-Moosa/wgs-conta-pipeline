# Input/output contracts

## Archive capture

`configs/archive.example.json` provides the effective configuration keys; the
capture CLI rejects missing keys, reversed dates and unsupported full platform
lists. Cutoff is an inclusive calendar date implemented as an exclusive following
midnight timestamp. `FIELDS`, `PLAN_COLUMNS` and `MANIFEST_COLUMNS` in `archive.py`
are executable schemas. Raw TSV parsing preserves ENA's literal tabular semantics.
Count-checks and drift failures block snapshot certification.

## Metadata classification

Input: raw ENA columns run_accession, study_accession, secondary_study_accession,
sample_accession, secondary_sample_accession, scientific_name, tax_id, country,
strain, isolate, sample_title, study_title, instrument_platform, instrument_model,
library_strategy, library_source, library_layout, first_public. Identifier aliases
must resolve atomically; conflicting/multi-sample records are excluded.
Output: the original fields plus `metadata_rules.DERIVED_FIELDS`, including
primary_eligible, exclusion_codes, direct_isolate_reason and canonical identities.
This is private local metadata, not an approved public roster.

The optional taxon deny-list JSON has `scientific_names` and/or `tax_ids` arrays.
No actual private taxon/cohort list is bundled. Define exclusions before outcomes.
Metadata rules cannot detect undocumented host-cell workflows; retain expert review.

## Reviewed selection frame

Required TSV columns: run_accession, canonical_bioproject, canonical_biosample,
scientific_name, instrument_model, instrument_platform, sample_country_iso3,
primary_eligible, country_compare_status. Optional first_public.
`country_compare_status=MATCH` must come from an actual independent country
cross-check; this package does not generate that evidence. Country ISO3 maps to
one M49 subregion per eligible BioSample and project. Multi-subregion/unmapped
projects are excluded, cell = project x scientific name x instrument x subregion.
At least the requested number of distinct BioSamples is required. Hamilton
allocation and HMAC ranking use only metadata plus the user's 32-byte seed.
A project selected with insufficient runs per BioSample causes a failure, not an
outcome-dependent replacement. Output is a private registry, not a published
accession/pseudonym crosswalk. No reserve/substitution policy is fabricated.

M49 TSV schema for selection: iso3, region_un, subregion_un. Classification also
needs country_un and a compatible ISO JSON object with key `3166-1` containing
alpha_3/name/alpha_2 records. These external reference tables are not bundled.

## Pipeline registry

Required TSV header `Accession` with atomic SRR/ERR/DRR identifier; optional Label,
Cohort, Active, Priority (integer), Preferred_source (auto/aws/ncbi/ena), Notes.
No real registry example is bundled. Unknown columns remain local, not exported.
Duplicates and invalid active/source values fail before execution.

## Configuration

`config.py::Settings` is the executable INI schema. Sections: paths, references,
acquisition, workers, assembly, cleanup; optional cache; release reference_identity.
All path fields resolve relative to the INI location. reference_identity requires
human_fasta, human_index_sha256 and kraken_identity_file. Worker and thread counts
must be positive. No original sampling seed is a default configuration field.

## Outputs

The master directory contains master_results.tsv, chromosome_counts.tsv (25 rows
per completed run/version), cohorts.tsv and commit_log.tsv. `master.py::CHROMS`
defines the expected reference set. worker.py::analyze_paired/analyze_ion define
all result fields. Primary fields:

- Analysis_branch and Quantification_unit: PAIRED_SHORTREAD/PAIRS or ION_SINGLE/READS.
- Total_pairs / Total_reads: normalized input depth with nonapplicable fields NA.
- Unique_strict_units, Unique_strict_fragments/reads: coordinate-deduplicated counts.
- Strict_per_million_units: count/depth x 1e6, in the appropriate unit.
- Human_positive: strict count >0, never a clinically validated human genotype.
- Assembly_check: PASS / BACTERIAL_MATCH_PRESENT / NO_EXACT_ASSEMBLY / FAILED /
  NOT_NEEDED_NEGATIVE. Missing/failed status is not a measured negative challenge.
- Bact_ge1_ID95: reads/pairs with at least one >=80% coverage, >=95% identity hit.
- Pipeline_version and Scientific_method_version: separate package and science IDs.

Detailed evidence can contain private strict-read identifiers and sequences;
never automatically publish it. `summarize` emits counts and statuses only,
requires exactly one result version/run, validates positivity against counts,
and keeps PAIRS/READS and failed/missing challenge states separate. Small-cell
review is still required. There is no automatic per-run publication export.

## Technical assembly failures

A nonzero samtools view exit, even after partial output, raises a technical error. The worker records Assembly_check=FAILED with bacterial-hit counts NA; it never treats missing/failed output as a completed negative challenge. Master/resume completion refers to completed strict processing: FAILED assembly rows remain explicitly FAILED when resumed, and summarize excludes them from completed/no-bacterial-match challenges. To retry adjudication, use a new isolated analysis workspace. Successful empty streams retain the original PASS criteria.
