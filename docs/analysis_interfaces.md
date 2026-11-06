# Boundary of the later manuscript-analysis repository

The core emits normalized strict counts, chromosome tables, source/read-filter
metadata, assembly statuses and private compact strict evidence. It does not
recompute manuscript analyses or publish accession-linked genetic observations.

| Analysis | Core support | Separate producer/interface needed |
|---|---|---|
| Project recurrence/ICC/MOR | component classifications + cohorts | cohort-specific statistical models, posterior conventions and full configs |
| Temporal and technical analyses | metadata + strict outputs | frozen project-balanced/technical selection and contrasts |
| Five-population affinity | retained strict candidate reads | complete five-population resource, observation extraction, likelihood and branch-specific bootstrap |
| SNP configurations | strict evidence interfaces | deterministic per-locus caller and original selection rules |
| Within-project state concordance | no pair enumeration or allele profiles | eligible clean-state intersections, duplicate exclusions, private comparisons |
| ClinVar/P/LP | no clinical module | frozen ClinVar source, candidate screen and read-pair placement challenge |
| SRAFEGUARD mitigation | original strict outputs can be independent inputs | PRE/POST audit, endpoint applicability and frozen affinity resource |
| Hostile/assembly/BUSCO | not implemented | separate filtering, assembly, annotation and completeness workflows |
| RNA splice evidence | not implemented | HISAT2/annotation-free validation, QNAME/junction units |
| Figures/Source Data | no manuscript artifacts | privacy-approved aggregates and exact final renderers |

Final affinity branches use AFR/AMR/EAS/EUR/SAS for both Delta AF ascertainment
and scoring. Figure 2a primary/sensitivity thresholds are .10/.30, reporting
eligibility >=20 clean loci, 1,000 within-run chromosome-qualified observed-1Mb
block bootstraps (scorer computation gate 10 loci/two blocks). ED5d uses .10/.20/.30
and 5,000 run-by-chromosome-qualified-block resamples, with no per-run 20-locus
reporting gate. SRAFEGUARD uses .30 and its own 1,000 within-run implementation.
These are different interfaces; no universal bootstrap implementation is implied.

The final SNP-recovery cohort uses five-population .10 markers and original
>=20 clean-locus/project-capping selection. Configuration thresholds and allele
concordance are distinct from population affinity and from diploid genotypes.
The original historical per-locus configuration caller and historical ED5
bootstrap producer remain unrecovered; maintained/replacement implementations
have their own validation provenance. Do not describe reconstructed producers
as the original programs.

Do not publish selected rosters, original seed values, accession/pseudonym
crosswalks, run-level allele vectors, clinical profiles or private mappings in
that later repository. Release executable algorithms/configuration schemas and
safe synthetic fixtures, plus approved aggregate data only. Reference resources
remain external with checksums, attribution and retrieval instructions.
