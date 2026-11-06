# Final-methods consistency and evidence strength

| Component | Implemented definition | Evidence/status |
|---|---|---|
| Human reference | CHM13v2 masked Y PAR plus rCRS, pinned FASTA | frozen configuration + retained preparation log; no new reference download |
| Mapper | minimap2 -ax sr; samtools sort | original production functions unchanged |
| Primary records | samtools -F 0x904 | original functions; synthetic actual-samtools smoke |
| Strict QC | MAPQ>=30, coverage>=80%, identity>=98% | unchanged strict QC predicates and boundary tests |
| Paired fragment | both mates pass, proper pair, same chromosome, outer-span dedup | original predicate; synthetic controls |
| Single Ion | one qualifying read, chromosome/start/end dedup ignoring strand | original predicate; separate synthetic tests |
| Depth normalization | normalized input pairs/read count after source filtering | original acquisition/statistics functions |
| Source read filter | biological PASS/CRITERIA retained, REJECT excluded; technical-index semantics preserved | retained acquisition module, packaging comment-only edit |
| Kraken | human taxon 9606 corroboration, not strict positivity | original scientific calls |
| Bacterial challenge | same BioSample/isolate; >=80% query coverage, >=95% identity in either mate/read | original assembly resolution + challenge; fake-process boundary test |
| Resume | result +25 chromosome rows and matching sum; compatible method | preserved master, synthetic commit/interruption tests; added conservative release-identity guard |
| Metadata eligibility | direct isolate evidence, explicit contradictions and identity resolution | later reviewed metadata module, not claimed original E producer |
| Cohort selection | metadata-only Hamilton/HMAC algorithm, user-supplied seed | later deterministic selector and synthetic order-invariance test; country cross-check external |
| Population SNPs | separate analysis package required; final five-population ascertainment and scoring | final migration provenance; no older four-population implementation bundled |

## Exact implementation semantics, not idealized formulas

For a CIGAR string, q consumes MI=X, alignment denominator a consumes MID=X,
reference span consumes MDN=X. Coverage is 100*q/len(SEQ); identity is
100*(a-NM)/a. Soft clips are excluded from q. The historical parser uses NM=0
when no NM tag is found; minimap2 normally supplies NM. These semantics, including
that fallback, have not been silently changed. Require the documented mapper
rather than importing arbitrary BAMs with missing NM or unsupported hard clips.
Primary filtering happens in samtools, not a second Python flag parser. Proper
pair flags are required without a new insert-size cutoff. Single-end strict
counts do not establish physical pair support. Coordinate deduplication is not
UMI/molecule deduplication and does not establish distinct contributors.

## Adjudication is not one universal Boolean

Historical E adjudicated positivity was strict-positive and Assembly_check !=
BACTERIAL_MATCH_PRESENT, including historical failed assembly checks. This rule
must not be mislabeled exact-assembly validation. Other manuscript populations
have endpoint-specific definitions. The public core retains the component fields
and aggregate statuses; it does not declare failed checks confirmed human origin.
Choosing a downstream cohort endpoint requires its own stated policy and denominator.

## Evidence and limits

Current paired analysis.py SHA-256:
`32163a3f4b986ab02b182a01d7106d5b2740225369c026b847b980f0b666f74f`. Historical source SHA-256 was
`eb07123bf5412346214797a30d4c6efc6c4c3bb8b9b24ae491343d3b442db122`.
Ion ion_analysis.py SHA-256:
`be356f0803d8a8de5753da2dff79aaddb5c80b19392e90c20cf84e5c62bd6a1b`.
Ion remains byte-identical to the validated original. Paired scientific predicates
are unchanged; only the required bacterial samtools-view error-status guard was added. Worker
science function ASTs are unchanged. Operational portability/security changes are
recorded in source_lineage.json and the separate review report.

The old standalone paired-only pipeline is superseded for multi-platform
orchestration. The later end-to-end package ports these scientific functions but
includes strict-core, native-import and reference-competition proposals with
separate scientific scope; they are not silently merged. The final manuscript
analyses consume frozen classifications and separately produced marker outputs.

No claim is made that passing synthetic tests reproduces actual historical run
counts. One bounded paired-end E. coli ENA integration and a real-database artificial
Kraken wrapper test passed (docs/validation.md). Full historical positive/negative
paired and Ion, SRA/AWS/Aspera and real-positive exact-assembly golden integration
remain untested; no universal validation is claimed.
