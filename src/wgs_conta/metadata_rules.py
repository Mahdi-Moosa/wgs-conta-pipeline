"""Outcome-blind, auditable ENA run metadata normalization rules."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import unicodedata
from pathlib import Path

PROJECT_RE = re.compile(r'PRJ(?:NA|EB|DB)\d+', re.I)
SAMPLE_RE = re.compile(r'SAM(?:N|EA|D)\d+', re.I)
PROJECT_ALIAS_RE = re.compile(r'(?:SRP|ERP|DRP)\d+', re.I)
SAMPLE_ALIAS_RE = re.compile(r'(?:SRS|ERS|DRS)\d+', re.I)
EMPTY = {'', 'na', 'n a', 'n/a', 'none', 'null', 'unknown', 'not known', 'not provided',
         'not available', 'not applicable', 'not specified', 'not reported', 'not recorded',
         'not determined', 'not collected', 'unidentified', 'none given', 'no data',
         'missing', '-', 'unspecified', 'unavailable', 'other', 'strain', 'isolate'}
SPECIMEN_ONLY = re.compile(r'^(?:urine|blood|stool|faeces|feces|sputum|saliva|serum|plasma|skin|tissue|wound|swab|nasal swab|rectal swab|throat swab|vaginal swab|environment|environmental sample|clinical sample|clinical specimen)$', re.I)
SUPPORTED = {'ILLUMINA': ('PAIRED', 'ILLUMINA'), 'BGISEQ': ('PAIRED', 'BGI_MGI'),
             'DNBSEQ': ('PAIRED', 'BGI_MGI'), 'ION_TORRENT': ('SINGLE', 'ION_TORRENT')}
LABEL_PATTERNS = [
    ('LABEL_METAGENOME', re.compile(r'\bmetagenom(?:e|ic|ics)\b', re.I)),
    ('LABEL_MICROBIOME', re.compile(r'\bmicrobiom(?:e|ic)\b', re.I)),
    ('LABEL_UNCULTURED', re.compile(r'\buncultured\b', re.I)),
    ('LABEL_SYNTHETIC', re.compile(r'\b(?:synthetic\s+(?:construct|genome|DNA)|artificial\s+construct)\b', re.I)),
    ('LABEL_MIXED_AGGREGATE', re.compile(r'\b(?:mixed\s+(?:culture|sample|material|community)|pooled\s+(?:samples?|isolates?)|aggregate\s+(?:sample|material))\b', re.I)),
]
CONTRADICTION_PATTERNS = [
    ('DIRECT_SPECIMEN_SEQUENCING', re.compile(r'\b(?:direct(?:ly)?\s+(?:sequenc\w*|DNA\s+extraction)\s+(?:from|of)\s+(?:clinical\s+)?(?:blood|stool|feces|faeces|sputum|urine|tissue|swab|specimen|sample)|(?:blood|stool|feces|faeces|sputum|urine|tissue|swab)\s+direct\s+sequenc\w*)\b', re.I)),
    ('CULTURE_FREE_ENRICHMENT', re.compile(r'\b(?:culture[- ]free|without\s+cultur(?:e|ing)|no\s+cultur(?:e|ing)|target(?:ed)?\s+enrichment\s+(?:directly\s+)?from\s+(?:clinical\s+)?(?:specimen|sample|blood|stool))\b', re.I)),
    ('HOST_CELL_DEPENDENT_WORKFLOW', re.compile(r'\b(?:host[- ]cell\s+(?:dependent|culture|propagation)|propagat\w*\s+in\s+(?:host|mammalian|human)\s+cells?)\b', re.I)),
]


def norm_key(value: str) -> str:
    value = unicodedata.normalize('NFKD', value or '')
    value = ''.join(c for c in value if not unicodedata.combining(c))
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9]+', ' ', value.casefold())).strip()


def meaningful(value: str) -> bool:
    return norm_key(value) not in EMPTY


def resolve_identity(primary: str, secondary: str, kind: str) -> tuple[str, str]:
    canonical = PROJECT_RE if kind == 'project' else SAMPLE_RE
    alias = PROJECT_ALIAS_RE if kind == 'project' else SAMPLE_ALIAS_RE
    found = set()
    for raw in (primary, secondary):
        value = (raw or '').strip().upper()
        if not meaningful(value):
            continue
        if re.search(r'[,;|]', value) or len(value.split()) != 1:
            return '', 'NON_ATOMIC'
        if canonical.fullmatch(value):
            found.add(value)
        elif not alias.fullmatch(value):
            return '', 'CONFLICT'
    if len(found) > 1:
        return '', 'CONFLICT'
    if not found:
        return '', 'MISSING'
    return next(iter(found)), 'RESOLVED'


class CountryReference:
    def __init__(self, m49_path: Path, iso_path: Path):
        self.by_iso = {}
        self.aliases = {}
        with open(m49_path, newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f, delimiter='\t'):
                iso = row['iso3']
                self.by_iso[iso] = row
                self.add_alias(row['country_un'], iso)
                self.add_alias(iso, iso)
        with open(iso_path, encoding='utf-8') as f:
            for row in json.load(f)['3166-1']:
                iso = row['alpha_3']
                if iso not in self.by_iso:
                    continue
                for field in ('name', 'official_name', 'common_name', 'alpha_2'):
                    if row.get(field):
                        self.add_alias(row[field], iso)
        manual = {'usa': 'USA', 'u s a': 'USA', 'united states': 'USA',
                  'uk': 'GBR', 'great britain': 'GBR', 'england': 'GBR',
                  'scotland': 'GBR', 'wales': 'GBR', 'northern ireland': 'GBR',
                  'russia': 'RUS', 'south korea': 'KOR', 'north korea': 'PRK',
                  'iran': 'IRN', 'vietnam': 'VNM', 'czech republic': 'CZE',
                  'ivory coast': 'CIV', 'dr congo': 'COD',
                  'democratic republic of congo': 'COD', 'republic of congo': 'COG',
                  'tanzania': 'TZA', 'bolivia': 'BOL', 'venezuela': 'VEN'}
        for name, iso in manual.items():
            if iso in self.by_iso:
                self.add_alias(name, iso)
        self.ambiguous = {'congo', 'korea', 'multiple countries', 'various countries',
                          'worldwide', 'global', 'international', 'asia', 'africa',
                          'europe', 'americas', 'oceania'}

    def add_alias(self, name: str, iso: str):
        self.aliases.setdefault(norm_key(name), set()).add(iso)

    def normalize(self, raw: str) -> dict:
        text = (raw or '').strip()
        result = dict(country_raw=raw or '', sample_country_normalized='',
                      sample_country_iso3='', sample_region_candidate='',
                      sample_geo_status='UNRESOLVED', sample_geo_source='ENA_READ_RUN_COUNTRY')
        if not meaningful(text):
            return result
        if any(x in text for x in (';', '|')):
            result['sample_geo_status'] = 'AMBIGUOUS'
            return result
        # ENA country often carries a locality after the first colon.
        candidate = text.split(':', 1)[0].strip()
        key = norm_key(candidate)
        if key in self.ambiguous:
            result['sample_geo_status'] = 'AMBIGUOUS'
            return result
        matches = self.aliases.get(key, set())
        if len(matches) != 1:
            result['sample_geo_status'] = 'AMBIGUOUS' if matches else 'UNRESOLVED'
            return result
        iso = next(iter(matches))
        ref = self.by_iso[iso]
        display = {'USA': 'United States', 'GBR': 'United Kingdom'}.get(iso, ref['country_un'])
        result.update(sample_country_normalized=display, sample_country_iso3=iso,
                      sample_region_candidate=ref['region_un'], sample_geo_status='RESOLVED')
        return result


class TaxonList:
    """Frozen explicit membership: exact stripped names OR exact tax IDs."""

    def __init__(self, path: Path):
        raw = path.read_bytes()
        data = json.loads(raw.decode('utf-8-sig'))
        if not isinstance(data, dict) or set(data) - {'scientific_names', 'tax_ids'}:
            raise ValueError('taxon list must be an object with scientific_names and/or tax_ids')
        for field, values in data.items():
            if not isinstance(values, list) or any(
                    not isinstance(v, str) or not v.strip() for v in values):
                raise ValueError(f'{field} must be a list of nonempty strings')
        self.names = frozenset(v.strip() for v in data.get('scientific_names', []))
        self.tax_ids = frozenset(v.strip() for v in data.get('tax_ids', []))
        if any(not re.fullmatch(r'[1-9][0-9]*', v) for v in self.tax_ids):
            raise ValueError('tax_ids must contain positive decimal ID strings without leading zeros')
        if not self.names and not self.tax_ids:
            raise ValueError('explicit taxon list must not be empty')
        self.provenance = dict(path=str(path.resolve()), sha256=hashlib.sha256(raw).hexdigest(),
                               scientific_names=sorted(self.names), tax_ids=sorted(self.tax_ids))

    def matches(self, row: dict) -> bool:
        return (row.get('scientific_name', '').strip() in self.names or
                row.get('tax_id', '').strip() in self.tax_ids)


def classify(row: dict, countries: CountryReference, *,
             taxon_allow_list: TaxonList | None = None,
             taxon_deny_list: TaxonList | None = None) -> dict:
    if taxon_allow_list is not None and taxon_deny_list is not None:
        raise ValueError('taxon allow-list and deny-list are mutually exclusive')
    project, project_status = resolve_identity(row.get('study_accession', ''), row.get('secondary_study_accession', ''), 'project')
    sample, sample_status = resolve_identity(row.get('sample_accession', ''), row.get('secondary_sample_accession', ''), 'sample')
    platform = row.get('instrument_platform', '').strip().upper()
    layout = row.get('library_layout', '').strip().upper()
    platform_family = SUPPORTED.get(platform, ('', ''))[1] if SUPPORTED.get(platform, ('', ''))[0] == layout else ''
    strain = row.get('strain', '')
    isolate = row.get('isolate', '')
    evidence = ('MEANINGFUL_STRAIN' if meaningful(strain) and not SPECIMEN_ONLY.fullmatch(strain.strip())
                else 'MEANINGFUL_ISOLATE' if meaningful(isolate) and not SPECIMEN_ONLY.fullmatch(isolate.strip()) else '')
    text = ' '.join(row.get(k, '') for k in ('scientific_name', 'sample_title', 'study_title', 'strain', 'isolate'))
    reasons = [code for code, pattern in LABEL_PATTERNS if pattern.search(text)]
    if taxon_allow_list is not None and not taxon_allow_list.matches(row):
        reasons.append('TAXON_NOT_ALLOWED')
    if taxon_deny_list is not None and taxon_deny_list.matches(row):
        reasons.append('TAXON_EXPLICITLY_DENIED')
    contradictions = [code for code, pattern in CONTRADICTION_PATTERNS if pattern.search(text)]
    reasons += contradictions
    if project_status != 'RESOLVED': reasons.append('BIOPROJECT_' + project_status)
    if sample_status != 'RESOLVED': reasons.append('BIOSAMPLE_' + sample_status)
    if not platform_family or row.get('library_strategy') != 'WGS' or row.get('library_source') != 'GENOMIC':
        reasons.append('UNSUPPORTED_LIBRARY_PLATFORM')
    scientific_name = row.get('scientific_name', '').strip()
    instrument_model = row.get('instrument_model', '').strip()
    if not scientific_name: reasons.append('SCIENTIFIC_NAME_MISSING')
    if not instrument_model: reasons.append('INSTRUMENT_MODEL_MISSING')
    if not evidence: reasons.append('NO_DIRECT_ISOLATE_EVIDENCE')
    reasons = sorted(set(reasons))
    geo = countries.normalize(row.get('country', ''))
    result = dict(row)
    result.update(canonical_bioproject=project, bioproject_status=project_status,
                  canonical_biosample=sample, biosample_status=sample_status,
                  instrument_platform_raw=row.get('instrument_platform', ''),
                  instrument_model_raw=row.get('instrument_model', ''),
                  platform_family=platform_family, **geo,
                  isolate_grade='A' if evidence else 'U', direct_isolate_reason=evidence,
                  contradiction_codes='|'.join(sorted(set(contradictions))),
                  exclusion_codes='|'.join(reasons), primary_eligible=not reasons,
                  scientific_name_canonical=scientific_name,
                  instrument_model_canonical=instrument_model,
                  first_public_year=row.get('first_public', '')[:4],
                  genus=scientific_name.split()[0] if scientific_name else '')
    return result


DERIVED_FIELDS = ['canonical_bioproject', 'bioproject_status', 'canonical_biosample',
                  'biosample_status', 'instrument_platform_raw', 'instrument_model_raw',
                  'platform_family', 'country_raw', 'sample_country_normalized',
                  'sample_country_iso3', 'sample_region_candidate', 'sample_geo_status',
                  'sample_geo_source', 'isolate_grade', 'direct_isolate_reason',
                  'contradiction_codes', 'exclusion_codes', 'primary_eligible',
                  'scientific_name_canonical', 'instrument_model_canonical',
                  'first_public_year', 'genus']
