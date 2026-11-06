"""
Deterministic accession-driven SRA cache resolver.

There is deliberately NO glob, rglob, os.walk, or directory enumeration.

For each active registry accession, directly test only documented paths.

Built-in layouts:
    <sra_root>/<ACC>.sra
    <sra_root>/<ACC>/<ACC>.sra
    <sra_root>/aws_cache/<ACC>/<ACC>.sra
    <sra_root>/ncbi_cache/<ACC>/<ACC>.sra

Explicit extra roots:
    <root>/<ACC>.sra
    <root>/<ACC>/<ACC>.sra

Extra roots come from:
    [cache] extra_sra_roots in the active config
    WGS_CONTA_EXTRA_SRA_ROOTS
    run_master.py --add-sra-cache PATH

External roots are never recursively searched.
"""

from __future__ import annotations

import configparser
import csv
import os
import re
import sys
from pathlib import Path
from typing import Iterable


ACC_RE = re.compile(r"^(?:SRR|ERR|DRR)\d+$")
TRUE = {"yes", "true", "1", "active", "on"}


def _arg_value(flag: str, default: str = "") -> str:
    prefix = flag + "="
    argv = sys.argv[1:]

    for i, token in enumerate(argv):
        if token.startswith(prefix):
            return token[len(prefix):]

        if token == flag and i + 1 < len(argv):
            return argv[i + 1]

    return default


def _split_roots(raw: str) -> list[Path]:
    """
    Accept comma, newline, or ':' separated roots on this Linux deployment.
    """
    if not raw:
        return []

    tokens = []

    for line in raw.replace(",", "\n").splitlines():
        tokens.extend(line.split(os.pathsep))

    out = []
    seen = set()

    for token in tokens:
        token = token.strip()

        if not token:
            continue

        p = Path(token).expanduser()
        key = str(p)

        if key not in seen:
            seen.add(key)
            out.append(p)

    return out


def _config_path() -> Path:
    raw = _arg_value("--config")

    if raw:
        return Path(raw).expanduser()

    return Path(os.environ.get("WGS_CONTA_CONFIG", "config.ini"))


def _accessions_path() -> Path:
    raw = _arg_value("--accessions")

    if raw:
        return Path(raw).expanduser()

    return Path(os.environ.get("WGS_CONTA_REGISTRY", "master/accessions.tsv"))


def active_registry_accessions() -> list[str]:
    """
    Read accession identifiers only from the authoritative registry.

    We deliberately use all ACTIVE accessions rather than trying to duplicate
    master_cli's compatibility/requeue logic here. Even thousands of active
    accessions cost only a few direct path stat() calls each.
    """
    path = _accessions_path()

    if not path.exists():
        return []

    out = []
    seen = set()

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        for row in reader:
            acc = (
                row.get("Accession")
                or row.get("accession")
                or row.get("run_accession")
                or ""
            ).strip()

            if not ACC_RE.fullmatch(acc):
                continue

            raw_active = str(
                row.get("Active", "yes")
            ).strip().lower()

            active = (
                raw_active == ""
                or raw_active in TRUE
            )

            if active and acc not in seen:
                seen.add(acc)
                out.append(acc)

    return out


def configured_extra_sra_roots() -> list[Path]:
    roots = []

    cfg_path = _config_path()

    if cfg_path.exists():
        cfg = configparser.ConfigParser()
        cfg.read(cfg_path)

        if cfg.has_section("cache"):
            roots.extend(
                _split_roots(
                    cfg.get(
                        "cache",
                        "extra_sra_roots",
                        fallback="",
                    )
                )
            )

    roots.extend(
        _split_roots(
            os.environ.get(
                "WGS_CONTA_EXTRA_SRA_ROOTS",
                "",
            )
        )
    )

    out = []
    seen = set()

    for p in roots:
        key = str(p)

        if key not in seen:
            seen.add(key)
            out.append(p)

    return out


def _usable(path: Path, accession: str) -> bool:
    try:
        return (
            path.name == f"{accession}.sra"
            and path.is_file()
            and path.stat().st_size > 0
        )
    except OSError:
        return False


def candidate_paths(
    sra_root: Path,
    accession: str,
) -> list[Path]:

    root = Path(sra_root)

    paths = [
        root / f"{accession}.sra",
        root / accession / f"{accession}.sra",
        root / "aws_cache" / accession / f"{accession}.sra",
        root / "ncbi_cache" / accession / f"{accession}.sra",
    ]

    for extra in configured_extra_sra_roots():
        paths.extend([
            extra / f"{accession}.sra",
            extra / accession / f"{accession}.sra",
        ])

    # Ordered deduplication.
    out = []
    seen = set()

    for p in paths:
        key = str(p)

        if key not in seen:
            seen.add(key)
            out.append(p)

    return out


def iter_deterministic_sra_files(
    sra_root: Path,
) -> Iterable[Path]:
    """
    Drop-in iterator replacing:

        root.rglob("*.sra")

    but accession-driven and nonrecursive.
    """
    yielded = set()

    for accession in active_registry_accessions():

        for path in candidate_paths(
            sra_root,
            accession,
        ):
            if not _usable(path, accession):
                continue

            try:
                key = str(path.resolve())
            except OSError:
                key = str(path)

            if key in yielded:
                break

            yielded.add(key)
            yield path
            break
