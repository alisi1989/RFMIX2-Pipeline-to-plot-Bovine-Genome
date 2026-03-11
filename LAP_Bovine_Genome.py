#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###    Example of usage: ###
# Caso base (tutti i chr1..22)
# python rfmix2_plot_pipeline_pro.py \
# --prefix /path/to/out/prefix_ \
# --output-dir out_bed

#  Con threads e feature su chr2
#  python rfmix2_plot_pipeline_pro.py \
#  --prefix results/sample_ \
#  --chr chr1 chr2 chr3 \
#  --output-dir out_bed \
#  --threads 4 \
#  --from-bp 135000000 --to-bp 135200000 --chromosome chr2

# Solo lista operazioni, senza scrivere file
# python rfmix2_plot_pipeline_pro.py \
# --prefix results/sample_ --output-dir out_bed --dry-run

"""
RFMix2 → unified ancestry BED generator (single-file, professionalized)

What it does (end‑to‑end):
1) Discovers and combines per‑chromosome RFMix2 .msp.tsv chunks into a single MSP per individual.
2) Converts MSP into haplotype BEDs (hap1 / hap2) with ancestry labels.
3) Merges hap BEDs into a final BED with rendering hints (geom_rect + colors),
   plus optional feature highlight lines (geom_line).

Key improvements over the original:
- Robust file discovery and individual ID parsing (regex‑based, supports chr and no‑chr tokens).
- Safer I/O (streaming copy for huge files, avoids readlines() on big MSPs).
- Pandas CSV parsing with comment handling and typed columns.
- Natural sorting by chromosome and position.
- Optional parallel processing across individuals (--threads).
- Proper logging (--debug / --quiet) instead of prints.
- Options: --require-all-chroms, --keep-temp, --dry-run, color config file (JSON).
- Resilient ancestry mapping: default to "ancestry{int}", or load from JSON mapping.

Supports up to ancestry0..14 (15 total ancestries) by default.

Dependencies: Python ≥3.8, pandas

Author: Alessandro Lisi
"""



from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

# -----------------------------------------------------------------------------
# Logging
# -----------------------------------------------------------------------------
LOGGER = logging.getLogger("rfmix2")


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
CHROM_TOKEN_RE = re.compile(r"^(?:chr)?(\w+)$", re.IGNORECASE)


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def strip_chr(token: str) -> str:
    """Remove 'chr' prefix if present (case‑insensitive)."""
    m = CHROM_TOKEN_RE.match(str(token))
    return m.group(1) if m else str(token)


def chrom_sort_key(chrom: str) -> Tuple[int, str]:
    """Natural chromosome sort: 1..22,X,Y,MT (case insensitive)."""
    c = strip_chr(chrom).upper()
    if c in {"X", "XX"}:
        return (1000, "X")
    if c in {"Y"}:
        return (1001, "Y")
    if c in {"MT", "M"}:
        return (1002, "MT")
    try:
        return (int(c), "")
    except ValueError:
        # fallback: non‑numeric contigs at the end, but keep stable order
        return (2000, c)


@dataclass
class Feature:
    chrom: str
    start_bp: int
    end_bp: int


# -----------------------------------------------------------------------------
# Discovery & Combination
# -----------------------------------------------------------------------------
MSP_SUFFIX = ".msp.tsv"




def _extract_individual_from_basename(basename: str, chrom: str) -> str:
    """Infer the 'individual' stem from a basename like '<prefix>_<sample>_<chr>.msp.tsv'.

    Strategy: remove the trailing '<sep>chrom' token and the suffix; then rstrip separators.
    This preserves the original prefix+sample stem, which is what the original script used
    as the individual key.
    """
    name = basename
    if name.endswith(MSP_SUFFIX):
        name = name[: -len(MSP_SUFFIX)]
    # Remove a single trailing chrom token preceded by optional separators
    # e.g., 'foo_bar_chr1' -> 'foo_bar', 'foo1chr1' -> 'foo1'
    chrom_escaped = re.escape(chrom)
    name = re.sub(rf"[._-]?{chrom_escaped}$", "", name)
    return name.rstrip("._-")


# ---------------------------
# Helper: scan small portion of an MSP to collect chrom tokens
# ---------------------------
def _scan_chroms_in_msp(path: Path, max_lines: int = 10000) -> set:
    """Scan up to `max_lines` non-comment lines of an MSP and return set of normalized chrom tokens found.
    Normalizes tokens via strip_chr() and uppercases (so 'chr1' and '1' map to '1')."""
    found = set()
    try:
        with path.open("r", encoding="utf-8", errors="ignore") as fh:
            count = 0
            for line in fh:
                if count >= max_lines:
                    break
                if not line or line.startswith("#"):
                    continue
                toks = line.rstrip("\n").split("\t")
                if not toks:
                    continue
                ch = toks[0]
                if ch is None:
                    continue
                ch_norm = strip_chr(ch).upper()
                found.add(ch_norm)
                count += 1
    except Exception:
        # on any read error, return empty set (safer to assume unknown)
        return set()
    return found


# ---------------------------
# Revised discover_inputs that accepts a single combined MSP file or the original glob behavior
# ---------------------------
def discover_inputs(prefix: str, chroms: Sequence[str]) -> Tuple[Dict[str, List[Path]], List[str]]:
    """Discover .msp.tsv files grouped by individual stem.

    Behavior:
    - If `prefix` points to an existing file ending with .msp.tsv, treat it as a single individual
      (use basename without suffix as the individual key).
    - Otherwise do the previous glob-based discovery of per-chromosome MSP files.
    """
    individuals: Dict[str, List[Path]] = {}

    prefix_path = Path(prefix)

    # Normalize requested chrom set (e.g., {'1','2',...,'22','X'})
    req_norm = {strip_chr(c).upper() for c in chroms}

    # If prefix is an explicit file, use it as a single individual (supports multi-chrom MSP)
    if prefix_path.is_file() and str(prefix_path).lower().endswith(MSP_SUFFIX):
        ind_key = prefix_path.with_suffix("").name  # basename without .msp.tsv
        individuals[ind_key] = [prefix_path]
        return individuals, [ind_key]

    # Otherwise fallback to the globbing behavior (same as before)
    search_dir = prefix_path.parent if str(prefix_path.parent) not in ("", ".") else Path(".")
    prefix_base = prefix_path.name

    all_hits = sorted(search_dir.glob(f"{prefix_base}*{MSP_SUFFIX}"))
    if not all_hits:
        return {}, []

    for fp in all_hits:
        if not fp.is_file():
            continue
        chrom_tok = _guess_chrom_from_filename(fp.name)
        chrom_norm = strip_chr(chrom_tok).upper()
        # If filename token doesn't look like a chrom (e.g., the file is a combined MSP with no trailing chrom),
        # we still accept it (it may be a multi-chrom MSP file). We'll include it only if scanning reveals requested chroms,
        # otherwise we skip it.
        if chrom_norm not in req_norm:
            # attempt to detect if this file is a multi-chrom MSP containing requested chroms
            scanned = _scan_chroms_in_msp(fp, max_lines=5000)
            if not scanned.intersection(req_norm):
                # no requested chroms found in the first chunk -> skip
                continue

        ind = _extract_individual_from_basename(fp.name, chrom_tok)
        individuals.setdefault(ind, []).append(fp)

    # Sort each individual's files by chrom natural order then by name
    for ind, files in list(individuals.items()):
        individuals[ind] = sorted(files, key=lambda p: (chrom_sort_key(_guess_chrom_from_filename(p.name)), p.name))

    unique_inds = sorted(individuals.keys())
    return individuals, unique_inds
    
    
    
def filter_headers_for_bed(headers: List[str]) -> List[str]:
    """Keep only informative comment headers for BED outputs.

    We drop the '#chm\t...' header line because in multi-sample MSPs it includes all sample columns
    and is misleading in per-sample BED outputs.
    """
    out: List[str] = []
    for h in headers:
        if h.startswith("#chm"):
            continue
        out.append(h)
    return out


def _guess_chrom_from_filename(basename: str) -> str:
    """Best effort to recover the chrom token from a filename that ends with '<chrom>.msp.tsv'."""
    core = basename
    if core.endswith(MSP_SUFFIX):
        core = core[: -len(MSP_SUFFIX)]
    # Grab last token after separators and hope it's the chrom (works for common RFMix2 names)
    last = re.split(r"[._-]", core)[-1]
    # normalize to chrX if it's purely numeric or already has chr
    if CHROM_TOKEN_RE.match(last):
        return last if last.lower().startswith("chr") else f"chr{last}"
    return last


# -----------------------------------------------------------------------------
# MSP → hap BEDs
# -----------------------------------------------------------------------------
BASE_MSP_COLS = [
    "chm", "spos", "epos", "sgpos", "egpos", "n_snps",
]


def parse_msp_columns_from_header(path: Path) -> List[str]:
    """Return column names from the '#chm\t...' header line.

    RFMix2 MSPs often include the column header as a commented line beginning with '#chm'.
    We normalize spaces to underscores (e.g., 'n snps' -> 'n_snps').
    """
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#chm"):
                cols = line.lstrip("#").rstrip("\n").split("\t")
                return [c.strip().replace(" ", "_") for c in cols]
    # Fallback: legacy single-individual format (no explicit header line)
    return BASE_MSP_COLS + ["ind1", "ind2"]


_HAP_COL_RE = re.compile(r"^(?P<sample>.+)\.(?P<hap>[01])$")


def infer_sample_hap_pairs(columns: Sequence[str]) -> Dict[str, Tuple[str, str]]:
    """Infer diploid sample -> (hap0_col, hap1_col) from columns.

    For multi-sample MSPs, hap columns are typically named like '<sample>.0' and '<sample>.1'.
    """
    pairs: Dict[str, Dict[str, str]] = {}
    for c in columns:
        m = _HAP_COL_RE.match(str(c))
        if not m:
            continue
        sample = m.group("sample")
        hap = m.group("hap")
        pairs.setdefault(sample, {})[hap] = str(c)

    out: Dict[str, Tuple[str, str]] = {}
    for sample, d in pairs.items():
        if "0" in d and "1" in d:
            out[sample] = (d["0"], d["1"])
    return out


def safe_sample_name(sample: str) -> str:
    """Make a sample name safe for filesystem paths."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", sample)


def read_msp(path: Path, columns: Sequence[str]) -> pd.DataFrame:
    """Load an MSP TSV with comment headers ('#') into a typed DataFrame.

    Supports both legacy single-individual MSPs (ind1/ind2) and multi-sample MSPs with many hap columns.
    """
    cols = list(columns)

    # Build dtype map: first 6 columns typed, remaining columns as nullable integers (ancestry codes)
    dtype_map: Dict[str, object] = {
        "chm": str,
        "spos": float,
        "epos": float,
        "sgpos": float,
        "egpos": float,
        "n_snps": int,
    }
    for c in cols[6:]:
        # ancestry codes; allow missing
        dtype_map[str(c)] = "Int16"

    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=cols,
        dtype=dtype_map,
        engine="c",
    )

    # Cast positions to integers if they are integral floats
    for col in ("spos", "epos"):
        df[col] = df[col].astype(float).round().astype(int)
    return df


def load_headers_from_msp(path: Path) -> List[str]:
    headers: List[str] = []
    with path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                headers.append(line.rstrip("\n"))
            else:
                break
    return headers


def msp_to_hap_beds(
    msp_path: Path,
    out_prefix: Path,
    ancestry_map: Dict[int, str],
    keep_headers: Optional[List[str]] = None,
) -> Dict[str, Tuple[Path, Path]]:
    """Convert an MSP to per-sample hap BEDs.

    Returns:
        Dict[sample] -> (hap1_bed_path, hap2_bed_path)

    Supports:
    - Legacy single-individual MSPs with columns ind1/ind2
    - Multi-sample MSPs with hap columns like '<sample>.0' and '<sample>.1'
    """
    cols = parse_msp_columns_from_header(msp_path)
    df = read_msp(msp_path, cols)

    # Convert ancestry integer codes to labels; if not numeric leave as is
    def map_series(s: pd.Series) -> pd.Series:
        codes = pd.to_numeric(s, errors="coerce")
        labels = codes.map(lambda x: ancestry_map.get(int(x), f"ancestry{int(x)}") if pd.notna(x) else s)
        out = labels.astype(object)
        mask_na = codes.isna()
        out.loc[mask_na] = s.loc[mask_na]
        return out

    results: Dict[str, Tuple[Path, Path]] = {}

    # Multi-sample mode: infer sample pairs from column names
    pairs = infer_sample_hap_pairs(cols)
    if pairs:
        for sample, (c0, c1) in sorted(pairs.items()):
            sample_safe = safe_sample_name(sample)
            # hap1 corresponds to .0, hap2 to .1 (consistent within our toolchain)
            hap_outputs: List[Tuple[str, str]] = [("hap1", c0), ("hap2", c1)]
            created: List[Path] = []

            for hap_name, col in hap_outputs:
                sub = df[["chm", "spos", "epos", col]].copy()
                sub[col] = map_series(sub[col])
                sub.columns = ["chm", "spos", "epos", "ancestry"]

                out_path = out_prefix.with_name(out_prefix.name + f"_{sample_safe}_{hap_name}.bed")
                with out_path.open("w", encoding="utf-8") as out:
                    if keep_headers:
                        out.write("\n".join(keep_headers) + "\n")
                    sub.to_csv(
                        out,
                        sep="\t",
                        header=False,
                        index=False,
                        quoting=csv.QUOTE_NONE,
                        escapechar="\\",
                        lineterminator="\n",
                    )
                LOGGER.info("Wrote %s", out_path)
                created.append(out_path)

            results[sample] = (created[0], created[1])

        return results

    # Legacy mode: expect ind1/ind2
    if "ind1" not in df.columns or "ind2" not in df.columns:
        raise ValueError(f"MSP appears to have no sample hap columns and is missing ind1/ind2: {msp_path}")

    hap_outputs_legacy: List[Tuple[str, str]] = [("hap1", "ind1"), ("hap2", "ind2")]
    created_legacy: List[Path] = []
    for hap_name, col in hap_outputs_legacy:
        sub = df[["chm", "spos", "epos", col]].copy()
        sub[col] = map_series(sub[col])
        sub.columns = ["chm", "spos", "epos", "ancestry"]

        out_path = out_prefix.with_name(out_prefix.name + f"_{hap_name}.bed")
        with out_path.open("w", encoding="utf-8") as out:
            if keep_headers:
                out.write("\n".join(keep_headers) + "\n")
            sub.to_csv(
                out,
                sep="\t",
                header=False,
                index=False,
                quoting=csv.QUOTE_NONE,
                escapechar="\\",
                lineterminator="\n",
            )
        LOGGER.info("Wrote %s", out_path)
        created_legacy.append(out_path)

    results["individual"] = (created_legacy[0], created_legacy[1])
    return results


# -----------------------------------------------------------------------------
# hap BEDs → final BED with colors & feature
# -----------------------------------------------------------------------------

def remove_chr_prefix(chrom: str) -> str:
    return strip_chr(str(chrom))


def _iter_bed_lines(bed_path: Path, geom_value: str) -> Iterable[Tuple[str, int, int, str]]:
    with bed_path.open("r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            toks = line.rstrip("\n").split("\t")
            if len(toks) < 4:
                continue
            chrom, s, e, ancestry = toks[0], toks[1], toks[2], toks[3]
            try:
                start = int(float(s))
                end = int(float(e))
            except ValueError:
                continue
            yield (remove_chr_prefix(chrom), start, end, ancestry)


def _final_line(chrom: str, start: int, end: int, feature: str, color: str, geom_value: str) -> str:
    return f"{chrom}\t{start}\t{end}\t{feature}\t{color}\t{geom_value}"


def build_color_map(args: argparse.Namespace, ancestry_labels: Iterable[str]) -> Dict[str, str]:
    # 1) From JSON file if provided
    if args.color_config:
        cfg_path = Path(args.color_config)
        with cfg_path.open("r", encoding="utf-8") as fh:
            cfg = json.load(fh)
        # Expect {"ancestry0": "#aabbcc", ...}
        return {str(k): str(v) for k, v in cfg.items()}

    # 2) From CLI defaults (ancestry0..14)
    colors: Dict[str, str] = {}
    for i in range(0, 15):
        colors[f"ancestry{i}"] = getattr(args, f"ancestry{i}")
    # Ensure any unseen labels at runtime still get some color (unknown)
    for lab in ancestry_labels:
        colors.setdefault(lab, args.unknown)
    return colors


def process_haps_to_final(
    hap1_path: Path,
    hap2_path: Path,
    colors: Dict[str, str],
    unknown_color: str,
    feature: Optional[Feature],
    final_output_path: Path,
    header_from: Optional[List[str]] = None,
) -> None:
    lines: List[str] = []

    for bed_path, gv in ((hap1_path, "1"), (hap2_path, "2")):
        for chrom, start, end, ancestry in _iter_bed_lines(bed_path, gv):
            color = colors.get(str(ancestry), unknown_color)
            lines.append(_final_line(chrom, start, end, "geom_rect", color, gv))

    if feature is not None:
        fchrom = remove_chr_prefix(feature.chrom)
        f1 = _final_line(fchrom, feature.start_bp, feature.end_bp, "geom_line", "#000000", "1")
        f2 = _final_line(fchrom, feature.start_bp, feature.end_bp, "geom_line", "#000000", "2")
        lines.extend([f1, f2])

    # Sort by chrom then start
    def sort_key(line: str) -> Tuple[Tuple[int, str], int]:
        toks = line.split("\t")
        ckey = chrom_sort_key(toks[0])
        spos = int(toks[1])
        return (ckey, spos)

    lines.sort(key=sort_key)

    with final_output_path.open("w", encoding="utf-8") as out:
        if header_from:
            out.write("\n".join(header_from) + "\n")
        out.write("\n".join(lines) + "\n")
    LOGGER.info("Wrote %s", final_output_path)


# -----------------------------------------------------------------------------
# Combine per-chrom MSPs → one MSP per individual
# -----------------------------------------------------------------------------

def combine_msp_for_individual(files: Sequence[Path], dest_path: Path, require_headers_from_first: bool = True) -> List[str]:
    """Concatenate MSPs skipping comment headers after the first file.

    Returns the header lines captured from the first file.
    """
    headers: List[str] = []
    wrote_header = False

    with dest_path.open("w", encoding="utf-8") as out:
        for i, fp in enumerate(files):
            LOGGER.debug("Concatenating %s (%d/%d)", fp, i + 1, len(files))
            with fp.open("r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if line.startswith("#"):
                        if not wrote_header:
                            headers.append(line.rstrip("\n"))
                            out.write(line)
                        # skip comment lines for subsequent files
                        if wrote_header:
                            continue
                    out.write(line)
            if not wrote_header:
                wrote_header = True
    return headers


# -----------------------------------------------------------------------------
# Orchestration per individual
# -----------------------------------------------------------------------------

def process_individual(
    individual: str,
    files: Sequence[Path],
    out_dir: Path,
    ancestry_map: Dict[int, str],
    colors: Dict[str, str],
    unknown_color: str,
    feature: Optional[Feature],
    keep_temp: bool,
) -> List[Tuple[str, Path]]:
    LOGGER.info("Processing individual: %s", individual)
    combined_path = out_dir / f"{individual}{MSP_SUFFIX}"

    # Combine MSPs
    headers = combine_msp_for_individual(files, combined_path)
    headers_bed = filter_headers_for_bed(headers)
    # Convert to hap beds (may yield many samples)
    out_prefix = combined_path.with_suffix("")
    sample_haps = msp_to_hap_beds(combined_path, out_prefix, ancestry_map, keep_headers=headers_bed)

    results: List[Tuple[str, Path]] = []
    temp_paths: List[Path] = [combined_path]

    for sample, (hap1, hap2) in sample_haps.items():
        # Final BED per sample
        sample_safe = safe_sample_name(sample)
        final_bed = out_dir / f"{sample_safe}.bed"
        process_haps_to_final(hap1, hap2, colors, unknown_color, feature, final_bed, header_from=headers_bed)
        results.append((sample, final_bed))
        temp_paths.extend([hap1, hap2])

    # Cleanup
    if not keep_temp:
        for p in temp_paths:
            try:
                p.unlink()
                LOGGER.debug("Removed %s", p)
            except Exception as e:
                LOGGER.warning("Could not remove %s: %s", p, e)

    return results


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine RFMix2 MSP chunks into final ancestry BEDs with colors and optional feature lines."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Discovery / I/O
    parser.add_argument("--prefix", required=True, help="Prefix for RFMix2 MSP files (may include path).")
    parser.add_argument(
        "--chr",
        nargs="+",
        default=[f"chr{i}" for i in range(1, 30)],
        help="List of chromosome tokens present in filenames (e.g., chr1 chr2 ... OR 1 2 ...)",
    )
    parser.add_argument("--output-dir", required=True, help="Output directory.")

    # Behavior toggles
    parser.add_argument("--require-all-chroms", action="store_true", help="Skip individuals missing any requested chromosome.")
    parser.add_argument("--keep-temp", action="store_true", help="Keep combined MSP and hap BEDs.")
    parser.add_argument("--dry-run", action="store_true", help="Only list planned operations, do not process.")
    parser.add_argument("--threads", type=int, default=1, help="Parallel individuals (I/O-bound, ThreadPool).")

    # Colors & mapping
    parser.add_argument("--color-config", type=str, help="JSON file mapping ancestry labels to colors.")
    parser.add_argument("--unknown", default="#808080", help="Color for unknown/unspecified ancestries.")
    for i, default in enumerate([
        "#a32e2e", "#0a0ae0", "#bfa004", "#d18311", "#22ba9d",
        "#839dfc", "#9a5dc1", "#26962b", "#707070", "#00cfff", "#790ee0",
        "#ff4d6d", "#2d6a4f", "#f77f00", "#4ea8de",
    ]):
        parser.add_argument(f"--ancestry{i}", default=default, help=f"Color for ancestry{i}")

    # Optional feature
    parser.add_argument("--from-bp", type=int, dest="from_bp", help="Feature start (bp)")
    parser.add_argument("--to-bp", type=int, dest="to_bp", help="Feature end (bp)")
    parser.add_argument("--chromosome", type=str, help="Chromosome for feature (e.g., chr2 or 2)")

    # Logging
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--debug", action="store_true", help="Verbose debug logs.")
    g.add_argument("--quiet", action="store_true", help="Only warnings and errors.")

    return parser.parse_args(argv)


def setup_logging(args: argparse.Namespace) -> None:
    level = logging.INFO
    if args.debug:
        level = logging.DEBUG
    elif args.quiet:
        level = logging.WARNING
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")


def build_ancestry_map() -> Dict[int, str]:
    # By default, labels are 'ancestry{i}' — mapping here is trivial, but present for clarity/extensibility.
    return {i: f"ancestry{i}" for i in range(0, 128)}  # generous upper bound


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    setup_logging(args)

    out_dir = Path(args.output_dir)
    ensure_dir(out_dir)

    # Normalize chrom tokens to include both '1' and 'chr1' variants when searching
    chroms: List[str] = []
    for tok in args.chr:
        tok = str(tok)
        if tok.lower().startswith("chr"):
            chroms.append(tok)
            chroms.append(strip_chr(tok))
        else:
            chroms.append(tok)
            chroms.append(f"chr{tok}")
    # Deduplicate while preserving order
    seen = set()
    chroms = [x for x in chroms if not (x in seen or seen.add(x))]

    individuals_map, individuals = discover_inputs(args.prefix, chroms)
    if not individuals:
        LOGGER.error("No individuals found for prefix '%s' and chromosomes %s", args.prefix, args.chr)
        return 2

    # Optionally filter individuals missing requested chroms
    if args.require_all_chroms:
        req_set = {c for c in chroms if str(c).lower().startswith("chr")}
        def has_all(files: Sequence[Path]) -> bool:
            present = { _guess_chrom_from_filename(f.name).lower() for f in files }
            return req_set.issubset(present)
        filtered = {ind: fs for ind, fs in individuals_map.items() if has_all(fs)}
        missing = sorted(set(individuals_map) - set(filtered))
        if missing:
            LOGGER.warning("Skipping %d individuals missing chromosomes: %s", len(missing), ", ".join(missing))
        individuals_map = filtered
        individuals = sorted(individuals_map)
        if not individuals:
            LOGGER.error("After filtering, no individuals have all requested chromosomes.")
            return 3

    LOGGER.info("Individuals to process: %s", ", ".join(individuals))

    # Dry run: list planned operations and exit
    if args.dry_run:
        for ind in individuals:
            files = individuals_map[ind]
            LOGGER.info("[DRY‑RUN] %s: %d MSP chunks → (per-sample BEDs in %s)", ind, len(files), out_dir)
        return 0

    ancestry_map = build_ancestry_map()

    # Build color map using discovered labels (pre‑populate with ancestry0..14)
    discovered_labels = {f"ancestry{i}" for i in range(0, 15)}
    colors = build_color_map(args, discovered_labels)

    feature: Optional[Feature] = None
    if args.from_bp is not None and args.to_bp is not None and args.chromosome is not None:
        feature = Feature(chrom=args.chromosome, start_bp=args.from_bp, end_bp=args.to_bp)

    results: List[Tuple[str, Path]] = []

    if args.threads and args.threads > 1:
        with ThreadPoolExecutor(max_workers=args.threads) as ex:
            futs = {
                ex.submit(
                    process_individual,
                    ind,
                    individuals_map[ind],
                    out_dir,
                    ancestry_map,
                    colors,
                    args.unknown,
                    feature,
                    args.keep_temp,
                ): ind
                for ind in individuals
            }
            for fut in as_completed(futs):
                try:
                    res_list = fut.result()
                    results.extend(res_list)
                except Exception as e:
                    LOGGER.error("Failed individual %s: %s", futs[fut], e)
                    return 4
    else:
        for ind in individuals:
            res = process_individual(
                ind,
                individuals_map[ind],
                out_dir,
                ancestry_map,
                colors,
                args.unknown,
                feature,
                args.keep_temp,
            )
            results.extend(res)

    LOGGER.info("Done. Generated %d BED files in %s", len(results), out_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
