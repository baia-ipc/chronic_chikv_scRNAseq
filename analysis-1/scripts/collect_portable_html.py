#!/usr/bin/env python3
"""
Collect portable HTML reports into a timestamped results directory.

Usage:
    collect_portable_html.py [options]

Arguments:
    None.

Options:
    -h --help                 Show this help message and exit.
    -V --version              Show this script version and exit.
    -v --verbose              Enable verbose logging.
    -q --quiet                Suppress warning output.
    --analysis-dir <path>     Analysis directory. Defaults to this script's
                              parent directory.
    --results-dir <path>      Results directory. Defaults to
                              <analysis-dir>/results.
    --prefix <prefix>         Output directory prefix
                              [default: portable_html].
    --no-index                Do not write the Markdown table to index.md.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from docopt import docopt
from loguru import logger

VERSION = "0.1.0"
LOCAL_HREF_SRC_RE = re.compile(
    r"""(?:href|src)=(["'])"""
    r"""(?P<ref>(?!data:|https?:|mailto:|#|javascript:)[^"'>]+)"""
    r"""\1""",
    re.IGNORECASE,
)


@dataclass(frozen=True)
class Report:
    source: Path
    output_name: str


def set_logger(verbose: bool, quiet: bool) -> None:
    logger.remove()
    if quiet:
        logger.add(lambda msg: None, level="DEBUG")
    elif verbose:
        logger.add(lambda msg: sys.stderr.write(msg), level="INFO")
    else:
        logger.add(lambda msg: sys.stderr.write(msg), level="WARNING")


def human_size(size: int) -> str:
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = float(size)
    for unit in units:
        if value < 1024 or unit == units[-1]:
            if unit == "B":
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024
    raise AssertionError("unreachable")


def markdown_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("|", "\\|")


def unique_output_dir(results_dir: Path, prefix: str) -> Path:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    candidate = results_dir / f"{prefix}_{stamp}"
    if not candidate.exists():
        return candidate

    index = 2
    while True:
        suffixed = results_dir / f"{prefix}_{stamp}_{index}"
        if not suffixed.exists():
            return suffixed
        index += 1


def local_href_src_refs(path: Path, limit: int = 5) -> list[str]:
    refs: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            for match in LOCAL_HREF_SRC_RE.finditer(line):
                refs.append(match.group("ref"))
                if len(refs) >= limit:
                    return refs
    return refs


def report_manifest() -> list[Report]:
    reports: list[Report] = [
        Report(
            Path("steps/03.V.integrated/rundir/UMAP_single_cell_types.portable.html"),
            "03a_UMAP_single_cell_types.html",
        ),
        Report(
            Path("steps/03.V.integrated/rundir/post_harmony_analysis.portable.html"),
            "03b_post_harmony_analysis.html",
        ),
        Report(
            Path("steps/04.V.prop_results/rundir/proportion_analysis.portable.html"),
            "04_proportion_analysis.html",
        ),
    ]

    for pair, prefix in (
        ("A", "05A"),
        ("C", "05C"),
        ("M6", "05S"),
        ("NC", "05N"),
    ):
        reports.append(
            Report(
                Path(
                    "steps/05.V.de_results/rundir/"
                    f"deseq_tables_and_vulcanos.{pair}.portable.html"
                ),
                f"{prefix}_deseq_tables_and_vulcanos.{pair}.html",
            )
        )

    reports.append(
        Report(
            Path("steps/06.V.pathways_results/rundir/viz_pathway_analysis.Reactome.GSEA.portable.html"),
            "06GR_viz_pathway_analysis.Reactome.GSEA.html",
        )
    )

    type_codes = (("circle", "C"), ("hierarchy", "H"), ("results", "R"))
    constant_codes = (("6m", "S"), ("A", "A"), ("C", "C"), ("NC", "N"))
    filter_codes = (("all", "a"), ("pfilt", "p"))
    collect_07v_all = False
    for report_type, type_code in type_codes:
        for constant, constant_code in constant_codes:
            comparison = "C_vs_NC" if constant in {"6m", "A"} else "A_vs_6m"
            for filtering, filter_code in filter_codes:
                if filtering == "all" and not collect_07v_all:
                    continue
                basename = (
                    f"cellchat_{report_type}.{constant}_{comparison}.{filtering}"
                )
                prefix = f"07{type_code}{constant_code}{filter_code}"
                reports.append(
                    Report(
                        Path(
                            "steps/07.V.cellchat_results/rundir/"
                            f"{basename}.portable.html"
                        ),
                        f"{prefix}_{basename}.html",
                    )
                )

    return reports


def link_reports(
    analysis_dir: Path,
    results_dir: Path,
    prefix: str,
    write_index: bool,
) -> int:
    reports = report_manifest()

    missing: list[Report] = []
    invalid_refs: list[tuple[Report, Path, list[str]]] = []
    for report in reports:
        source = analysis_dir / report.source
        if not source.exists():
            missing.append(report)
            continue
        refs = local_href_src_refs(source)
        if refs:
            invalid_refs.append((report, source, refs))

    if invalid_refs:
        logger.error(
            "Portable HTML collection aborted: {} report(s) contain external "
            "local href/src references.",
            len(invalid_refs),
        )
        for report, source, refs in invalid_refs:
            logger.error("Report: {} -> {}", source, report.output_name)
            for ref in refs:
                logger.error("  local href/src: {}", ref)
        logger.error(
            "Re-render these reports so all portable HTML assets are embedded "
            "as data URIs."
        )
        return 1

    output_dir = unique_output_dir(results_dir, prefix)
    output_dir.mkdir(parents=True)

    collected: list[tuple[Report, Path, os.stat_result]] = []
    total_size = 0

    for report in reports:
        source = analysis_dir / report.source
        target = output_dir / report.output_name
        if not source.exists():
            continue
        os.link(source, target)
        stat = source.stat()
        total_size += stat.st_size
        collected.append((report, target, stat))

    lines = [
        f"# Portable HTML collection: `{output_dir}`",
        "",
        f"- Collected: {len(collected)}",
        f"- Missing: {len(missing)}",
        f"- Total size: {human_size(total_size)}",
        "",
    ]

    if missing:
        lines.extend(["## Missing Portable HTML", ""])
        for report in missing:
            lines.append(f"- `{analysis_dir / report.source}`")
        lines.append("")

    lines.extend(
        [
            "## Collected Reports",
            "",
            "| File | Compiled | Size |",
            "|---|---:|---:|",
        ]
    )
    for _, target, stat in sorted(collected, key=lambda item: item[1].name):
        compiled = datetime.fromtimestamp(stat.st_mtime).strftime(
            "%Y-%m-%d %H:%M:%S"
        )
        lines.append(
            "| "
            f"`{markdown_escape(target.name)}` | "
            f"{compiled} | "
            f"{human_size(stat.st_size)} |"
        )
    lines.append(f"| **Total** |  | **{human_size(total_size)}** |")
    lines.append("")

    report_text = "\n".join(lines)
    print(report_text)

    if write_index:
        (output_dir / "index.md").write_text(report_text, encoding="utf-8")

    if missing:
        logger.warning("Missing {} expected portable HTML reports", len(missing))

    return 1 if missing else 0


if __name__ == "__main__":
    args = docopt(__doc__, version=VERSION)
    set_logger(verbose=args["--verbose"], quiet=args["--quiet"])

    if args["--analysis-dir"] is None:
        analysis_dir = Path(__file__).resolve().parents[1]
    else:
        analysis_dir = Path(args["--analysis-dir"]).resolve()
    if args["--results-dir"] is None:
        results_dir = analysis_dir / "results"
    else:
        results_dir = Path(args["--results-dir"]).resolve()

    sys.exit(
        link_reports(
            analysis_dir=analysis_dir,
            results_dir=results_dir,
            prefix=args["--prefix"],
            write_index=not args["--no-index"],
        )
    )
