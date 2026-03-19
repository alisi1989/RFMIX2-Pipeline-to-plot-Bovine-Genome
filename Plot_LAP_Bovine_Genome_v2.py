#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Alessandro Lisi & Michael C. Campbell
# Modified for bovine genome (29 chromosomes)
#
# Requires rsvg-convert installed for SVG -> PDF conversion.
# macOS: brew install librsvg
# Linux: apt install librsvg2-bin

import argparse
import xml.etree.ElementTree as ET
import os
import subprocess
from typing import Dict, List, Tuple


def parse_svg_filename(filename: str) -> str:
    """Ensure the input SVG filename ends with .svg."""
    if not filename.endswith(".svg"):
        filename += ".svg"
    return filename


def parse_output_paths(out_arg: str) -> Tuple[str, str]:
    """Return (output_svg_path, output_pdf_path) from a user-provided -O argument."""
    if out_arg.lower().endswith(".pdf"):
        base = os.path.splitext(out_arg)[0]
        return base + ".svg", out_arg
    if out_arg.lower().endswith(".svg"):
        base = os.path.splitext(out_arg)[0]
        return out_arg, base + ".pdf"
    return out_arg + ".svg", out_arg + ".pdf"


def parse_subpopulation_header(bed_file: str) -> Tuple[List[str], List[int]]:
    """
    Parse a header like:
    #Subpopulation order/codes: AF_Taurine=0 EU_Taurine=1 Indicine=2

    Returns:
        subpop_names: ordered list of subpopulation names
        subpop_codes: ordered list of subpopulation numeric codes
    """
    subpop_names = []
    subpop_codes = []

    with open(bed_file, "r") as f:
        first_line = f.readline().strip()

    if first_line.startswith("#Subpopulation order/codes:"):
        try:
            items = first_line.split(":", 1)[1].strip().split()
            for item in items:
                if "=" not in item:
                    continue
                name, code = item.split("=", 1)
                subpop_names.append(name)
                subpop_codes.append(int(code))
        except Exception:
            # If parsing fails, fall back later to generic labels
            return [], []

    return subpop_names, subpop_codes


def insert_colored_regions(
    svg_file: str,
    bed_file: str,
    output_svg_file: str,
    output_pdf_file: str,
    individual_name: str,
    subpop_names: List[str],
):
    # ---------------------------
    # Layout parameters
    # ---------------------------
    x_base = 130.900
    x_step = 63.0
    hap_width = 21.0
    top_y = 176.50
    max_svg_height = 666.300

    # ---------------------------
    # Bovine chromosome lengths (1..29)
    # ---------------------------
    chromosome_lengths = [
        158534110, 136231102, 121005158, 120000601, 120089316, 117806340, 110682743,
        113319770, 105454467, 103308737, 106982474, 87216183, 83472345, 82403003,
        85007780, 81013979, 73167244, 65820629, 63449741, 71974595, 69862954, 60773035,
        52498615, 62317253, 42350435, 51992305, 45612108, 45940150, 51098607
    ]

    # ---------------------------
    # Build chromosome coordinates
    # ---------------------------
    chromosome_coordinates = {}
    max_length = max(chromosome_lengths)
    baseline = top_y + max_svg_height

    for i, length in enumerate(chromosome_lengths, start=1):
        scaled_height = (length / max_length) * max_svg_height
        y_top = baseline - scaled_height
        x_chr = x_base + (i - 1) * x_step

        chromosome_coordinates[f"chromosome{i}_hap1"] = {
            "x": round(x_chr, 3),
            "y": round(y_top, 3),
            "width": hap_width,
            "height": round(scaled_height, 3),
        }
        chromosome_coordinates[f"chromosome{i}_hap2"] = {
            "x": round(x_chr + hap_width, 3),
            "y": round(y_top, 3),
            "width": hap_width,
            "height": round(scaled_height, 3),
        }

    # ---------------------------
    # Parse SVG
    # ---------------------------
    original_svg_tree = ET.parse(svg_file)
    original_svg_root = original_svg_tree.getroot()

    # Title
    title_element = ET.Element(
        "text",
        {"x": "950", "y": "40", "fill": "black", "font-size": "45"}
    )
    title_element.text = individual_name
    original_svg_root.append(title_element)

    rectangle_elements = []
    line_elements = []

    # Keep unique colors in the order they appear in the BED file
    observed_colors = []
    seen_colors = set()

    # ---------------------------
    # Read BED and create SVG elements
    # ---------------------------
    with open(bed_file, "r") as bed:
        for line in bed:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            chrom_raw, start_s, end_s, line_type, color, hap_s = parts[:6]

            # Save unique colors in order of first appearance
            if color not in seen_colors:
                seen_colors.add(color)
                observed_colors.append(color)

            # Normalize chromosome token (accepts "chr1" or "1")
            chrom_token = chrom_raw.strip()
            if chrom_token.lower().startswith("chr"):
                chrom_token = chrom_token[3:]

            try:
                chromosome = int(chrom_token)
                start = int(float(start_s))
                end = int(float(end_s))
                haplotype = int(hap_s)
            except Exception:
                continue

            chromosome_id = f"chromosome{chromosome}_hap{haplotype}"
            if chromosome_id not in chromosome_coordinates:
                continue

            chromosome_data = chromosome_coordinates[chromosome_id]
            chromosome_length = chromosome_lengths[chromosome - 1]

            x = chromosome_data["x"]
            width = chromosome_data["width"]

            # Map genomic coordinates onto the chromosome height
            start_y = chromosome_data["y"] + chromosome_data["height"] - (end / chromosome_length) * chromosome_data["height"]
            end_y = chromosome_data["y"] + chromosome_data["height"] - (start / chromosome_length) * chromosome_data["height"]

            if line_type == "geom_line":
                offset = 15
                line_element_start = ET.Element("line", {
                    "x1": str(x - offset),
                    "y1": str(start_y),
                    "x2": str(x + width + offset),
                    "y2": str(start_y),
                    "stroke": color,
                    "stroke-width": "2",
                    "stroke-dasharray": "10"
                })
                line_element_end = ET.Element("line", {
                    "x1": str(x - offset),
                    "y1": str(end_y),
                    "x2": str(x + width + offset),
                    "y2": str(end_y),
                    "stroke": color,
                    "stroke-width": "2",
                    "stroke-dasharray": "10"
                })
                line_elements.append(line_element_start)
                line_elements.append(line_element_end)
            else:
                rect_h = max(1.0, end_y - start_y)
                rect_element = ET.Element("rect", {
                    "x": str(x),
                    "y": str(start_y),
                    "width": str(width),
                    "height": str(rect_h),
                    "fill": color
                })
                rectangle_elements.append(rect_element)

    # Insert rectangles first (behind)
    for elem in reversed(rectangle_elements):
        original_svg_root.insert(0, elem)

    # Insert lines on top of rectangles
    for elem in reversed(line_elements):
        original_svg_root.insert(0, elem)

    # ---------------------------
    # Legend
    # ---------------------------
    # The legend uses the subpopulation names from the header and the colors
    # observed in the BED file in order of first appearance.
    if not subpop_names:
        subpop_names = [f"ancestry{i}" for i in range(len(observed_colors))]

    # If there are more names than observed colors, use fallback colors
    fallback_palette = [
        "#a32e2e", "#0a0ae0", "#bfa004", "#d18311", "#22ba9d",
        "#839dfc", "#9a5dc1", "#26962b", "#707070", "#00cfff",
        "#790ee0", "#ff4d6d", "#2d6a4f", "#f77f00", "#4ea8de"
    ]

    y_offset = 100
    legend_x_rect = 1800
    legend_x_text = 1860

    for idx, name in enumerate(subpop_names):
        if idx < len(observed_colors):
            color = observed_colors[idx]
        else:
            color = fallback_palette[idx % len(fallback_palette)]

        rect_element = ET.Element("rect", {
            "x": str(legend_x_rect),
            "y": str(y_offset),
            "width": "45",
            "height": "25",
            "fill": color
        })
        original_svg_root.append(rect_element)

        text_element = ET.Element("text", {
            "x": str(legend_x_text),
            "y": str(y_offset + 20),
            "fill": "black",
            "font-size": "26"
        })
        text_element.text = name
        original_svg_root.append(text_element)

        y_offset += 30

    # ---------------------------
    # Save SVG
    # ---------------------------
    os.makedirs(os.path.dirname(output_svg_file) or ".", exist_ok=True)
    original_svg_tree.write(output_svg_file)

    # ---------------------------
    # Convert SVG -> PDF
    # ---------------------------
    os.makedirs(os.path.dirname(output_pdf_file) or ".", exist_ok=True)
    subprocess.run(
        ["rsvg-convert", "-f", "pdf", "-o", output_pdf_file, output_svg_file],
        check=True
    )
    print(f"File converted to PDF: {output_pdf_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Insert colored regions from a BED file into an SVG file (bovine-ready)."
    )
    parser.add_argument(
        "-B",
        type=parse_svg_filename,
        default="hg38",
        help="Input SVG base file (without extension if using default)."
    )
    parser.add_argument("-I", type=str, required=True, help="Input BED file")
    parser.add_argument(
        "-O",
        type=str,
        required=True,
        help="Output path: .pdf, .svg, or basename/prefix"
    )

    args = parser.parse_args()

    individual_name = os.path.basename(args.I).split(".")[0]
    out_svg, out_pdf = parse_output_paths(args.O)

    subpop_names, subpop_codes = parse_subpopulation_header(args.I)

    insert_colored_regions(
        args.B,
        args.I,
        out_svg,
        out_pdf,
        individual_name,
        subpop_names
    )


if __name__ == "__main__":
    main()