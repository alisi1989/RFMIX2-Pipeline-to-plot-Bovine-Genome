#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Authors: Alessandro Lisi & Michael C. Campbell
# Modificato per supportare genoma bovino (29 cromosomi) — Alessandro Lisi (assist)
#
# Nota: assicurati di avere rsvg-convert installato per la conversione da SVG a PDF.
# macOS: brew install librsvg
# Linux: apt install librsvg2-bin  (o equivalente)

import argparse
import xml.etree.ElementTree as ET
import os
import subprocess
from typing import Dict

def insert_colored_regions(svg_file: str, bed_file: str, output_svg_file: str, output_pdf_file: str, individual_name: str, ancestries: Dict[str, str]):

    # ---------------------------
    # Parametri di layout (tweak in base al tuo SVG con Illustrator)
    # ---------------------------
    # Questi valori definiscono la posizione orizzontale di partenza (x_base),
    # lo spazio orizzontale tra coppie di cromosomi (x_step), la larghezza di
    # ogni blocco hap (hap_width), la y di partenza top_y e l'altezza massima
    # nella quale viene scalato il cromosoma più lungo (max_svg_height).
    #
    # Se aggiusti l'SVG in Illustrator, qui puoi cambiare x_base/x_step/top_y
    # per far combaciare i rettangoli sui tuoi cromosomi nell'SVG.
    x_base = 130.900      # x coordinate per chromosome1_hap1 (puoi adattare)
    x_step = 63.0         # passo orizzontale tra successive coppie chr
    hap_width = 21.0      # larghezza di ciascuna hap block (mantieni se vuoi)
    top_y = 176.50        # y superiore comune
    max_svg_height = 666.300  # altezza alla quale mappare il cromosoma più lungo

    # ---------------------------
    # Lunghezze cromosomiche bovine (1..29) — valori forniti dall'utente
    # ---------------------------
    chromosome_lengths = [
        158534110, 136231102, 121005158, 120000601, 120089316, 117806340, 110682743,
        113319770, 105454467, 103308737, 106982474, 87216183, 83472345, 82403003,
        85007780, 81013979, 73167244, 65820629, 63449741, 71974595, 69862954, 60773035,
        52498615, 62317253, 42350435, 51992305, 45612108, 45940150, 51098607
    ]

    # ---------------------------
    # Costruisci dinamicamente le coordinate (chromosome{i}_hap1 / hap2)
    # Allineamento bottom-anchored: tutti i cromosomi condividono lo stesso baseline
    # (top_y + max_svg_height). Per ogni cromosoma posizioniamo il top del box
    # come baseline - scaled_height così i cromosomi più corti restano "a terra".
    # ---------------------------
    chromosome_coordinates = {}
    max_length = max(chromosome_lengths)
    baseline = top_y + max_svg_height
    for i, length in enumerate(chromosome_lengths, start=1):
        scaled_height = (length / max_length) * max_svg_height
        # y_top in modo che il fondo (baseline) sia comune
        y_top = baseline - scaled_height
        x_chr = x_base + (i - 1) * x_step
        chromosome_coordinates[f'chromosome{i}_hap1'] = {
            'x': round(x_chr, 3),
            'y': round(y_top, 3),
            'width': hap_width,
            'height': round(scaled_height, 3),
        }
        chromosome_coordinates[f'chromosome{i}_hap2'] = {
            'x': round(x_chr + hap_width, 3),
            'y': round(y_top, 3),
            'width': hap_width,
            'height': round(scaled_height, 3),
        }

    # ---------------------------
    # Parsing e inserimento
    # ---------------------------
    original_svg_tree = ET.parse(svg_file)
    original_svg_root = original_svg_tree.getroot()

    # Aggiungi il nome dell'individuo come titolo
    title_element = ET.Element('text', {'x': "950", 'y': "40", 'fill': "black", 'font-size': "45"})
    title_element.text = individual_name
    original_svg_root.insert(0, title_element)

    rectangle_elements = []
    line_elements = []
    found_ancestries = set()

    # Leggi il BED e crea gli elementi SVG
    with open(bed_file, 'r') as bed:
        for line in bed:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split('\t')
            if len(parts) < 6:
                continue

            chrom_raw, start_s, end_s, line_type, color, hap_s = parts[:6]

            # Normalizza il token cromosoma (accetta "chr1" o "1")
            chrom_token = chrom_raw.strip()
            if chrom_token.lower().startswith("chr"):
                chrom_token = chrom_token[3:]
            # tentiamo cast a int — se fallisce saltiamo
            try:
                chromosome = int(chrom_token)
            except ValueError:
                # se non è numerico, ignora (puoi estendere per contig named)
                continue

            try:
                start = int(float(start_s))
                end = int(float(end_s))
                haplotype = int(hap_s)
            except Exception:
                continue

            chromosome_id = f'chromosome{chromosome}_hap{haplotype}'
            if chromosome_id not in chromosome_coordinates:
                # il cromosoma non è nella mappa (ad es. contig diverso): skip
                continue

            chromosome_data = chromosome_coordinates[chromosome_id]
            chromosome_length = chromosome_lengths[chromosome - 1]

            x = chromosome_data['x']
            width = chromosome_data['width']
            # calcola y relativo proporzionale alla lunghezza reale
            start_y = chromosome_data['y'] + chromosome_data['height'] - (end / chromosome_length) * chromosome_data['height']
            end_y = chromosome_data['y'] + chromosome_data['height'] - (start / chromosome_length) * chromosome_data['height']
            if line_type == 'geom_line':
                offset = 15
                line_element_start = ET.Element('line', {
                    'x1': str(x - offset),
                    'y1': str(start_y),
                    'x2': str(x + width + offset),
                    'y2': str(start_y),
                    'stroke': color,
                    'stroke-width': '2',
                    'stroke-dasharray': '10'
                })
                line_element_end = ET.Element('line', {
                    'x1': str(x - offset),
                    'y1': str(end_y),
                    'x2': str(x + width + offset),
                    'y2': str(end_y),
                    'stroke': color,
                    'stroke-width': '2',
                    'stroke-dasharray': '10'
                })
                original_svg_root.insert(0, line_element_start)
                original_svg_root.insert(0, line_element_end)
                found_ancestries.add(("geom_line", color))
            else:
                # geom_rect
                # assicurati che height sia positivo
                rect_h = max(1.0, end_y - start_y)
                rect_element = ET.Element('rect', {
                    'x': str(x),
                    'y': str(start_y),
                    'width': str(width),
                    'height': str(rect_h),
                    'fill': color
                })
                rectangle_elements.append(rect_element)

                # registra ancestry matching colore se esiste nella mappa ancestries
                for ancestry_name, ancestry_color in ancestries.items():
                    if color == ancestry_color:
                        found_ancestries.add((ancestry_name, ancestry_color))

    # inserisci i rettangoli (dietro titoli/linee)
    for elem in reversed(rectangle_elements):
        original_svg_root.insert(0, elem)

    # Crea la legenda con le ancestries trovate
    if found_ancestries:
        y_offset = 100
        legend_x_rect = 1800
        legend_x_text = 1860
        for name, color in sorted(found_ancestries):
            # Se è geom_line registrato come tuple ("geom_line", color) gestiscilo diversamente
            if name == "geom_line":
                # linea di legenda
                line_leg = ET.Element('line', {
                    'x1': str(legend_x_rect),
                    'y1': str(y_offset + 8),
                    'x2': str(legend_x_rect + 25),
                    'y2': str(y_offset + 8),
                    'stroke': color,
                    'stroke-width': '2',
                    'stroke-dasharray': '4'
                })
                original_svg_root.insert(0, line_leg)
                text_el = ET.Element('text', {'x': str(legend_x_text), 'y': str(y_offset + 12), 'fill': "black", 'font-size': "16"})
                text_el.text = "feature"
                original_svg_root.insert(0, text_el)
                y_offset += 30
                continue

            # rettangolo legenda
            rect_element = ET.Element('rect', {'x': str(legend_x_rect), 'y': str(y_offset), 'width': "45", 'height': "25", 'fill': color})
            original_svg_root.insert(0, rect_element)

            text_element = ET.Element('text', {'x': str(legend_x_text), 'y': str(y_offset + 20), 'fill': "black", 'font-size': "26"})
            text_element.text = name
            original_svg_root.insert(0, text_element)

            y_offset += 30

    # Salva SVG modificato
    os.makedirs(os.path.dirname(output_svg_file) or '.', exist_ok=True)
    original_svg_tree.write(output_svg_file)

    # Converti SVG -> PDF
    os.makedirs(os.path.dirname(output_pdf_file) or '.', exist_ok=True)
    subprocess.run(["rsvg-convert", "-f", "pdf", "-o", output_pdf_file, output_svg_file], check=True)
    print(f"File convertito in PDF: {output_pdf_file}")


def parse_svg_filename(filename: str) -> str:
    if not filename.endswith('.svg'):
        filename += '.svg'
    return filename


def parse_output_paths(out_arg: str):
    """Return (output_svg_path, output_pdf_path) from a user-provided -O argument."""
    if out_arg.lower().endswith('.pdf'):
        base = os.path.splitext(out_arg)[0]
        return base + '.svg', out_arg
    if out_arg.lower().endswith('.svg'):
        base = os.path.splitext(out_arg)[0]
        return out_arg, base + '.pdf'
    return out_arg + '.svg', out_arg + '.pdf'


def main():
    parser = argparse.ArgumentParser(description='Insert colored regions from a BED file into an SVG file (bovine-ready).')
    parser.add_argument('-B', type=parse_svg_filename, default='hg38', help='Input SVG base file (senza estensione se usi il default).')
    parser.add_argument('-I', type=str, required=True, help='Input BED file')
    parser.add_argument('-O', type=str, required=True, help='Output path: .pdf, .svg, or basename/prefix')

    args = parser.parse_args()

    # Default colors for ancestries (ancestry0..14)
    default_ancestry_colors = [
        ("ancestry0", "#a32e2e"),
        ("ancestry1", "#0a0ae0"),
        ("ancestry2", "#bfa004"),
        ("ancestry3", "#d18311"),
        ("ancestry4", "#22ba9d"),
        ("ancestry5", "#839dfc"),
        ("ancestry6", "#9a5dc1"),
        ("ancestry7", "#26962b"),
        ("ancestry8", "#707070"),
        ("ancestry9", "#00cfff"),
        ("ancestry10", "#790ee0"),
        ("ancestry11", "#ff4d6d"),
        ("ancestry12", "#2d6a4f"),
        ("ancestry13", "#f77f00"),
        ("ancestry14", "#4ea8de"),
    ]

    # Prova a leggere una possibile riga header che definisce i nomi di subpop
    ancestries = {name: color for name, color in default_ancestry_colors}
    with open(args.I, 'r') as bed_file:
        first_line = bed_file.readline().strip()
        if first_line.startswith("#Subpopulation order/codes:"):
            # formato: "#Subpopulation order/codes: NAME1=0 NAME2=1 ..."
            try:
                ancestry_names = first_line.split(":", 1)[1].strip().split()
                ancestry_map = {}
                for item in ancestry_names:
                    if "=" not in item:
                        continue
                    name, code = item.split("=")
                    ancestry_map[int(code)] = name
                ancestries = {ancestry_map.get(i, f'ancestry{i}'): color for i, (name, color) in enumerate(default_ancestry_colors)}
            except Exception:
                # fallback ai default colors
                ancestries = {name: color for name, color in default_ancestry_colors}

    individual_name = os.path.basename(args.I).split('.')[0]
    out_svg, out_pdf = parse_output_paths(args.O)
    insert_colored_regions(args.B, args.I, out_svg, out_pdf, individual_name, ancestries)


if __name__ == "__main__":
    main()
