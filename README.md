RFMIX2 → BOVINE: Professionalized plotting pipeline (RFMix2 → BED → SVG/PDF)
===============================================================================

Authors:
  Alessandro Lisi (alisi@usc.edu)
  Michael C. Campbell (mc44680@usc.edu)

Affiliation:
  Human Evolutionary Genomics Laboratory
  Department of Biological Sciences
  University of Southern California


OVERVIEW
-------------------------------------------------------------------------------

This repository provides a fully redesigned and production-ready pipeline to:

  1. Process RFMix v2 (.msp.tsv) output files.
  2. Combine per-chromosome outputs.
  3. Generate per-haplotype BED files.
  4. Create unified ancestry BED files with color and feature encoding.
  5. Insert colored local ancestry segments into a bovine genome SVG template.
  6. Export high-resolution, vector-quality PDF chromosomal paintings.

This version is specifically adapted for the bovine genome (29 autosomes)
(Bos taurus / Bos indicus).

Compared to the original human-based version, this release includes:

  • Robust file discovery
  • Streaming-safe I/O for large MSP files
  • Natural chromosome sorting (1..29)
  • Multi-sample MSP support
  • Parallel processing
  • JSON color configuration
  • Optional feature highlighting
  • Illustrator-ready SVG/PDF output
  • Professionalized logging and error handling


===============================================================================
PART I — RUNNING RFMIX2
===============================================================================

1. PREPARE INPUT DATASET
-------------------------------------------------------------------------------

Target and reference datasets must:

  • Be phased
  • Be in VCF or BCF format
  • Contain consistent SNP ordering
  • Use the same reference genome assembly


2. CREATE A REFERENCE SAMPLE MAP FILE
-------------------------------------------------------------------------------

The sample map file must:

  • Be tab-delimited
  • Contain ONLY reference individuals
  • Match VCF sample order exactly

Example:

<pre><code>
ET_Abergelle01    Indicine
ET_Abergelle02    Indicine
FR_Holstein01     Taurine
FR_Holstein02     Taurine
</code></pre>


3. GENETIC MAP FILE
-------------------------------------------------------------------------------

Must contain all chromosomes together (1–29 for bovine).

Columns (tab-delimited):

<pre><code>
chr    pos    cM
</code></pre>

Example:

<pre><code>
chr1    55550     0
chr1    82571     0.080572
chr1    285245    0.439456
</code></pre>


4. EXECUTE RFMIX2
-------------------------------------------------------------------------------

Example for chromosomes 1–29:

<pre><code>
for i in {1..29};
do
  rfmix \
    -f Target_Phased.vcf.gz \
    -r Reference_Phased.vcf.gz \
    -m Reference_Map.txt \
    -g genetic_map.txt \
    -o Output/sample_ \
    --chromosome=${i}
done
</code></pre>

RFMix2 generates:

  • *.rfmix.Q       → global ancestry
  • *.msp.tsv       → local ancestry (CRF output)
  • *.sis.tsv
  • *.tsv           → marginal probabilities

The bovine LAP pipeline uses *.msp.tsv files.


===============================================================================
PART II — BOVINE LOCAL ANCESTRY PIPELINE
===============================================================================

Main scripts:

  • LAP_Bovine_Genome.py
  • Plot_LAP_Bovine_Genome.py


-------------------------------------------------------------------------------
SCRIPT 1: LAP_Bovine_Genome.py
-------------------------------------------------------------------------------

Purpose:
  Convert RFMix2 MSP files into final ancestry BED files.

Pipeline flow:

  1. Discover MSP chunks per chromosome.
  2. Combine them into a single MSP per individual.
  3. Convert MSP → hap1 / hap2 BED files.
  4. Merge hap BEDs into final rendering BED.
  5. Optionally insert feature highlight lines.

Key Improvements:

  • Accepts per-chromosome MSPs OR single combined MSP.
  • Supports multi-sample MSP files.
  • Parallel processing via --threads.
  • JSON-based color configuration.
  • Optional feature highlighting.
  • Automatic chromosome natural sorting.
  • Safe streaming for large files.


Basic usage:

<pre><code>
python LAP_Bovine_Genome.py \
  --prefix results/sample_ \
  --chr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29 \
  --output-dir out_bed
</code></pre>

Parallel version:

<pre><code>
python LAP_Bovine_Genome.py \
  --prefix results/sample_ \
  --output-dir out_bed \
  --threads 4
</code></pre>

Feature highlight:

<pre><code>
python LAP_Bovine_Genome.py \
  --prefix results/sample_ \
  --output-dir out_bed \
  --from-bp 135000000 \
  --to-bp 135200000 \
  --chromosome chr2
</code></pre>


Output:
  • sampleA.bed
  • sampleB.bed
  • etc.


-------------------------------------------------------------------------------
SCRIPT 2: Plot_LAP_Bovine_Genome.py
-------------------------------------------------------------------------------

Purpose:
  Insert colored BED regions into a bovine chromosome SVG template
  and export a vector PDF.

Requirements:

  Python ≥ 3.8
  pandas
  rsvg-convert (librsvg)

Install librsvg:

macOS:
<pre><code>
brew install librsvg
</code></pre>

Linux:
<pre><code>
sudo apt-get install -y librsvg2-bin
</code></pre>


Basic usage:

<pre><code>
python Plot_LAP_Bovine_Genome.py \
  -B base_bovine_template.svg \
  -I out_bed/sampleA.bed \
  -O figures/sampleA.pdf
</code></pre>


What the plotting script does:

  • Reads BED file (geom_rect / geom_line entries)
  • Scales rectangles proportionally to bovine chromosome lengths (1–29)
  • Bottom-aligns chromosomes
  • Inserts SVG rectangle elements
  • Builds legend automatically
  • Converts SVG → PDF using rsvg-convert


-------------------------------------------------------------------------------
BOVINE GENOME SUPPORT
-------------------------------------------------------------------------------

The script includes bovine chromosome lengths (1–29).
To adapt to a different assembly, edit:

  chromosome_lengths = [ ... ]

inside Plot_LAP_Bovine_Genome.py.


-------------------------------------------------------------------------------
FULL WORKFLOW EXAMPLE
-------------------------------------------------------------------------------

Step 1: Generate final BED

<pre><code>
python LAP_Bovine_Genome.py \
  --prefix results/sample_ \
  --chr chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
        chr20 chr21 chr22 chr23 chr24 chr25 chr26 chr27 chr28 chr29 \
  --output-dir out_bed \
  --threads 4
</code></pre>

Step 2: Plot

<pre><code>
python Plot_LAP_Bovine_Genome.py \
  -B base_bovine_template.svg \
  -I out_bed/sampleA.bed \
  -O figures/sampleA.pdf
</code></pre>


-------------------------------------------------------------------------------
TROUBLESHOOTING
-------------------------------------------------------------------------------

"No individuals found":
  → Check prefix and chromosome tokens.

Missing ancestry colors:
  → Use --color-config JSON.

Plot misalignment:
  → Adjust x_base, x_step, top_y inside plotting script.

Large MSP files:
  → Ensure sufficient disk space.
  → Use --keep-temp for debugging.


-------------------------------------------------------------------------------
OUTPUT QUALITY
-------------------------------------------------------------------------------

  • Fully vector PDF
  • Illustrator-compatible
  • Suitable for publication
  • High resolution
  • Editable legend and rectangles


-------------------------------------------------------------------------------
CONTACT
-------------------------------------------------------------------------------

Alessandro Lisi
alisi@usc.edu

Michael C. Campbell
mc44680@usc.edu


-------------------------------------------------------------------------------
VERSION
-------------------------------------------------------------------------------

v1.0 — Bovine genome professional release

  • Multi-sample MSP support
  • Parallel processing
  • Feature highlighting
  • JSON color mapping
  • Robust SVG → PDF vector export
