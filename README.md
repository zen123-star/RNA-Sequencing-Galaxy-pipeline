 # Reference-based RNA-Seq Data Analysis Pipeline

A complete walkthrough of an RNA-Seq differential gene expression analysis pipeline run on *Drosophila melanogaster* data — from raw FASTQ reads to functional pathway enrichment.

---

## Table of Contents

- [Project Summary](#project-summary)
- [Biological Background](#biological-background)
- [Samples & Experimental Design](#samples--experimental-design)
- [Pipeline Overview](#pipeline-overview)
- [Step 1 — Quality Control](#step-1--quality-control)
- [Step 2 — Spliced Mapping to Reference Genome](#step-2--spliced-mapping-to-reference-genome)
- [Step 3 — Library Strandness Estimation](#step-3--library-strandness-estimation)
- [Step 4 — Read Counting per Gene](#step-4--read-counting-per-gene)
- [Step 5 — Differential Gene Expression Analysis](#step-5--differential-gene-expression-analysis)
- [Step 6 — Functional Enrichment Analysis](#step-6--functional-enrichment-analysis)
- [Key Results](#key-results)
- [File Formats Used](#file-formats-used)
- [Tools & References](#tools--references)

---

## Project Summary

This pipeline identifies **differentially expressed (DE) genes** between *Pasilla (PS) gene-depleted* and normal *Drosophila melanogaster* cells using RNA-Seq data. Starting from raw sequencing reads, the analysis performs quality filtering, genome mapping, read quantification, statistical DE testing, and finally functional enrichment to understand the biological impact of PS depletion.

**Organism:** *Drosophila melanogaster*
**Reference genome:** dm6
**Original study:** Brooks et al. 2011 — *Conservation of an RNA regulatory map between Drosophila and mammals*

---

## Biological Background

The **Pasilla (PS)** gene encodes a splicing regulator in *Drosophila*, homologous to the mammalian **Nova-1** and **Nova-2** proteins. In the original study, PS was depleted by RNA interference (RNAi) to identify which genes and molecular pathways it regulates.

- **Treated samples:** Pasilla gene knocked down by RNAi
- **Untreated samples:** Normal expression of Pasilla

Comparing RNA-Seq reads between these two conditions reveals which genes change expression when PS is absent.

---

## Samples & Experimental Design

| Sample ID | Condition | Library Type |
|---|---|---|
| GSM461176 | Untreated | Single-end |
| GSM461177 | Untreated | Paired-end |
| GSM461178 | Untreated | Paired-end |
| GSM461179 | Treated (PS depleted) | Single-end |
| GSM461180 | Treated (PS depleted) | Paired-end |
| GSM461181 | Treated (PS depleted) | Paired-end |
| GSM461182 | Untreated | Single-end |

**4 untreated replicates, 3 treated replicates.**
Two treated and two untreated samples are paired-end; the rest are single-end. This mix is accounted for as a second factor in the statistical model.

---

## Pipeline Overview

```
Raw FASTQ reads
      │
      ▼
 Quality Control          ← Falco, MultiQC, Cutadapt
      │
      ▼
 Spliced Mapping          ← STAR (dm6 genome + GTF annotation)
      │
      ▼
 Mapping QC               ← MultiQC (STAR logs), IGV / JBrowse2
      │
      ▼
 Strandness Estimation    ← Infer Experiment / STAR counts / pyGenomeTracks
      │
      ▼
 Read Counting            ← featureCounts
      │
      ▼
 Differential Expression  ← DESeq2 (multi-factor)
      │
      ▼
 Annotation & Filtering   ← DESeq2 annotator, Filter tool
      │
      ▼
 Visualization            ← heatmap2 (normalized counts + Z-scores)
      │
      ▼
 Functional Enrichment    ← goseq (GO + KEGG), Pathview
```

---

## Step 1 — Quality Control

### Purpose
Raw reads contain sequencing errors and adapter contamination. This step removes low-quality bases and short/poor-quality reads before mapping.

---

### 1.1 Per-sample Quality Check

**Tool:** Falco (FastQC-compatible quality report generator)

**Input:**
```
GSM461177_1.fastqsanger    ← forward reads
GSM461177_2.fastqsanger    ← reverse reads
GSM461180_1.fastqsanger
GSM461180_2.fastqsanger
```

**Output:**
- `RawData` — HTML quality report per sample (one per file)

**What to check in the report:**
- Per-base sequence quality (should be mostly green / Q > 28)
- Adapter content (flags if adapters are present)
- Sequence length distribution
- Read length (e.g., 36 bp for this dataset)

---

### 1.2 Aggregated QC Report

**Tool:** MultiQC

**Input:** Falco reports (RawData collection)

**Output:**
- `MultiQC Report` — Combined HTML report across all samples

**What to look for:**
- Quality drops toward the 3' end of reads (common in Illumina data)
- Overrepresented sequences or adapter contamination
- Consistent quality across all samples

---

### 1.3 Trimming and Filtering

**Tool:** Cutadapt

**Input:**
```
2 PE fastqs    ← paired-end collection (GSM461177 + GSM461180)
```

**Parameters:**
| Parameter | Value | Reason |
|---|---|---|
| Quality cutoff (R1) | 20 | Remove bases with Phred score < 20 |
| Minimum read length | 20 bp | Discard reads too short to map reliably |
| Input type | Paired-end collection | Processes both mates together |

**Output:**
- `Reads` — Trimmed paired FASTQ collection (used in next step)
- `Report` — Trimming statistics per sample

> **Note:** Paired-end samples are trimmed together in a single run. This ensures both mates of a pair are treated consistently — if one read is discarded, its partner is too.

---

## Step 2 — Spliced Mapping to Reference Genome

### Why spliced mapping?

In eukaryotes, mature mRNA has introns removed (spliced out). A standard aligner would fail to map reads that span two exons because the sequence doesn't exist as-is in the genomic DNA. **Splice-aware aligners** identify exon-exon junction reads and correctly map them across intron gaps.

```
Genome:    ████ exon1 ████----intron----████ exon2 ████
Read:                    ████████████████
                         (spans junction → needs spliced mapping)
```

---

### 2.1 Genome Mapping

**Tool:** STAR (Spliced Transcripts Alignment to a Reference)

**Input:**
```
Cutadapt output → Reads (trimmed paired-end collection)
Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz    ← gene annotation
```

**Reference genome:** dm6 (built-in)

**Key parameters:**
| Parameter | Value | Reason |
|---|---|---|
| Read type | Paired-end (as collection) | Matches our library prep |
| Reference genome | dm6 Full | Drosophila reference |
| GTF annotation | BDGP6.32.109 | Defines splice sites |
| Junction overhang length | 36 | Read length − 1 (for 37 bp reads) |
| Per gene read counts | GeneCounts | Generates counts alongside mapping |
| Coverage output | Bedgraph (strand 1 & 2) | Used for strandness estimation |

**Output:**
| File | Description |
|---|---|
| `mapped.bam` | Aligned reads in binary format |
| `log` | Mapping statistics (used for QC) |
| `reads per gene` | Raw gene-level counts (3 columns: unstranded, forward, reverse) |
| `Coverage Uniquely mapped strand 1` | Forward strand coverage bedgraph |
| `Coverage Uniquely mapped strand 2` | Reverse strand coverage bedgraph |

---

### 2.2 Mapping Quality Check

**Tool:** MultiQC (on STAR logs)

**Input:** STAR log collection

**Output:** `MultiQC Report` — mapping statistics summary

**Key metrics to evaluate:**

| Metric | Expected | Action if not met |
|---|---|---|
| % uniquely mapped reads | ~80% | < 70% → investigate contamination |
| % multi-mapped reads | < 10% | High multi-mapping → check genome version |
| % unmapped reads | < 20% | High unmapped → check trimming or wrong genome |

> For this dataset: ~80% of reads mapped uniquely for both samples — acceptable.

---

### 2.3 Visual Inspection of Alignments

**Tool:** IGV or JBrowse2

**Input:** `mapped.bam` collection

**Region to inspect:** chr4:540,000–560,000

**What to observe:**
- **Grey coverage peaks** — read depth across the region; higher = more expression
- **Connecting arcs/lines between reads** — these are spliced reads spanning introns; the arc represents the skipped intron
- **Sashimi plots** — arcs with numbers represent splice junctions; numbers = how many reads support each junction
- IGV:
<img width="1136" height="233" alt="igv_panel" src="https://github.com/user-attachments/assets/860d1b09-7a3a-4aa2-aff5-c7d1b945464a" />

<img width="1136" height="680" alt="igv_snapshot" src="https://github.com/user-attachments/assets/45947ec2-bfb7-4b85-8055-6d14dfb5c886" />
Sashmi plot:
<img width="930" height="733" alt="Sashimi" src="https://github.com/user-attachments/assets/9e58975d-6c1c-499f-99d5-ec77f1377575" />
Jbrowse:
<img width="1246" height="975" alt="jbrowse" src="https://github.com/user-attachments/assets/2453899e-adef-4d5f-bc4a-6f44a96cd4f5" />


---

## Step 3 — Library Strandness Estimation

### Why strandness matters

Some library preparation protocols preserve information about which strand the original RNA came from. This affects how reads are counted per gene — particularly for reads in regions where two genes on opposite strands overlap.

**Three possible outcomes:**

| Library type | Description |
|---|---|
| Unstranded | Reads map to genes on both strands equally |
| Stranded forward | Reads map predominantly to the same strand as the gene |
| Stranded reverse | Reads map predominantly to the opposite strand of the gene |

---

### Method: Infer Experiment (RSeQC)

**Tool:** Infer Experiment

**Input:**
```
mapped.bam             ← STAR output (BAM collection)
BED12 annotation file  ← converted from GTF using Convert GTF to BED12 tool
```

**Parameters:**
- Number of reads sampled: 200,000

**Output:** Text file per sample with fractions:
```
Fraction of reads explained by "1++,1--,2+-,2-+"  →  forward strand
Fraction of reads explained by "1+-,1-+,2++,2--"  →  reverse strand
```

**Interpretation:**
- If both fractions are close (~0.5 each) → **Unstranded**
- If one fraction dominates (> 0.8) → **Stranded** (forward or reverse)

> For GSM461177: both fractions ≈ 0.5 → **Unstranded library**

---

### Alternative Method: MultiQC on STAR GeneCounts

STAR outputs read counts for all three strandness scenarios simultaneously. By comparing the totals:
- The scenario with the **highest assigned read count** matches your library type

---

## Step 4 — Read Counting per Gene

### Purpose

Quantify how many reads (fragments) map to each annotated gene. This produces the raw count matrix used for differential expression testing.

---

### 4.1 Counting with featureCounts

**Tool:** featureCounts

**Input:**
```
mapped.bam (STAR output, BAM collection)
Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz    ← gene annotation
```

**Parameters:**
| Parameter | Value |
|---|---|
| Strandness | Unstranded (determined in Step 3) |
| Feature type | exon |
| Gene identifier | gene_id |
| Count paired reads as | 1 fragment |
| Minimum mapping quality | 10 |
| Output format | Gene-ID + read-count (tab-delimited) |
| Create gene-length file | Yes |

**Output:**
| File | Description |
|---|---|
| `Counts` | Tab-delimited: Gene ID → read count per sample |
| `Feature lengths` | Gene ID → length in bp (required for goseq later) |
| `Summary` | Assignment statistics (% assigned, unassigned) |

**Format of counts file:**
```
GeneID          GSM461177
FBgn0000003     0
FBgn0000008     44
FBgn0000014     0
FBgn0000015     3
...
```

> **Tip:** An assignment rate below 50% is a warning sign. Causes include wrong genome/annotation version mismatch, incorrect strandness setting, or poor mapping quality.

---

### 4.2 Why We Can't Compare Raw Counts Directly

Before differential expression testing, we need to understand why simple count comparison is misleading:

| Bias | Cause | Effect |
|---|---|---|
| Sequencing depth | Different samples have different total reads | Deeper samples have higher counts for all genes |
| Gene length | Longer genes produce more reads | Long genes appear more expressed |
| Library composition | One condition may express unique high-abundance genes | Remaining genes appear falsely downregulated |

DESeq2 (next step) corrects for all of these during normalization.

---

## Step 5 — Differential Gene Expression Analysis

### 5.1 Importing All Count Files

For the DESeq2 step, pre-computed count files for all 7 samples are used:

```
GSM461176_untreat_single_featureCounts.counts
GSM461177_untreat_paired_featureCounts.counts
GSM461178_untreat_paired_featureCounts.counts
GSM461179_treat_single_featureCounts.counts
GSM461180_treat_paired_featureCounts.counts
GSM461181_treat_paired_featureCounts.counts
GSM461182_untreat_single_featureCounts.counts
```

---

### 5.2 Running DESeq2

**Tool:** DESeq2

**Input:** 7 count files (above)

**Experimental design — two factors:**

| Factor | Level 1 | Level 2 |
|---|---|---|
| Treatment (primary) | treated → GSM461179, 180, 181 | untreated → GSM461176, 177, 178, 182 |
| Sequencing (covariate) | PE → GSM461177, 178, 180, 181 | SE → GSM461176, 179, 182 |

> Including the sequencing type as a second factor removes technical variation caused by library preparation method, making the treatment effect estimate cleaner.

**Additional setting:** Beta priors enabled (`Yes`)

**Output files:**

| Output | Description |
|---|---|
| `DESeq2 result` | Per-gene statistics table (see columns below) |
| `Normalized counts` | Normalized expression values for all genes × all samples |
| `Plots` | PCA, distance heatmap, dispersion plot, p-value histogram, MA plot |

**DESeq2 result table columns:**

| Column | Description |
|---|---|
| GeneID | Flybase gene ID (e.g., FBgn0003360) |
| Base mean | Mean normalized count across all 7 samples |
| log₂(FC) | Log₂ fold change: treated vs untreated (positive = up in treated) |
| StdErr | Standard error of the log₂FC estimate |
| Wald-Stat | Test statistic |
| P-value | Raw p-value |
| P-adj | Benjamini-Hochberg adjusted p-value (FDR-corrected) |

---

### 5.3 Interpreting DESeq2 Diagnostic Plots

**PCA plot:**
- PC1 (largest source of variation) separates **treated vs untreated** — confirms treatment has a strong effect
- PC2 separates **paired-end vs single-end** — confirms the second factor was necessary

**Sample distance heatmap:**
- Samples cluster by treatment condition
- Treated samples are more similar to each other than to untreated, and vice versa

**MA plot:**
- X-axis: average expression (log scale)
- Y-axis: log₂ fold change
- Blue dots: genes with adjusted p-value < 0.1
- Shows that significant DE genes are spread across all expression levels

---

### 5.4 Annotation of Results

**Tool:** Annotate DESeq2/DEXSeq output tables

**Input:**
```
DESeq2 result file
Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
```

**Output:** Annotated results table — adds columns:

| Added Column | Description |
|---|---|
| Chromosome | Chromosomal location |
| Start / End | Genomic coordinates |
| Strand | + or − strand |
| Feature | Feature type |
| Gene symbol | Human-readable name (e.g., `ps`, `Sxl`) |

---

### 5.5 Filtering DE Genes

Two sequential filters applied to the annotated results:

**Filter 1 — Statistical significance:**
```
Condition: adjusted p-value < 0.05
Column:    c7 < 0.05
```

**Filter 2 — Biological relevance (fold change threshold):**
```
Condition: |log₂FC| > 1  (equivalent to FC > 2 or FC < 0.5)
Column:    abs(c3) > 1
```

**Output:** `Genes with significant adj p-value & abs(log2(FC)) > 1`
→ ~113 genes pass both filters

---

### 5.6 Visualization — Heatmap of Normalized Counts

**Tool:** heatmap2

**Input:**
```
Normalized counts for the most DE genes    ← joined from DESeq2 normalized counts
                                              + filtered DE gene list (113 genes × 7 samples)
```

**Settings:**
| Parameter | Value |
|---|---|
| Data transformation | Log2(value + 1) |
| Clustering | Enabled |
| Axis labels | Label columns (samples) only |
| Color scheme | 2-color gradient |

**Output:** Heatmap PNG
<img width="840" height="843" alt="image" src="https://github.com/user-attachments/assets/d7e409fb-7a51-4c13-a218-0db7b1f1b271" />


- X-axis: 7 samples
- Y-axis: 113 DE genes (clustered)
- Treated samples cluster together; untreated cluster together
- Expression patterns split clearly along the treatment axis

---

### 5.7 Visualization — Z-score Heatmap

**Tool:** heatmap2

**Input:** Same normalized counts table

**Settings:**
| Parameter | Value |
|---|---|
| Data transformation | Plot as-is |
| Z-score computation | On rows (per gene) |
| Clustering | Enabled |
| Color scheme | 3-color gradient |

**What Z-score shows:**
For each gene, the Z-score is how many standard deviations a sample's expression is from that gene's average across all samples. This removes the effect of absolute expression level and highlights relative patterns.

```
Z-score formula:  z = (x_ij − mean_i) / sd_i
```
<img width="832" height="838" alt="image" src="https://github.com/user-attachments/assets/59364dfb-7e7c-43ff-bd2f-50f8a14de394" />

---

## Step 6 — Functional Enrichment Analysis

### 6.1 Gene Ontology (GO) Analysis

**Purpose:** Identify which biological functions, cellular components, or molecular processes are significantly enriched among the DE genes.

**Tool:** goseq (corrects for gene-length bias inherent to RNA-Seq data)

---

**Prepare input 1 — DE gene boolean list:**

Starting from the DESeq2 result file:

1. **Compute** — add column: `bool(float(c7) < 0.05)` → True/False per gene
2. **Cut** — keep columns c1 (Gene ID) and c8 (True/False)
3. **Change Case** — convert Gene IDs to uppercase

**Output:** `Gene IDs and differential expression`
```
FBGN0000008    False
FBGN0000014    False
FBGN0003360    True
...
```

---

**Prepare input 2 — Gene length file:**

From featureCounts feature lengths output:

1. **Extract** the first dataset from the length collection
2. **Change Case** — convert Gene IDs to uppercase

**Output:** `Gene IDs and length`
```
FBGN0000008    1059
FBGN0000014    3310
FBGN0003360    2031
...
```

---

**Run goseq:**

**Input:**
```
Gene IDs and differential expression
Gene IDs and length
```

**Parameters:**
| Parameter | Value |
|---|---|
| Genome | Fruit fly (dm6) |
| Gene ID format | Ensembl Gene ID |
| Categories | GO: Cellular Component, Biological Process, Molecular Function |
| Output top GO terms plot | Yes |
| Extract DE genes per category | Yes |

**Output files:**

| File | Description |
|---|---|
| `Ranked category list` | All GO terms with over/under-representation statistics |
| `Top GO terms plot` | Bar chart: top 10 over-represented GO terms |
| `DE genes for categories` | Which DE genes belong to which GO terms |

**Output table columns:**

| Column | Description |
|---|---|
| category | GO term ID |
| over_rep_pval | p-value for over-representation |
| under_rep_pval | p-value for under-representation |
| numDEInCat | DE genes in this category |
| numInCat | All genes in this category |
| term | GO term description |
| ontology | BP / MF / CC |
| p.adjust.over_represented | BH-adjusted p-value (use this for significance) |
| p.adjust.under_represented | BH-adjusted p-value |

---

### 6.2 KEGG Pathway Analysis

**Tool:** goseq (same inputs as GO analysis)

**Parameters:** Same as above, but select `KEGG` as the category instead of GO.

**Output:**
| File | Description |
|---|---|
| `Ranked category list` | All KEGG pathways with statistics |
| `DE genes for categories` | DE genes mapped to each pathway |

**Notable pathways in this dataset:**

| KEGG ID | Pathway | Result |
|---|---|---|
| dme00010 | Glycolysis / Gluconeogenesis | Over-represented (adj. p < 0.05) |
| dme03040 | Spliceosome | Under-represented (biologically relevant: PS is a splicing regulator) |

---

### 6.3 KEGG Pathway Visualization

**Tool:** Pathview

**Purpose:** Overlay log₂ fold change values onto KEGG pathway diagrams to visualize which genes in a pathway are up/downregulated.

**Input 1 — Pathway IDs to plot:**
```
00010
03040
```

**Input 2 — Gene expression data:**
- Columns: Gene ID + log₂FC
- Source: `Genes with significant adj p-value` → Cut columns c1, c3

**Parameters:**
| Parameter | Value |
|---|---|
| Species | Fly |
| Gene data format | Ensembl Gene ID |
| Output format | KEGG native |
| Plot on same layer | Yes |

**Output:** One pathway image per ID
<img width="663" height="912" alt="image" src="https://github.com/user-attachments/assets/7cd94a78-4b32-450b-b62d-48f4fde35df4" />
<img width="840" height="760" alt="image" src="https://github.com/user-attachments/assets/36bc2ac9-1a1c-499b-a5c9-0a99f25a0bc7" />

**Interpretation of pathway images:**
- Each **box** = a gene/enzyme in the pathway
- **Green box** = downregulated in treated samples
- **Red box** = upregulated in treated samples
- **Grey box** = gene present but not significantly DE
- **White box** = gene not detected in our dataset

---

## Key Results

| Analysis Stage | Result |
|---|---|
| Read mapping rate | ~80% of reads mapped uniquely (both test samples) |
| Library strandness | Unstranded |
| Total genes tested | ~14,000 |
| Genes with adj p < 0.05 | Several hundred |
| Genes with adj p < 0.05 AND \|log₂FC\| > 1 | ~113 |
| Pasilla gene (FBgn0261552) | Downregulated in treated samples ✓ (validates experiment) |
| Top over-represented GO terms | RNA binding, splicing, nucleic acid processing |
| Significant KEGG pathway | dme00010 — Glycolysis / Gluconeogenesis |

---

## File Formats Used

| Format | Extension | Description |
|---|---|---|
| FASTQ | `.fastqsanger` | Raw sequencing reads + quality scores |
| BAM | `.bam` | Binary file of read-to-genome alignments |
| GTF | `.gtf.gz` | Gene annotation (exon positions, gene IDs) |
| BED12 | `.bed` | 12-column annotation format (for Infer Experiment) |
| Counts | `.counts` | Tab-delimited: Gene ID → integer read count |
| Bedgraph | `.bedgraph` | Per-base coverage across the genome |

---

## Tools & References

| Tool | Purpose | Citation |
|---|---|---|
| Falco | Read quality reports | — |
| MultiQC | Aggregate QC reports | Ewels et al. 2016 |
| Cutadapt | Adapter trimming | Marcel 2011 |
| STAR | Spliced genome mapping | Dobin et al. 2013 |
| featureCounts | Read counting | Liao et al. 2013 |
| Infer Experiment (RSeQC) | Strandness estimation | Wang et al. 2012 |
| DESeq2 | Differential expression | Love et al. 2014 |
| heatmap2 | Heatmap visualization | — |
| goseq | GO & KEGG enrichment | Young et al. 2010 |
| Pathview | KEGG pathway overlay | Luo & Brouwer 2013 |

**Dataset:** Brooks AN et al. (2011). *Conservation of an RNA regulatory map between Drosophila and mammals.* Genome Research 21:193–202. GEO accession: GSE18508.

---

*Pipeline completed as part of the Galaxy Training Network RNA-Seq tutorial (CC-BY 4.0). Tutorial PURL: https://gxy.io/GTN:T00295*

