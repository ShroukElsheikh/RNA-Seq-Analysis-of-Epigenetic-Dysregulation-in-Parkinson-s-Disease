# 🧬 Parkinson's Disease RNA-Seq Analysis Project
### Mitochondrial Stress-Induced H4K12 Hyperacetylation Dysregulates Transcription in Parkinson's Disease

---

## 📌 Project Overview

This project presents a complete bioinformatics analysis pipeline applied to RNA-Seq data from the **MitoPark mouse model of Parkinson's Disease (PD)**, conducted as part of two complementary courses:

- **Omics Course** — Downstream analysis (Differential Expression, Visualization)
- **NGS Course** — Upstream analysis (Quality Control, Trimming, Alignment, Quantification)

The study is based on the following published paper:

> Huang M, Jin H, Anantharam V, Kanthasamy A, Kanthasamy AG.
> **"Mitochondrial stress-induced H4K12 hyperacetylation dysregulates transcription in Parkinson's disease."**
> *Frontiers in Cellular Neuroscience*, 2024. DOI: [10.3389/fncel.2024.1422362](https://doi.org/10.3389/fncel.2024.1422362)

---

## 🎯 Biological Background

Parkinson's Disease (PD) is a progressive neurodegenerative disorder characterized by the irreversible loss of dopaminergic neurons in the substantia nigra. This project investigates the epigenetic mechanism underlying PD pathogenesis, specifically the role of **H4K12 acetylation (H4K12ac)** — a histone modification driven by mitochondrial stress — in dysregulating gene transcription.

The MitoPark mouse model used in this study carries a conditional knockout of the mitochondrial transcription factor **TFAM** specifically in dopaminergic neurons, mimicking the progressive mitochondrial impairment observed in human PD.

---

## 🐭 Experimental Design

| Sample | Condition | SRR Accession |
|--------|-----------|---------------|
| Control_1 | Wild-type C57BL/6 | SRR15194154 |
| Control_2 | Wild-type C57BL/6 | SRR15194155 |
| Control_3 | Wild-type C57BL/6 | SRR15194156 |
| PD_1 | MitoPark Transgenic | SRR15194157 |
| PD_2 | MitoPark Transgenic | SRR15194158 |
| PD_3 | MitoPark Transgenic | SRR15194159 |

- **Organism:** *Mus musculus* (Mouse)
- **Sequencing type we used:** Single-end RNA-Seq
- **Platform:** Illumina HiSeq X10
- **GEO Accession:** GSE180408

---

## 🔬 Analysis Pipeline

The full analysis was divided into two phases across two courses:

---

### 📘 Phase 1 — Upstream Analysis (NGS Course)

> 🔗 **NGS Course Presentation:** [View on Canva](https://canva.link/yjthpql15c2268h)

The upstream analysis covers all steps from raw sequencing data to gene-level quantification.

#### Steps performed:

**1. Quality Control (FastQC)**
- Ran FastQC on all 6 raw FASTQ files to assess read quality
- Evaluated per-base sequence quality, adapter content, GC content, and sequence duplication levels
- Generated HTML reports for each sample before trimming

**2. Read Trimming (Trimmomatic)**
- Trimmed adapter sequences and low-quality bases using Trimmomatic SE mode
- Parameters used:
  - `ILLUMINACLIP`: TruSeq3-SE adapters
  - `LEADING:3` `TRAILING:3`
  - `SLIDINGWINDOW:4:15`
  - `MINLEN:36`
- Generated trimmed FASTQ files for all 6 samples

**3. Post-trimming Quality Control (FastQC)**
- Re-ran FastQC on all trimmed files to confirm quality improvement
- Compared adapter content and per-base quality before and after trimming

**4. Transcript Quantification (Kallisto)**
- Built a Kallisto index from the mouse reference transcriptome
- Quantified transcript-level expression (TPM) for all 6 samples
- Generated `abundance.tsv` files for downstream analysis

#### Tools Used:
```
FastQC v0.11+
Trimmomatic v0.40
Kallisto v0.46+
```

---

### 📗 Phase 2 — Downstream Analysis (Omics Course)

> 🔗 **Omics Course Presentation:** [View on Canva](https://canva.link/iiqytos7m21agyi)

The downstream analysis covers differential expression analysis and biological interpretation of results.

#### Steps performed:

**1. TPM Matrix Construction**
- Merged TPM values from all 6 Kallisto output files into a single expression matrix
- Rows represent transcripts, columns represent samples

**2. Differential Expression Analysis (limma)**
- Applied log2 transformation to TPM values
- Used limma's linear model with empirical Bayes moderation
- Applied strict filtering thresholds:
  - `|LogFC| ≥ 1.5`
  - `adj.P.Val < 0.05`
- Identified **6 strictly significant DEGs**

**3. Visualization**
- Generated an **Enhanced Volcano Plot** showing the distribution of all genes by fold change and significance
- Generated an **Annotated Heatmap** of the top significant DEGs across all samples
- Generated an **Density PLot** for data before and after processing

**4. Biological Interpretation**
- Mapped significant DEGs to their gene functions
- Connected findings to the paper's central axis:
  **Mitochondrial impairment → H4K12ac deposition → Transcriptional dysregulation → Neurodegeneration**

#### Tools Used:
```
R (limma, EnhancedVolcano, writexl)
Python / Bash scripting
```

---

## 📊 Key Results

### Significant DEGs Identified (adj.P.Val < 0.05 & |LogFC| ≥ 1.5)

| Transcript ID | Direction | Associated Gene | Biological Role in PD |
|---------------|-----------|-----------------|----------------------|
| ENSMUST00000113331.8 | ⬆️ Upregulated | *Mmp12* | Neuroinflammation, ECM degradation |
| ENSMUST00000155633.8 | ⬆️ Upregulated | *Dffb* | DNA fragmentation, apoptosis |
| ENSMUST00000073691.5 | ⬆️ Upregulated | *Lima1* | Dysregulated in PD |
| ENSMUST00000133949.2 | ⬆️ Upregulated | *Atf3* | Antioxidant stress response failure |
| ENSMUST00000146616.2 | ⬆️ Upregulated | *Hvcn1* | ROS generation, neuroinflammation |
| ENSMUST00000187347.7 | ⬆️ Upregulated | *Ly6c1* | Dysregulated in PD |

### Volcano Plot
The volcano plot confirms a small but highly specific set of DEGs meeting strict significance thresholds, with ENSMUST00000113331.8 showing the highest fold change in PD samples. The majority of genes show biologically meaningful changes that do not reach statistical significance, consistent with the limited sample size.

### Heatmap
The heatmap demonstrates a perfect separation between Control and PD samples based on DEG expression profiles. All significant DEGs show consistently low expression (blue) in Control samples and high expression (red) in PD samples, supporting the paper's finding that H4K12ac-driven transcriptional activation promotes neuroinflammatory and apoptotic pathways in PD dopaminergic neurons.

---

## 📁 Repository Structure

```
Shrouk_NGS_Assignment/
│
├── fastq_files/                  # Raw downloaded FASTQ files (6 samples)
├── fastqc_before/                # FastQC reports before trimming
├── trimmed_files/                # Trimmed FASTQ files (Trimmomatic output)
├── trimmomatic_logs/             # Trimming statistics per sample
├── fastqc_after/                 # FastQC reports after trimming
├── kallisto_results/             # Kallisto quantification results
│   ├── SRR15194154_results/
│   ├── SRR15194155_results/
│   ├── SRR15194156_results/
│   ├── SRR15194157_results/
│   ├── SRR15194158_results/
│   └── SRR15194159_results/
├── merged_tpm_matrix.tsv         # Combined TPM expression matrix
├── metadata.txt                  # Sample condition annotations
├── Filtered_Significant_DEGs.xlsx # Excel file of significant DEGs
├── Volcano_Plot_Enhanced.png     # Volcano plot output
└── Heatmap_Annotated.png         # Heatmap output
```

---

## 🔗 Presentations

| Course | Presentation Link | Content |
|--------|-------------------|---------|
| 🔬 NGS Course | [View Presentation](https://canva.link/yjthpql15c2268h) | Upstream analysis: QC, trimming, alignment, quantification |
| 🧬 Omics Course | [View Presentation](https://canva.link/iiqytos7m21agyi) | Downstream analysis: DEG analysis, volcano plot, heatmap, biological interpretation |

---

## 📚 Reference

Huang M, Jin H, Anantharam V, Kanthasamy A, Kanthasamy AG.
*Mitochondrial stress-induced H4K12 hyperacetylation dysregulates transcription in Parkinson's disease.*
Front. Cell. Neurosci. 2024;18:1422362.
DOI: [10.3389/fncel.2024.1422362](https://doi.org/10.3389/fncel.2024.1422362)
GEO Dataset: [GSE180408](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180408)

---

## 👩‍🔬 Author

- **Judi Yousri 231002400**
- **Merolla Raafat  231001407** 
- **Shrouk ElSkeikh  231002025**
- **Mohamed Magdy 231000660**
- School of Biotechnology — Nile University
- Bioinformatics Program

