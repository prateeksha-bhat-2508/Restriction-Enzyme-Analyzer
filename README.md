## How to Run

### Requirements

- Python 3.10+
- pip

### Installation

```bash
# Clone the repository
git clone https://github.com/ShreeramHegde01/restriction-enzyme-analyzer.git
cd restriction-enzyme-analyzer

# Create virtual environment
python -m venv bioenv

# Activate (Windows)
bioenv\Scripts\activate

# Activate (macOS/Linux)
source bioenv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run app.py
```

The app opens at `http://localhost:8501`

### Dependencies

| Package | Purpose |
|---------|---------|
| streamlit | Web interface |
| biopython | Restriction analysis, NCBI access |
| plotly | Interactive charts |
| pandas | Data tables |

---
# REA — Restriction Enzyme Analyzer

A Biopython-powered web tool for restriction enzyme analysis, cut site detection, and virtual gel electrophoresis simulation.

---

## Table of Contents

1. [How to Run](#how-to-run)
2. [Biological Significance](#biological-significance)
3. [Real-World Applications](#real-world-applications)
4. [Biological Examples](#biological-examples)
5. [How Biopython Is Used](#how-biopython-is-used)
6. [Features](#features)
7. [Workflow](#workflow)
8. [Challenges and Solutions](#challenges-and-solutions)
9. [Limitations and Future Work](#limitations-and-future-work)
10. [References](#references)
---

## Biological Significance

### What Are Restriction Enzymes?

Restriction enzymes (restriction endonucleases) are bacterial proteins that recognize specific DNA sequences and cut the double helix at precise locations. They are part of the bacterial immune system, protecting against viral infections by destroying foreign DNA.

### Why Are They Important?

The discovery of restriction enzymes in the 1970s revolutionized molecular biology, earning Werner Arber, Daniel Nathans, and Hamilton O. Smith the **1978 Nobel Prize in Physiology or Medicine**. These molecular scissors became the foundation of:

- **Recombinant DNA technology** — Creating genetically modified organisms
- **Gene cloning** — Inserting genes into vectors for expression
- **DNA fingerprinting** — Forensic identification and paternity testing
- **Genetic disease diagnosis** — Detecting mutations that alter cut sites

### Recognition Sites and Overhangs

Each restriction enzyme recognizes a specific **palindromic sequence** (reads the same on both strands 5'→3'):

| Enzyme | Recognition Site | Overhang Type |
|--------|-----------------|---------------|
| EcoRI | G↓AATTC | 5' overhang (4 nt) |
| BamHI | G↓GATCC | 5' overhang (4 nt) |
| HindIII | A↓AGCTT | 5' overhang (4 nt) |
| EcoRV | GAT↓ATC | Blunt end |
| PstI | CTGCA↓G | 3' overhang (4 nt) |

The type of overhang (sticky ends vs blunt ends) determines which DNA fragments can be joined together during cloning.

---

## Real-World Applications

### 1. Gene Cloning and Expression

Scientists use restriction enzymes to cut both a gene of interest and a plasmid vector at compatible sites, then ligate them together. This is how:
- **Insulin** is produced in bacteria for diabetic patients
- **Growth hormone** is manufactured for medical use
- **Vaccines** are developed using recombinant proteins

### 2. Disease Diagnosis (RFLP Analysis)

Single nucleotide mutations can create or destroy restriction sites. This allows genetic testing for:
- **Sickle cell anemia** — The HbS mutation destroys a DdeI site in the β-globin gene
- **Huntington's disease** — RFLP markers linked to the HD gene
- **Cystic fibrosis** — Detection of ΔF508 mutation

### 3. Forensic DNA Analysis

Before modern STR profiling, RFLP analysis of Variable Number Tandem Repeats (VNTRs) was the standard for DNA fingerprinting in criminal investigations and paternity testing.

### 4. Plasmid Verification

After cloning, scientists digest plasmids with restriction enzymes and run the fragments on a gel to verify:
- Correct insert size
- Proper orientation
- Plasmid integrity

---

## Biological Examples

### Example 1: BRCA1 Gene Analysis

The **BRCA1 gene** (Breast Cancer susceptibility gene 1) is a tumor suppressor located on chromosome 17. Mutations in BRCA1 significantly increase breast and ovarian cancer risk.

Using REA, researchers can:
- Fetch the BRCA1 mRNA sequence (NM_007294) from NCBI
- Identify restriction sites for cloning BRCA1 into expression vectors
- Plan diagnostic restriction digests to detect known mutations

### Example 2: pUC19 Plasmid

**pUC19** is one of the most widely used cloning vectors in molecular biology. Its Multiple Cloning Site (MCS) contains recognition sequences for many common restriction enzymes clustered in a small region, enabling flexible gene insertion.

REA includes a pUC19 MCS sample to demonstrate how multiple enzymes cut within this important cloning region.

### Example 3: Lambda Phage DNA

**Bacteriophage λ (lambda)** DNA is commonly used as a size standard in gel electrophoresis. When cut with HindIII, it produces fragments of known sizes (23.1, 9.4, 6.6, 4.4, 2.3, 2.0, 0.56 kb) used to calibrate gels.

---

## How Biopython Is Used

REA leverages three key Biopython modules:

### 1. Bio.Restriction — Core Analysis Engine

The `Bio.Restriction` module contains a database of **600+ restriction enzymes** with their recognition sites, cut positions, and overhang information.

- `RestrictionBatch` — Searches multiple enzymes simultaneously
- `enzyme.site` — Returns the recognition sequence
- `enzyme.ovhg` — Returns overhang length (+5', -3', 0=blunt)

### 2. Bio.Entrez — NCBI Database Access

The `Bio.Entrez` module enables live fetching of sequences from NCBI GenBank:

- Fetch any sequence by accession number (e.g., NM_007294 for BRCA1)
- Retrieve full GenBank records with annotations
- Access the world's largest nucleotide database

### 3. Bio.SeqIO — Sequence File Parsing

The `Bio.SeqIO` module handles multiple file formats:

- **FASTA** — Simple sequence format (.fasta, .fa, .fna)
- **GenBank** — Annotated sequences (.gb, .gbk)
- Automatic format detection and parsing

---

## Features

| Feature | Description |
|---------|-------------|
| Multiple Input Methods | Paste sequence, upload files (FASTA/GenBank), samples, or NCBI fetch |
| 18 Restriction Enzymes | EcoRI, HindIII, BamHI, NotI, XhoI, and 13 more common enzymes |
| Cut Site Map | Interactive visualization of cut positions on DNA backbone |
| Virtual Gel | Simulated agarose gel electrophoresis with DNA ladder |
| Fragment Analysis | Size distribution of restriction fragments |
| GC Content Profile | Sliding window analysis of GC percentage |
| Enzyme Details | Recognition sites, overhang types, cut counts |

---

## Workflow

```
┌────────────────────────────────────────────────────────────┐
│                      INPUT                                  │
│  Paste | Upload File | Sample | NCBI Fetch                 │
└────────────────────────┬───────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│                 VALIDATION                                  │
│  • Check for valid DNA characters (ATGC)                   │
│  • Minimum sequence length (20 bp)                         │
│  • Clean whitespace and formatting                         │
└────────────────────────┬───────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│              BIOPYTHON ANALYSIS                             │
│  • Create Bio.Seq object                                   │
│  • Build RestrictionBatch with selected enzymes            │
│  • Search for all cut sites (linear DNA mode)              │
│  • Extract: positions, fragments, overhangs                │
└────────────────────────┬───────────────────────────────────┘
                         │
                         ▼
┌────────────────────────────────────────────────────────────┐
│               VISUALIZATION                                 │
│  • Cut Site Map (Plotly scatter)                           │
│  • Virtual Gel (log-scale migration)                       │
│  • Fragment Distribution (bar chart)                       │
│  • GC Profile (line chart)                                 │
│  • Enzyme Summary Table                                    │
└────────────────────────────────────────────────────────────┘
```

---


## Challenges and Solutions

### Challenge 1: Enzyme Object vs String Keys

**Problem:** Biopython returns enzyme objects as dictionary keys, not strings, making display and comparison difficult.

**Solution:** Created a mapping dictionary (`COMMON_ENZYMES`) linking string names to enzyme classes, and used `str(enzyme)` for consistent conversion.

### Challenge 2: Gel Electrophoresis Visualization

**Problem:** DNA fragment sizes span 50–5000+ bp. Linear scaling clusters small fragments together with no visual separation.

**Solution:** Applied logarithmic transformation (`math.log10(bp)`) to simulate real gel migration physics, where smaller fragments migrate farther.

### Challenge 3: Multi-Record FASTA Files

**Problem:** NCBI downloads often contain multiple sequences. `SeqIO.read()` fails with "more than one record" error.

**Solution:** Switched to `SeqIO.parse()` which returns an iterator, then take the first record and notify users if multiple sequences exist.

### Challenge 4: NCBI Connection Handling

**Problem:** NCBI API requires email, has rate limits, and can timeout unpredictably.

**Solution:** Wrapped fetch calls in try/except blocks, cached successful fetches in `st.session_state`, and displayed clear error messages.

### Challenge 5: Empty Results Handling

**Problem:** When no enzyme cuts the sequence, visualization functions crashed with empty data.

**Solution:** Added checks before all visualizations to display informative messages instead of blank/broken charts.

### Challenge 6: Plotly Color Format

**Problem:** Plotly doesn't accept 8-character hex colors for transparency (e.g., `#22C55E22`).

**Solution:** Converted to proper `rgba()` format: `rgba(34,197,94,0.15)`.

---

## Limitations and Future Work

### Current Limitations

- **Linear DNA only** — No circular plasmid support
- **Limited enzyme panel** — 18 enzymes from 600+ available
- **No methylation sensitivity** — Assumes unmethylated DNA
- **Single sequence** — No batch processing

### Future Improvements

- Circular plasmid mode with appropriate visualization
- Expanded enzyme search and filtering
- Double digest optimization for cloning
- Methylation-sensitive site prediction
- Batch sequence analysis

---

## References

1. Roberts, R.J. (2005). How restriction enzymes became the workhorses of molecular biology. *PNAS*, 102(17), 5905-5908.

2. Pingoud, A., & Jeltsch, A. (2001). Structure and function of type II restriction endonucleases. *Nucleic Acids Research*, 29(18), 3705-3727.

3. [Biopython Documentation](https://biopython.org/wiki/Documentation)

4. [REBASE — Restriction Enzyme Database](http://rebase.neb.com/)

---

*Biopython Project — Academic Demonstration*