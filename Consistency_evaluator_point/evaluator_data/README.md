# Evaluator Data Directory

This directory contains the input data files used by the Consistency Request Evaluators (Sequences used for all 3 Evaluators are the same).

## Contents

### `all_consistency_data.csv`

The primary input file containing genomic sequences for consistency evaluation. This file includes:
- **150 original sequences** from ATAC-seq peaks (with reverse complements)
- **150 mononucleotide-shuffled sequences** (with reverse complements)
- **150 dinucleotide-shuffled sequences** (with reverse complements)

**Total: 900 sequences** (450 forward + 450 reverse complement)

### Additional Files

- `sequences.fasta` - Original ATAC-seq peak sequences (100kb centered on peak midpoints)
- `mononucleotide_shuffled.fa` - Mononucleotide-shuffled control sequences (generated via BiasAway)
- `dinucleotide_shuffled.fa` - Dinucleotide-shuffled control sequences (generated via BiasAway)
- `ENCFF121CAA.bed` - **ENCODE ATAC-seq data** for iPSC cells
- `sequence_design.py` - Code to generate `all_consistency_data.csv`

## File Format

The input file is a **tab-separated values (TSV)** file with the following structure:

```tsv
id	seq_name	sequence
1	chr1:12345-112345	ATCGATCGATCGATCG...
2	chr1:12345-112345_RC	CGATCGATCGATCGAT...
3	chr1:12345-112345_mononuc_shuffled	GCTAGCTAGCTAGCTA...
4	chr1:12345-112345_mononuc_shuffled_RC	TAGCTAGCTAGCTAGC...
```

### Column Specifications

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | `id` | Integer | Unique identifier for each sequence (used as index) |
| 2 | `seq_name` | String | Sequence identifier with genomic coordinates and type |
| 3 | `sequence` | String | DNA sequence (100,000 bp, A/T/C/G nucleotides) |

### Naming Conventions

| Suffix/Tag | Description | Count |
|------------|-------------|-------|
| `chrX:start-end` | Original ATAC-seq peak sequence | 150 |
| `_RC` | Reverse complement of any sequence | 450 |
| `_mononuc_shuffled` | Mononucleotide-shuffled control | 150 |
| `_dinuc_shuffled` | Dinucleotide-shuffled control | 150 |

### Important Requirements

1. **Paired Sequences**: Each forward sequence must have a corresponding reverse complement
2. **RC Identification**: Reverse complement sequences must include 'RC' in their `seq_name`
3. **Tab Separation**: Columns must be separated by tabs, not spaces or commas
4. **Header Row**: First row must contain column headers
5. **Sequence Length**: All sequences are 100,000 bp
6. **Valid Nucleotides**: Sequences contain only standard DNA bases (A, T, C, G, N)

## Data Generation Pipeline (optional)

The consistency evaluator dataset was created using a multi-step pipeline (`sequence_design.py`):

### Step 1: Extract ATAC-seq Peak Sequences

Original sequences were extracted from **ENCODE ATAC-seq data** for iPSC cells:
- **Source**: `ENCFF121CAA.bed` (ENCODE accession)
- **Cell Type**: Induced pluripotent stem cells (iPSC)
- **Species**: Homo sapiens (hg38 reference genome)

**Process:**
1. Load ATAC-seq peak coordinates from BED file
2. Remove duplicate peaks
3. Randomly sample 150 unique peaks
4. Calculate peak midpoint
5. Extract 100,000 bp centered on each peak midpoint (±50,000 bp)
6. Generate reverse complement for each sequence

### Step 2: Generate Shuffled Controls using Biasaway webserver

**Mononucleotide Shuffling:**
- Preserves nucleotide composition (GC content)
- Randomizes sequence order

**Dinucleotide Shuffling:**
- Preserves dinucleotide frequencies
- Maintains some local sequence properties

**Tool Used:** [BiasAway](https://biasaway.readthedocs.io/en/latest/)
- Web server or command-line tool
- Generates background sequences for motif analysis
- Maintains statistical properties of original sequences

### Step 3: Combine and Format

All sequences combined into final dataset:
```
150 original sequences + 150 RC
150 mononuc shuffled + 150 RC  
150 dinuc shuffled + 150 RC
─────────────────────────────
900 total sequences
```

### Prerequisites
```bash
pip install pandas matplotlib seaborn pyfaidx pysam biopython numpy
```

### Required Files
1. ATAC-seq peaks: `ENCFF121CAA.bed`
2. Reference genome: `hg38.fa` (and index)
3. BiasAway tool (web or command-line)

### Generation Steps

**Step 1: Extract original sequences**
```python
# Uncomment in the generation script
create_fasta()
```
This creates `sequences.fasta` with 150 random 100kb peak-centered sequences.

**Step 2: Generate shuffled sequences**

Use BiasAway to create control sequences:

**Via Web Server:**
1. Go to https://biasaway.uio.no/
2. Upload `sequences.fasta`
3. Select "K-mer shuffling" with kmer "1" for `mononucleotide_shuffled.fa` and kmer "2" for `dinucleotide_shuffled.fa`
4. Download results

**Step 3: Create final dataset**
```python
# Uncomment in the generation script
create_evaluator_datafile()
```
This combines all sequences and creates `all_consistency_data.csv`.

