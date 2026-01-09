# TCRAssembler

A comprehensive Python package for assembling T-cell receptor (TCR) sequences from IMGT templates or full contig sequences, with automatic restriction site removal, sequence validation, and formatting for molecular cloning.

## Table of Contents

- [Overview](#overview)
- [Assembly Logic](#assembly-logic)
- [Templates Used](#templates-used)
- [Input File Requirements](#input-file-requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Scripts and Functions](#scripts-and-functions)
- [Output Structure](#output-structure)
- [License](#license)
- [Citation](#citation)

## Overview

This package provides two complementary methods for assembling TCR sequences:

1. **IMGT Template-Based Assembly**: Assembles TCR sequences from IMGT gene templates and CDR3 sequences
2. **Contig-Based Assembly**: Assembles TCR sequences directly from full contig sequences

Both methods include:
- Automatic restriction site removal (AscI, PstI) via silent mutations
- Sequence validation (frame checking, stop codon detection)
- Proper sequence formatting (uppercase for restriction sites, CDR3s, P2A linker)
- TRB-P2A-TRA final assembly with Kozak sequence

## Assembly Logic

### IMGT Template-Based Assembly (`assemble_IMGT_tcr.py`)

This method assembles TCR sequences from IMGT gene templates:

1. **TRB Assembly**:
   - Aligns TRBV template with TRB CDR3 sequence (finds overlap)
   - Merges TRBJ sequence with CDR3 (finds J overlap)
   - Optionally replaces CDR1/CDR2 if provided (creates IMGT and Seq variants)
   - Removes restriction sites (AscI: `GGCGCGCC`, PstI: `CTGCAG`) via silent mutations
   - Adds complete TRBC constant region
   - Validates frame and checks for stop codons

2. **TRA Assembly**:
   - Aligns TRAV template with TRA CDR3 sequence
   - Merges TRAJ sequence with CDR3
   - Keeps "at" ending from TRAJ if present
   - Removes restriction sites via silent mutations
   - Attaches TRAC constant region (with frame adjustment if needed)
   - Validates frame and checks for stop codons

3. **Final Assembly**:
   - Combines TRB and TRA with P2A linker: `AscI + Kozak + TRB + P2A + TRA + PstI`
   - Formats sequence with proper casing

### Contig-Based Assembly (`assemble_from_contig.py`)

This method uses full contig sequences directly:

1. **TRB Contig Processing**:
   - Aligns full TRB contig to TRBV template to find V gene start position
   - Extracts contig sequence from V gene start
   - Finds J-C junction and replaces incomplete TRBC with complete template
   - Removes restriction sites via silent mutations

2. **TRA Contig Processing**:
   - Aligns full TRA contig to TRAV template to find V gene start position
   - Extracts contig sequence from V gene start
   - Keeps TRAJ-TRAC junction from contig (preserves natural frame)
   - Ensures TRAC length matches template
   - Removes restriction sites via silent mutations

3. **Final Assembly**:
   - Combines TRB and TRA with P2A linker
   - Formats sequence with proper casing
   - Validates sequences (translation and stop codon check)

### Restriction Site Removal

Both methods use a sophisticated silent mutation strategy:

1. **Individual codon mutations**: Tries to mutate one codon at a time
2. **Paired codon mutations**: Mutates two adjacent codons simultaneously
3. **Overlapping codon mutations**: Handles sites spanning multiple codons

All mutations preserve the amino acid sequence (synonymous mutations).

### Sequence Formatting

The final assembled sequences follow specific casing rules:

- **Uppercase**:
  - ASCI_SITE: `GGCGCGCC`
  - KOZAK_SEQ: `GGCCACC`
  - P2A sequence: Full P2A linker
  - PSTI_SITE: `CTGCAG`
  - CDR3 sequences (TRA and TRB)
  - Mutated codons (from restriction site removal)

- **Lowercase**: All other sequence regions

## Templates Used

The package requires IMGT gene templates downloaded from the [IMGT database](https://www.imgt.org/). Templates are organized in the `Templates/` directory:

### Required Template Files

- **TRBV templates** (`Templates/TRBV/*.fasta`): TRB variable gene templates
  - Format: `fastaNConstructs` (includes L-PART1+V-EXON)
  - Example: `TRBV1-1.fasta`, `TRBV2-1.fasta`

- **TRBJ templates** (`Templates/TRBJ/*.fasta`): TRB joining gene templates
  - Format: `fastaNAllP` (J-REGION only)
  - Example: `TRBJ1-1.fasta`, `TRBJ2-1.fasta`

- **TRBC templates** (`Templates/TRBC/*.fasta`): TRB constant region templates
  - Example: `TRBC1.fasta`, `TRBC2.fasta`

- **TRAV templates** (`Templates/TRAV/*.fasta`): TRA variable gene templates
  - Format: `fastaNConstructs` (includes L-PART1+V-EXON)
  - Example: `TRAV1-1.fasta`, `TRAV26-2.fasta`
  - Includes /DV variants: `TRAV38-2_DV8.fasta`

- **TRAJ templates** (`Templates/TRAJ/*.fasta`): TRA joining gene templates
  - Format: `fastaNAllP` (J-REGION only)
  - Example: `TRAJ1.fasta`, `TRAJ45.fasta`

- **TRAC template** (`Templates/TRAC/TRAC_head.fasta`): TRA constant region
  - Length: 18 nucleotides

- **P2A linker** (`Templates/codon_optimized_P2A.fasta`): P2A self-cleaving peptide sequence

### Template Download

Templates can be downloaded using the provided script:

```bash
python Script/download_imgt_sequences.py
```

This script downloads all required templates from IMGT and organizes them into the appropriate folders.

## Input File Requirements

### For IMGT Template-Based Assembly

Required CSV columns:
- `tcr_id`: Unique TCR identifier
- `TRA_v_gene`: TRA variable gene name (e.g., "TRAV1-1", "TRAV26-2")
- `TRA_j_gene`: TRA joining gene name (e.g., "TRAJ1", "TRAJ45")
- `TRA_cdr3_nt`: TRA CDR3 nucleotide sequence
- `TRB_v_gene`: TRB variable gene name (e.g., "TRBV1-1", "TRBV5-4")
- `TRB_j_gene`: TRB joining gene name (e.g., "TRBJ1-1", "TRBJ2-1")
- `TRB_c_gene`: TRB constant gene name ("TRBC1" or "TRBC2")
- `TRB_cdr3_nt`: TRB CDR3 nucleotide sequence

Optional columns (for CDR1/CDR2 replacement):
- `TRA_CDR1_nn`, `TRA_CDR1_aa`: TRA CDR1 nucleotide and amino acid sequences
- `TRA_CDR2_nn`, `TRA_CDR2_aa`: TRA CDR2 nucleotide and amino acid sequences
- `TRB_CDR1_nn`, `TRB_CDR1_aa`: TRB CDR1 nucleotide and amino acid sequences
- `TRB_CDR2_nn`, `TRB_CDR2_aa`: TRB CDR2 nucleotide and amino acid sequences

### For Contig-Based Assembly

Required CSV columns:
- All columns from IMGT template-based assembly, plus:
- `TRA_full_contig`: Full TRA contig sequence from sequencing
- `TRB_full_contig`: Full TRB contig sequence from sequencing
- `TRA_contig_id`: TRA contig identifier (optional)
- `TRB_contig_id`: TRB contig identifier (optional)

### Example Input File

See `Input/sample_input.csv` for a complete example.

## Installation

### Prerequisites

- Python 3.6 or higher
- pip (Python package manager)

### Install Dependencies

```bash
pip install -r Script/requirements.txt
```

Required packages:
- `pandas >= 1.5.0`
- `requests >= 2.31.0` (for template downloading)
- `beautifulsoup4 >= 4.12.0` (for template downloading)
- `lxml >= 4.9.0` (for template downloading)

## Usage

### Master Script (Recommended)

The master script automatically selects the appropriate assembly method:

```bash
# Use default input file
python Script/assemble_tcr.py

# Specify input file
python Script/assemble_tcr.py -i Input/my_data.csv

# Specify input file and output directory
python Script/assemble_tcr.py -i Input/my_data.csv -o Output/my_results

# Process only first 10 TCRs
python Script/assemble_tcr.py -i Input/my_data.csv -n 10

# Force specific assembly method
python Script/assemble_tcr.py -i Input/my_data.csv -m imgt      # IMGT only
python Script/assemble_tcr.py -i Input/my_data.csv -m contig    # Contig only
python Script/assemble_tcr.py -i Input/my_data.csv -m both       # Both methods
```

### Individual Scripts

You can also run the assembly scripts directly:

```bash
# IMGT template-based assembly
python Script/assemble_IMGT_tcr.py

# Contig-based assembly
python Script/assemble_from_contig.py

# Test assembly (with detailed validation)
python Script/test_assemble_IMGT_tcr.py
```

### Command-Line Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--input` | `-i` | `Input/sample_input.csv` | Path to input CSV file |
| `--output` | `-o` | `Output/{input_filename}` | Path to output directory |
| `--num-tcrs` | `-n` | `None` (all) | Number of TCRs to process |
| `--method` | `-m` | Auto-detect | Assembly method: `imgt`, `contig`, or `both` |

## Scripts and Functions

### Main Assembly Scripts

#### `assemble_tcr.py` (Master Script)
- **Purpose**: Main entry point that automatically routes to appropriate assembly method
- **Features**:
  - Auto-detects assembly method based on input file structure
  - Supports explicit method selection (`--method` flag)
  - Organizes outputs in subdirectories (`Contig_tcr/`, `IMGT_tcr/`)
  - Automatically runs test script after IMGT assembly

#### `assemble_IMGT_tcr.py`
- **Purpose**: IMGT template-based TCR sequence assembly
- **Features**:
  - Assembles from V, J, C gene templates and CDR3 sequences
  - Supports CDR1/CDR2 replacement (creates IMGT and Seq variants)
  - Restriction site removal via silent mutations
  - Sequence validation (frame, stop codons)
  - Outputs multiple FASTA files per TCR

#### `assemble_from_contig.py`
- **Purpose**: Contig-based TCR sequence assembly
- **Features**:
  - Uses full contig sequences directly
  - Aligns contigs to V gene templates to find start positions
  - Completes incomplete constant regions
  - Preserves natural TRAJ-TRAC junctions
  - Comprehensive validation with translation and stop codon checks
  - Generates validation logs

#### `test_assemble_IMGT_tcr.py`
- **Purpose**: Test version with detailed validation output
- **Features**:
  - Same functionality as `assemble_IMGT_tcr.py`
  - Detailed logging to file
  - Step-by-step assembly process output
  - Useful for debugging and validation

### Utility Scripts

#### `download_imgt_sequences.py`
- **Purpose**: Download IMGT gene templates
- **Features**:
  - Downloads TRAV, TRBV, TRAJ, TRBJ templates from IMGT database
  - Organizes templates into appropriate folders
  - Handles /DV gene variants
  - Skips already downloaded templates
  - Respectful rate limiting

#### `compare_sequences.py`
- **Purpose**: Compare sequences from different assembly methods
- **Features**:
  - Compares IMGT vs contig-based assemblies
  - Nucleotide and amino acid alignment
  - Mutation analysis (synonymous/non-synonymous)
  - Detailed alignment visualization
  - Comprehensive mutation summary

### Key Functions

#### Sequence Assembly Functions
- `assemble_trb_sequence()`: Assembles TRB sequence from templates/CDR3
- `assemble_tra_sequence()`: Assembles TRA sequence from templates/CDR3
- `assemble_trb_p2a_tra()`: Combines TRB and TRA with P2A linker

#### Restriction Site Functions
- `check_restriction_sites()`: Finds all occurrences of restriction sites
- `make_silent_mutation()`: Makes silent mutations to remove sites
- `remove_restriction_sites()`: Removes all restriction sites from sequence

#### Alignment Functions
- `find_overlap()`: Finds overlap between two sequences
- `merge_sequences()`: Merges sequences with overlap handling
- `align_v_to_contig()`: Aligns V gene template to contig sequence

#### Validation Functions
- `translate_sequence()`: Translates nucleotide to amino acid
- `check_stop_codons()`: Checks for premature stop codons
- `validate_frame()`: Validates sequence is in frame

#### Formatting Functions
- `format_trb_tra_sequence()`: Formats TRB/TRA sequences with proper casing
- `format_final_assembly()`: Formats final TRB-P2A-TRA assembly

## Output Structure

### IMGT Template-Based Assembly Output

```
Output/{input_filename}/IMGT_tcr/
├── {tcr_id}/
│   ├── TRB_IMGT_CDRs.fasta          # TRB with IMGT CDR1/CDR2
│   ├── TRB_Seq_CDRs.fasta           # TRB with Seq CDR1/CDR2
│   ├── TRA_IMGT_CDRs.fasta          # TRA with IMGT CDR1/CDR2
│   ├── TRA_Seq_CDRs.fasta           # TRA with Seq CDR1/CDR2
│   ├── TRB_P2A_TRA_IMGT.fasta       # Final: TRB(IMGT) + P2A + TRA(IMGT)
│   └── TRB_P2A_TRA_Seq.fasta       # Final: TRB(Seq) + P2A + TRA(Seq)
└── test_assembly_log_*.txt          # Test assembly log (if test run)
```

### Contig-Based Assembly Output

```
Output/{input_filename}/Contig_tcr/
├── {tcr_id}/
│   ├── TRB_P2A_TRA_original_contig.fasta    # Final assembled sequence
│   └── {tcr_id}_validation_log.txt         # Validation log
└── final_validation_summary.txt              # Summary of all TCRs
```

### FASTA File Format

Each FASTA file contains:
- Header: `>{tcr_id}_TRB_P2A_TRA`
- Sequence: Formatted nucleotide sequence with proper casing

## License

This software is provided for academic and research purposes. For use in scientific publications, please cite this package appropriately.

### License Type

This package is released under a permissive license for academic and research use. Users are free to:

- Use the software for research and academic purposes
- Modify the code for their own use
- Distribute the software with proper attribution

### Attribution

If you use this package in your research, please cite:

```
TCRAssembler
[Yuexin Xu/Fred Hutchinson Cancer Center]
[2026]
```

### IMGT Database

This package uses templates from the IMGT database. When using IMGT data, please acknowledge:

> "This work uses data from IMGT, the international ImMunoGeneTics information system® (http://www.imgt.org) (Lefranc et al., Nucleic Acids Res., 1999, 27, 209-212)."

## Citation

If you use this package in a scientific publication, please cite:

```bibtex
@software{tcr_assembler,
  title = {TCRAssembler: A Python Package for TCR Sequence Assembly},
  author = {[Yuexin Xu]},
  year = {2026},
  url = {https://github.com/[your-username]/TCRAssembler}
}
```

## Acknowledgments

- IMGT database for providing TCR gene templates
- Contributors and testers of this package

## Support

For issues, questions, or contributions, please open an issue on the GitHub repository.

---

**Note**: This package is designed for research use. Always validate assembled sequences experimentally before use in molecular cloning applications.

