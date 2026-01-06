#!/usr/bin/env python3
"""
Script to compare sequences from two different paths.

For each TCR with the same TCR_id in both paths, compares nucleotide and amino acid
sequences, creates aligned visualizations, and reports synonymous and non-synonymous
mutations with locations.
"""

import os
import sys
import argparse
from datetime import datetime
from pathlib import Path

# Genetic code for translation
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


def read_fasta(file_path):
    """Read FASTA file and return only the first sequence (uppercase, no gaps)."""
    if not os.path.exists(file_path):
        return None
    sequence = ''
    first_sequence_started = False
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if first_sequence_started:
                    break
                first_sequence_started = True
            elif first_sequence_started:
                sequence += line.upper()
    return sequence.replace('-', '').replace(' ', '')


def translate(seq):
    """Translate nucleotide sequence to amino acid."""
    seq = seq.upper()
    aa_seq = ''
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if len(codon) == 3:
            aa_seq += GENETIC_CODE.get(codon, 'X')
    return aa_seq


def get_codon_at_position(seq, pos):
    """Get the codon that contains the given position."""
    codon_start = (pos // 3) * 3
    if codon_start + 3 <= len(seq):
        return seq[codon_start:codon_start + 3]
    return None


def is_synonymous_mutation(seq1, seq2, pos):
    """
    Check if a mutation at position pos is synonymous (silent).
    Returns (is_synonymous, codon1, codon2, aa1, aa2)
    """
    codon1 = get_codon_at_position(seq1, pos)
    codon2 = get_codon_at_position(seq2, pos)
    
    if not codon1 or not codon2:
        return None, codon1, codon2, None, None
    
    aa1 = GENETIC_CODE.get(codon1, 'X')
    aa2 = GENETIC_CODE.get(codon2, 'X')
    
    is_syn = (aa1 == aa2) and (aa1 != 'X' and aa2 != 'X')
    return is_syn, codon1, codon2, aa1, aa2


def find_all_mutations(seq1_nt, seq2_nt):
    """
    Find all mutations and classify them as synonymous or non-synonymous.
    Returns: (synonymous_mutations, non_synonymous_mutations)
    """
    synonymous = []
    non_synonymous = []
    
    min_len = min(len(seq1_nt), len(seq2_nt))
    
    # Check each position for differences
    for pos in range(min_len):
        if seq1_nt[pos] != seq2_nt[pos]:
            is_syn, codon1, codon2, aa1, aa2 = is_synonymous_mutation(seq1_nt, seq2_nt, pos)
            
            if is_syn is None:
                continue
            
            mutation_info = {
                'position': pos,
                'nucleotide': (seq1_nt[pos], seq2_nt[pos]),
                'codon': (codon1, codon2),
                'amino_acid': (aa1, aa2)
            }
            
            if is_syn:
                synonymous.append(mutation_info)
            else:
                non_synonymous.append(mutation_info)
    
    return synonymous, non_synonymous


def align_sequences(seq1, seq2):
    """
    Perform pairwise alignment of two sequences (simple gap-free alignment for now).
    Returns aligned sequences with gaps as hyphens.
    """
    # For now, assume sequences are already aligned (no gaps)
    # If lengths differ, pad the shorter one with hyphens
    max_len = max(len(seq1), len(seq2))
    aligned_seq1 = seq1 + '-' * (max_len - len(seq1))
    aligned_seq2 = seq2 + '-' * (max_len - len(seq2))
    return aligned_seq1, aligned_seq2


def calculate_similarity(seq1, seq2):
    """Calculate similarity percentage between two aligned sequences."""
    if len(seq1) != len(seq2):
        return 0.0, 0, 0
    
    matches = 0
    total = len(seq1)
    
    for i in range(total):
        if seq1[i] != '-' and seq2[i] != '-' and seq1[i] == seq2[i]:
            matches += 1
    
    similarity = (matches / total * 100) if total > 0 else 0.0
    return similarity, matches, total


def format_alignment_chris_style(seq1_nt, seq2_nt, seq1_aa, seq2_aa, seq1_name="Seq_1", seq2_name="Seq_2", line_length=60):
    """
    Format alignment in the style shown in the Chris document:
    - Title with sequence names
    - Amino acid sequence for Seq1 (above)
    - Nucleotide sequence for Seq1 with positions
    - Alignment line with | for matches (aligned with nucleotides)
    - Nucleotide sequence for Seq2 with positions (gaps as -)
    - Amino acid sequence for Seq2 (below)
    """
    # Align sequences
    aligned_seq1_nt, aligned_seq2_nt = align_sequences(seq1_nt, seq2_nt)
    
    # Calculate similarity
    similarity, matches, total = calculate_similarity(aligned_seq1_nt, aligned_seq2_nt)
    
    lines = []
    # Title with sequence names
    lines.append(f"Similarity : {matches}/{total} ({similarity:.2f} %).")
    lines.append(f"Sequence 1: {seq1_name}")
    lines.append(f"Sequence 2: {seq2_name}")
    lines.append("")
    
    max_len = len(aligned_seq1_nt)
    
    for i in range(0, max_len, line_length):
        end = min(i + line_length, max_len)
        
        # Extract chunks
        seq1_nt_chunk = aligned_seq1_nt[i:end]
        seq2_nt_chunk = aligned_seq2_nt[i:end]
        
        # Format with position numbers on the right
        seq1_start = i + 1  # 1-based positions
        seq1_end = end
        seq2_start = i + 1
        seq2_end = end
        
        # Adjust seq2 positions if there are leading/trailing gaps
        seq2_actual_start = seq2_start
        seq2_actual_end = seq2_end
        
        # Count leading gaps in seq2
        leading_gaps = 0
        for char in seq2_nt_chunk:
            if char == '-':
                leading_gaps += 1
            else:
                break
        
        # Count trailing gaps in seq2
        trailing_gaps = 0
        for char in reversed(seq2_nt_chunk):
            if char == '-':
                trailing_gaps += 1
            else:
                break
        
        # Adjust positions for gaps
        if leading_gaps > 0:
            seq2_actual_start = max(1, seq2_start - leading_gaps)
        if trailing_gaps > 0:
            seq2_actual_end = max(seq2_actual_start, seq2_end - trailing_gaps)
        
        # Format nucleotide sequences with position numbers
        seq1_nt_line = f"  {seq1_nt_chunk} {seq1_start}-{seq1_end}"
        seq2_nt_line = f"  {seq2_nt_chunk} {seq2_actual_start}-{seq2_actual_end}"
        
        # Create alignment line with | for matches, aligned with nucleotides
        # The alignment line must match the exact length of the nucleotide sequence part
        alignment_line = ''
        for j in range(len(seq1_nt_chunk)):
            pos = i + j
            if pos < len(aligned_seq1_nt) and pos < len(aligned_seq2_nt):
                if aligned_seq1_nt[pos] != '-' and aligned_seq2_nt[pos] != '-' and aligned_seq1_nt[pos] == aligned_seq2_nt[pos]:
                    alignment_line += '|'
                else:
                    alignment_line += ' '
            else:
                alignment_line += ' '
        
        # Pad alignment line to match the nucleotide line format (with position numbers)
        # The alignment line should align with the nucleotide sequence, not the position numbers
        alignment_line_padded = f"  {alignment_line}"
        
        # Calculate amino acid sequences for this chunk
        # Translate nucleotide chunks to amino acids, handling gaps
        seq1_aa_display = ''
        seq2_aa_display = ''
        
        # Process seq1 nucleotide chunk codon by codon
        for j in range(0, len(seq1_nt_chunk) - 2, 3):
            codon1 = seq1_nt_chunk[j:j+3]
            if len(codon1) == 3 and '-' not in codon1:
                aa = GENETIC_CODE.get(codon1, 'X')
                seq1_aa_display += aa + '  '  # Space between amino acids
            elif len(codon1) == 3:
                seq1_aa_display += '   '  # Gap in nucleotide sequence
        
        # Process seq2 nucleotide chunk codon by codon
        for j in range(0, len(seq2_nt_chunk) - 2, 3):
            codon2 = seq2_nt_chunk[j:j+3] if j < len(seq2_nt_chunk) else '---'
            if len(codon2) == 3 and '-' not in codon2:
                aa = GENETIC_CODE.get(codon2, 'X')
                seq2_aa_display += aa + '  '  # Space between amino acids
            elif len(codon2) == 3:
                seq2_aa_display += '   '  # Gap in nucleotide sequence
        
        # Trim trailing spaces
        seq1_aa_display = seq1_aa_display.rstrip()
        seq2_aa_display = seq2_aa_display.rstrip()
        
        # Write amino acid sequence for Seq1
        lines.append(f"  {seq1_aa_display}")
        
        # Write nucleotide sequence for Seq1 (without file name)
        lines.append(seq1_nt_line)
        
        # Write alignment line (aligned with nucleotides)
        lines.append(alignment_line_padded)
        
        # Write nucleotide sequence for Seq2 (without file name)
        lines.append(seq2_nt_line)
        
        # Write amino acid sequence for Seq2
        lines.append(f"  {seq2_aa_display}")
        
        lines.append("")
    
    return "\n".join(lines)


def format_alignment(seq1, seq2, seq1_name="Seq1", seq2_name="Seq2", line_length=60, start_pos=0):
    """
    Format two sequences side-by-side with differences highlighted.
    Returns formatted string in a clear alignment format.
    (Kept for backward compatibility, but we'll use format_alignment_chris_style for the report)
    """
    lines = []
    max_len = max(len(seq1), len(seq2))
    
    for i in range(0, max_len, line_length):
        end = min(i + line_length, max_len)
        seq1_chunk = seq1[i:end] if i < len(seq1) else ''
        seq2_chunk = seq2[i:end] if i < len(seq2) else ''
        
        # Pad shorter sequence
        if len(seq1_chunk) < len(seq2_chunk):
            seq1_chunk += ' ' * (len(seq2_chunk) - len(seq1_chunk))
        elif len(seq2_chunk) < len(seq1_chunk):
            seq2_chunk += ' ' * (len(seq1_chunk) - len(seq2_chunk))
        
        # Create difference indicator line
        diff_line = ''
        for j in range(len(seq1_chunk)):
            pos = i + j
            if pos < min(len(seq1), len(seq2)):
                if seq1[pos] != seq2[pos]:
                    diff_line += '^'
                else:
                    diff_line += ' '
            else:
                diff_line += '|'  # Indicate length difference
        
        # Format with labels and position numbers
        pos_label = f"{start_pos + i:>6d}"
        lines.append(f"{pos_label}  {seq1_name:20s} {seq1_chunk}")
        lines.append(f"{' ' * 8}{seq2_name:20s} {seq2_chunk}")
        lines.append(f"{' ' * 8}{'':20s} {diff_line}")
        lines.append("")
    
    return "\n".join(lines)


def find_fasta_files(directory, filename_pattern=None):
    """
    Find all FASTA files in a directory structure.
    Returns dict: {tcr_id: file_path}
    
    Args:
        directory: Directory to search
        filename_pattern: Optional pattern to match (e.g., "TRB_P2A_TRA_Seq.fasta")
                         If None, finds first .fasta file in each TCR directory
    """
    fasta_files = {}
    
    if not os.path.exists(directory):
        return fasta_files
    
    # Check if directory contains subdirectories (TCR folders)
    items = os.listdir(directory)
    has_subdirs = any(os.path.isdir(os.path.join(directory, item)) for item in items if item != '.DS_Store')
    
    if has_subdirs:
        # Structure: directory/tcr_id/*.fasta
        for item in items:
            if item == '.DS_Store':
                continue
            tcr_dir = os.path.join(directory, item)
            if os.path.isdir(tcr_dir):
                # If filename_pattern is specified, look for exact match first
                if filename_pattern:
                    pattern_path = os.path.join(tcr_dir, filename_pattern)
                    if os.path.exists(pattern_path):
                        fasta_files[item] = pattern_path
                        continue
                
                # Otherwise, find first .fasta file in this directory
                for file in sorted(os.listdir(tcr_dir)):
                    if file.endswith('.fasta'):
                        fasta_files[item] = os.path.join(tcr_dir, file)
                        break
    else:
        # Structure: directory/*.fasta (files directly in directory)
        for file in items:
            if file.endswith('.fasta'):
                if filename_pattern and file != filename_pattern:
                    continue
                tcr_id = os.path.splitext(file)[0]
                fasta_files[tcr_id] = os.path.join(directory, file)
    
    return fasta_files


def compare_tcr_sequences(tcr_id, seq1_path, seq2_path):
    """Compare sequences for a single TCR."""
    seq1_nt = read_fasta(seq1_path)
    seq2_nt = read_fasta(seq2_path)
    
    # Extract file names from paths
    seq1_filename = os.path.basename(seq1_path)
    seq2_filename = os.path.basename(seq2_path)
    
    if seq1_nt is None:
        return {
            'tcr_id': tcr_id,
            'status': 'missing_seq1',
            'seq1_nt': None,
            'seq2_nt': seq2_nt,
            'seq1_filename': seq1_filename,
            'seq2_filename': seq2_filename,
            'synonymous': [],
            'non_synonymous': []
        }
    
    if seq2_nt is None:
        return {
            'tcr_id': tcr_id,
            'status': 'missing_seq2',
            'seq1_nt': seq1_nt,
            'seq2_nt': None,
            'seq1_filename': seq1_filename,
            'seq2_filename': seq2_filename,
            'synonymous': [],
            'non_synonymous': []
        }
    
    # Translate to amino acids
    seq1_aa = translate(seq1_nt)
    seq2_aa = translate(seq2_nt)
    
    # Find mutations
    synonymous, non_synonymous = find_all_mutations(seq1_nt, seq2_nt)
    
    status = 'identical' if not synonymous and not non_synonymous else 'different'
    
    return {
        'tcr_id': tcr_id,
        'status': status,
        'seq1_nt': seq1_nt,
        'seq2_nt': seq2_nt,
        'seq1_aa': seq1_aa,
        'seq2_aa': seq2_aa,
        'seq1_filename': seq1_filename,
        'seq2_filename': seq2_filename,
        'synonymous': synonymous,
        'non_synonymous': non_synonymous
    }


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='Compare sequences from two paths and generate alignment report'
    )
    parser.add_argument('seq1_path', help='Path to first set of sequences')
    parser.add_argument('seq2_path', help='Path to second set of sequences')
    parser.add_argument('output_name', nargs='?', default=None,
                       help='Optional output filename (default: seq1_folder_seq2_folder)')
    parser.add_argument('--seq1_file', default=None,
                       help='Optional: specific FASTA filename to look for in seq1_path (e.g., TRB_P2A_TRA_Seq.fasta)')
    parser.add_argument('--seq2_file', default=None,
                       help='Optional: specific FASTA filename to look for in seq2_path (e.g., TRB_P2A_TRA_original_contig.fasta)')
    
    args = parser.parse_args()
    
    seq1_path = args.seq1_path
    seq2_path = args.seq2_path
    
    # Determine output filename
    if args.output_name:
        output_name = args.output_name
    else:
        seq1_folder = os.path.basename(os.path.normpath(seq1_path))
        seq2_folder = os.path.basename(os.path.normpath(seq2_path))
        output_name = f"{seq1_folder}_{seq2_folder}"
    
    if not output_name.endswith('.txt'):
        output_name += '.txt'
    
    # Add alignment_ prefix to filename
    if not output_name.startswith('alignment_'):
        output_name = f"alignment_{output_name}"
    
    output_dir = "Test_result"
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, output_name)
    
    print("=" * 80)
    print("TCR Sequence Comparison")
    print("=" * 80)
    print(f"Sequence 1 path: {seq1_path}")
    print(f"Sequence 2 path: {seq2_path}")
    print(f"Output file: {log_file}")
    print()
    
    # Find all FASTA files in both paths
    print("Scanning for FASTA files...")
    if args.seq1_file:
        print(f"  Looking for '{args.seq1_file}' in seq1_path")
    if args.seq2_file:
        print(f"  Looking for '{args.seq2_file}' in seq2_path")
    seq1_files = find_fasta_files(seq1_path, args.seq1_file)
    seq2_files = find_fasta_files(seq2_path, args.seq2_file)
    
    print(f"Found {len(seq1_files)} FASTA files in seq1_path")
    print(f"Found {len(seq2_files)} FASTA files in seq2_path")
    
    # Find common TCR IDs
    common_tcrs = sorted(set(seq1_files.keys()) & set(seq2_files.keys()))
    print(f"Found {len(common_tcrs)} common TCR IDs to compare")
    print()
    
    if not common_tcrs:
        print("No common TCR IDs found to compare.")
        return
    
    # Compare each TCR
    results = []
    identical_count = 0
    different_count = 0
    missing_count = 0
    
    for tcr_id in common_tcrs:
        print(f"Comparing {tcr_id}...")
        result = compare_tcr_sequences(tcr_id, seq1_files[tcr_id], seq2_files[tcr_id])
        results.append(result)
        
        if result['status'] == 'identical':
            identical_count += 1
            print(f"  ✓ Identical")
        elif result['status'] == 'different':
            different_count += 1
            syn_count = len(result['synonymous'])
            non_syn_count = len(result['non_synonymous'])
            print(f"  ✗ Different: {syn_count} synonymous, {non_syn_count} non-synonymous mutations")
        else:
            missing_count += 1
            print(f"  ⚠ {result['status']}")
    
    print()
    print("=" * 80)
    print("Generating comparison report...")
    
    # Generate report
    with open(log_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("TCR Sequence Comparison Report\n")
        f.write("=" * 80 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"\n")
        f.write(f"Sequence 1 path: {seq1_path}\n")
        f.write(f"Sequence 2 path: {seq2_path}\n")
        f.write(f"\n")
        f.write(f"Total TCRs compared: {len(common_tcrs)}\n")
        f.write(f"Identical: {identical_count}\n")
        f.write(f"Different: {different_count}\n")
        f.write(f"Missing files: {missing_count}\n")
        f.write(f"\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
        
        # Write detailed results for each TCR
        for result in results:
            tcr_id = result['tcr_id']
            f.write(f"\n{'=' * 80}\n")
            f.write(f"TCR ID: {tcr_id}\n")
            f.write(f"{'=' * 80}\n")
            
            if result['status'] in ['missing_seq1', 'missing_seq2']:
                f.write(f"Status: {result['status']}\n")
                continue
            
            f.write(f"Status: {result['status'].upper()}\n")
            f.write(f"\n")
            
            # Sequence lengths
            seq1_nt = result['seq1_nt']
            seq2_nt = result['seq2_nt']
            seq1_aa = result['seq1_aa']
            seq2_aa = result['seq2_aa']
            
            f.write("Sequence Lengths:\n")
            f.write(f"  Nucleotide - Seq1: {len(seq1_nt)} bp, Seq2: {len(seq2_nt)} bp\n")
            f.write(f"  Amino acid - Seq1: {len(seq1_aa)} aa, Seq2: {len(seq2_aa)} aa\n")
            f.write(f"\n")
            
            # Mutation summary
            syn_count = len(result['synonymous'])
            non_syn_count = len(result['non_synonymous'])
            f.write("Mutation Summary:\n")
            f.write(f"  Synonymous mutations: {syn_count}\n")
            f.write(f"  Non-synonymous mutations: {non_syn_count}\n")
            f.write(f"\n")
            
            # Detailed mutation lists
            if result['synonymous']:
                f.write("Synonymous Mutations (Silent):\n")
                f.write("  Position | Nucleotide | Codon      | Amino Acid\n")
                f.write("  " + "-" * 60 + "\n")
                for mut in result['synonymous']:
                    pos = mut['position']
                    nt1, nt2 = mut['nucleotide']
                    codon1, codon2 = mut['codon']
                    aa1, aa2 = mut['amino_acid']
                    f.write(f"  {pos:8d} | {nt1} -> {nt2:3s} | {codon1} -> {codon2:3s} | {aa1} -> {aa2:3s} (synonymous)\n")
                f.write(f"\n")
            
            if result['non_synonymous']:
                f.write("Non-Synonymous Mutations:\n")
                f.write("  Position | Nucleotide | Codon      | Amino Acid\n")
                f.write("  " + "-" * 60 + "\n")
                for mut in result['non_synonymous']:
                    pos = mut['position']
                    nt1, nt2 = mut['nucleotide']
                    codon1, codon2 = mut['codon']
                    aa1, aa2 = mut['amino_acid']
                    f.write(f"  {pos:8d} | {nt1} -> {nt2:3s} | {codon1} -> {codon2:3s} | {aa1} -> {aa2:3s}\n")
                f.write(f"\n")
            
            # Aligned sequences visualization in Chris document style
            f.write("=" * 80 + "\n")
            f.write("ALIGNED SEQUENCES (Nucleotide and Amino Acid):\n")
            f.write("=" * 80 + "\n")
            f.write("\n")
            # Use actual file names
            seq1_name = result.get('seq1_filename', 'Seq_1')
            seq2_name = result.get('seq2_filename', 'Seq_2')
            f.write(format_alignment_chris_style(
                seq1_nt, seq2_nt, 
                seq1_aa, seq2_aa,
                seq1_name=seq1_name, 
                seq2_name=seq2_name,
                line_length=60
            ))
            f.write("\n")
        
        # Write comprehensive summary of all mutations at the end
        f.write("\n")
        f.write("=" * 80 + "\n")
        f.write("COMPREHENSIVE MUTATION SUMMARY\n")
        f.write("=" * 80 + "\n")
        f.write("\n")
        
        # Collect all mutations by type
        all_synonymous = []
        all_non_synonymous = []
        
        for result in results:
            if result['status'] == 'different':
                tcr_id = result['tcr_id']
                for mut in result['synonymous']:
                    mut_copy = mut.copy()
                    mut_copy['tcr_id'] = tcr_id
                    all_synonymous.append(mut_copy)
                for mut in result['non_synonymous']:
                    mut_copy = mut.copy()
                    mut_copy['tcr_id'] = tcr_id
                    all_non_synonymous.append(mut_copy)
        
        # Write non-synonymous mutations summary
        f.write("NON-SYNONYMOUS MUTATIONS:\n")
        f.write("-" * 80 + "\n")
        if all_non_synonymous:
            f.write("TCR ID                                    | Position | Nucleotide | Codon      | Amino Acid\n")
            f.write("-" * 80 + "\n")
            for mut in sorted(all_non_synonymous, key=lambda x: (x['tcr_id'], x['position'])):
                tcr_id = mut['tcr_id']
                pos = mut['position']
                nt1, nt2 = mut['nucleotide']
                codon1, codon2 = mut['codon']
                aa1, aa2 = mut['amino_acid']
                f.write(f"{tcr_id:40s} | {pos:8d} | {nt1} -> {nt2:3s} | {codon1} -> {codon2:3s} | {aa1} -> {aa2:3s}\n")
        else:
            f.write("No non-synonymous mutations found.\n")
        f.write("\n")
        
        # Write synonymous mutations summary
        f.write("SYNONYMOUS MUTATIONS (SILENT):\n")
        f.write("-" * 80 + "\n")
        if all_synonymous:
            f.write("TCR ID                                    | Position | Nucleotide | Codon      | Amino Acid\n")
            f.write("-" * 80 + "\n")
            for mut in sorted(all_synonymous, key=lambda x: (x['tcr_id'], x['position'])):
                tcr_id = mut['tcr_id']
                pos = mut['position']
                nt1, nt2 = mut['nucleotide']
                codon1, codon2 = mut['codon']
                aa1, aa2 = mut['amino_acid']
                f.write(f"{tcr_id:40s} | {pos:8d} | {nt1} -> {nt2:3s} | {codon1} -> {codon2:3s} | {aa1} -> {aa2:3s} (synonymous)\n")
        else:
            f.write("No synonymous mutations found.\n")
        f.write("\n")
        f.write("=" * 80 + "\n")
    
    print(f"Comparison report saved to: {log_file}")
    print()
    print("=" * 80)
    print("Summary:")
    print(f"  Total TCRs compared: {len(common_tcrs)}")
    print(f"  Identical: {identical_count}")
    print(f"  Different: {different_count}")
    print(f"  Missing files: {missing_count}")
    print("=" * 80)


if __name__ == "__main__":
    main()
