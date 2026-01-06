#!/usr/bin/env python3
"""
Test script for TCR sequence assembly - processes first 5 TCRs with detailed validation output.
"""

import os
import re
import sys
from datetime import datetime
import pandas as pd
from pathlib import Path

# Input and output paths
INPUT_FILE = "Input/hERV/hERV_TCR_synthesize_input_122225_with_CDR12.csv"
TEMPLATES_DIR = "Templates"
OUTPUT_DIR = "Output_test"

# Number of TCRs to process (None = process all)
MAX_TCRS = 16  # Set to a number (e.g., 10) to process only the first N TCRs, or None to process all

# Restriction sites
ASCI_SITE = "GGCGCGCC"
PSTI_SITE = "CTGCAG"

# Kozak sequence
KOZAK_SEQ = "GGCCACC"

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
                    # We've already started reading a sequence, stop here
                    break
                # This is the first sequence header, start reading
                first_sequence_started = True
            elif first_sequence_started:
                # We're in the first sequence, add to sequence
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


def find_max_overlap(seq1_end, seq2_start, min_overlap=3):
    """
    Find maximum overlap between the end of seq1 and the beginning of seq2.
    Checks if the beginning of seq2 appears anywhere at the end of seq1.
    Returns: (overlap_length, position_in_seq1) where position is where overlap starts in seq1
    """
    seq1_end = seq1_end.upper()
    seq2_start = seq2_start.upper()
    
    max_overlap = min(len(seq1_end), len(seq2_start))
    overlap_len = 0
    overlap_pos = -1
    
    # Try different overlap lengths, starting from maximum
    for test_overlap in range(max_overlap, min_overlap - 1, -1):
        seq2_prefix = seq2_start[:test_overlap]
        # First check if this prefix of seq2 appears at the very end of seq1
        if seq1_end.endswith(seq2_prefix):
            overlap_len = test_overlap
            overlap_pos = len(seq1_end) - test_overlap
            break
        # Also check if seq2 prefix appears anywhere in seq1_end
        # (in case there's a match that's not at the very end)
        pos_in_seq1 = seq1_end.find(seq2_prefix)
        if pos_in_seq1 != -1:
            # Found the prefix in seq1, use this overlap
            overlap_len = test_overlap
            overlap_pos = pos_in_seq1
            break
    
    return overlap_len, overlap_pos


def align_and_merge_v_cdr3(v_seq, cdr3_nt):
    """
    Align TRBV/TRAV with CDR3 and merge.
    Aligns the end of V with the beginning of CDR3, finds max overlap,
    keeps V up to the overlap point (in lowercase), then appends the FULL CDR3 sequence.
    Returns: (merged_sequence, overlap_length, position_in_v)
    """
    v_seq = v_seq.upper()
    cdr3_nt = cdr3_nt.upper()
    
    # Find maximum overlap between end of V and beginning of CDR3
    overlap_len, overlap_pos = find_max_overlap(v_seq, cdr3_nt, min_overlap=3)
    
    if overlap_len > 0:
        # Keep V sequence up to the overlap point (in lowercase), then append FULL CDR3
        # overlap_pos is where the overlap starts in V
        merged = v_seq[:overlap_pos].lower() + cdr3_nt.upper()
        pos = overlap_pos
    else:
        # No overlap found, append full CDR3 (V in lowercase)
        pos = len(v_seq)
        merged = v_seq.lower() + cdr3_nt.upper()
    
    return merged, overlap_len, pos


def align_and_merge_cdr3_j(v_cdr3_merged, j_seq, cdr3_nt_original):
    """
    Align CDR3 (at end of merged V+CDR3) with J sequence and merge.
    CDR3 comes first, then J gene. Some head nucleotides of J will be replaced by CDR3.
    Finds overlap between end of CDR3 and any position in J, preserves FULL CDR3,
    then appends the remaining J sequence after the overlap.
    Returns: (merged_sequence, overlap_length)
    
    Args:
        v_cdr3_merged: The merged V+CDR3 sequence (CDR3 is at the end, in uppercase)
        j_seq: The J gene sequence
        cdr3_nt_original: The original CDR3 nucleotide sequence (for overlap detection)
    """
    v_cdr3_merged = v_cdr3_merged.upper()
    j_seq = j_seq.upper()
    cdr3_nt_original = cdr3_nt_original.upper()
    
    # Find overlap between end of CDR3 and any position in J sequence
    # We need to check if the end of CDR3 appears anywhere in J
    cdr3_end_len = min(len(cdr3_nt_original), len(j_seq))
    overlap_len = 0
    overlap_pos_in_j = -1
    
    # Try different overlap lengths, starting from maximum
    for test_overlap in range(cdr3_end_len, 3 - 1, -1):
        cdr3_end = cdr3_nt_original[-test_overlap:]
        # Check if this CDR3 end appears in J sequence
        pos_in_j = j_seq.find(cdr3_end)
        if pos_in_j != -1:
            # Found overlap! The overlap starts at pos_in_j in J
            overlap_len = test_overlap
            overlap_pos_in_j = pos_in_j
            break
    
    if overlap_len > 0:
        # Preserve FULL V+CDR3 merged sequence, then append J sequence starting after the overlap
        # The overlapping part of J (from start to overlap_pos_in_j + overlap_len) is replaced by CDR3
        remaining_j = j_seq[overlap_pos_in_j + overlap_len:].lower()
        merged = v_cdr3_merged + remaining_j
    else:
        # No overlap found, append full J in lowercase
        merged = v_cdr3_merged + j_seq.lower()
    
    return merged, overlap_len


def check_restriction_sites(seq):
    """Check for restriction sites and return positions (all occurrences)."""
    seq = seq.upper()
    sites = {}
    
    # Check AscI - find all occurrences
    pos = 0
    asci_positions = []
    while True:
        pos = seq.find(ASCI_SITE, pos)
        if pos == -1:
            break
        asci_positions.append(pos)
        pos += 1
    if asci_positions:
        sites['AscI'] = asci_positions
    
    # Check PstI - find all occurrences
    pos = 0
    psti_positions = []
    while True:
        pos = seq.find(PSTI_SITE, pos)
        if pos == -1:
            break
        psti_positions.append(pos)
        pos += 1
    if psti_positions:
        sites['PstI'] = psti_positions
    
    return sites


def make_silent_mutation(seq, site_pos, site_seq):
    """
    Make a silent mutation to disrupt restriction site.
    Translates the region spanning the restriction site, gets amino acids,
    and tries alternative codons that preserve amino acids but break the site.
    Tries single codon mutations first, then multiple codon mutations if needed.
    
    Args:
        seq: The sequence (uppercase)
        site_pos: Position of the restriction site
        site_seq: The restriction site sequence (e.g., "CTGCAG" for PstI)
    
    Returns: (new_sequence, mutation_info) where mutation_info is (position, old_codon, new_codon) or None
    """
    site_len = len(site_seq)
    
    # Find all codons that actually contain any part of the restriction site
    # A codon contains the restriction site if it overlaps with [site_pos, site_pos + site_len)
    overlapping_codons = []
    
    # Check all possible codon positions that could overlap with the restriction site
    # Start from the codon that contains site_pos, and check codons until we've covered the entire site
    first_codon_start = (site_pos // 3) * 3
    last_codon_start = ((site_pos + site_len - 1) // 3) * 3
    
    # Check each codon from first_codon_start to last_codon_start (inclusive)
    for codon_start in range(first_codon_start, last_codon_start + 1, 3):
        codon_end = codon_start + 3
        if codon_end > len(seq):
            break
        
        # Check if this codon overlaps with the restriction site
        if codon_end > site_pos and codon_start < site_pos + site_len:
            current_codon = seq[codon_start:codon_end]
            if len(current_codon) == 3:
                aa = GENETIC_CODE.get(current_codon, 'X')
                if aa != 'X' and aa != '*':
                    alt_codons = [c for c, a in GENETIC_CODE.items() if a == aa and c != current_codon]
                    if alt_codons:
                        overlapping_codons.append((codon_start, current_codon, alt_codons))
    
    if not overlapping_codons:
        return seq, None
    
    # Strategy 1: Try mutating each codon individually
    for codon_start_absolute, current_codon, alt_codons in overlapping_codons:
        for alt_codon in alt_codons:
            test_seq = seq[:codon_start_absolute] + alt_codon + seq[codon_start_absolute + 3:]
            
            # Check if the specific site at site_pos is gone (not just any occurrence)
            # Account for potential position shifts due to mutation
            check_start = max(0, site_pos - 3)
            check_end = min(len(test_seq), site_pos + site_len + 3)
            site_gone = site_seq not in test_seq[check_start:check_end].upper()
            
            if site_gone:
                other_site = ASCI_SITE if site_seq == PSTI_SITE else PSTI_SITE
                # Check if we're creating the other restriction site in the mutated region
                other_site_in_region = other_site in test_seq[check_start:check_end].upper()
                if not other_site_in_region:
                    return test_seq, (codon_start_absolute, current_codon, alt_codon)
    
    # Strategy 2: Try mutating pairs of codons simultaneously
    if len(overlapping_codons) >= 2:
        for i, (codon1_start, codon1_old, codon1_alts) in enumerate(overlapping_codons):
            for j, (codon2_start, codon2_old, codon2_alts) in enumerate(overlapping_codons):
                if i >= j:  # Avoid duplicates
                    continue
                for alt1 in codon1_alts:
                    for alt2 in codon2_alts:
                        # Apply both mutations
                        test_seq = seq
                        # Apply mutation 1
                        test_seq = test_seq[:codon1_start] + alt1 + test_seq[codon1_start + 3:]
                        # Apply mutation 2 (adjust position if needed)
                        if codon2_start > codon1_start:
                            test_seq = test_seq[:codon2_start] + alt2 + test_seq[codon2_start + 3:]
                        else:
                            test_seq = test_seq[:codon2_start] + alt2 + test_seq[codon2_start + 3:]
                        
                        # Check if the specific site at site_pos is gone
                        check_start = max(0, site_pos - 3)
                        check_end = min(len(test_seq), site_pos + site_len + 3)
                        site_gone = site_seq not in test_seq[check_start:check_end].upper()
                        
                        if site_gone:
                            other_site = ASCI_SITE if site_seq == PSTI_SITE else PSTI_SITE
                            other_site_in_region = other_site in test_seq[check_start:check_end].upper()
                            if not other_site_in_region:
                                # Return the first mutation (we'll track both in the calling function)
                                return test_seq, (codon1_start, codon1_old, alt1)
    
    # Strategy 3: Try mutating all overlapping codons simultaneously
    if len(overlapping_codons) >= 2:
        # Try all combinations (limit to avoid combinatorial explosion)
        import itertools
        max_combinations = 100  # Limit combinations to try
        combinations_tried = 0
        
        for combo in itertools.product(*[alts for _, _, alts in overlapping_codons]):
            if combinations_tried >= max_combinations:
                break
            combinations_tried += 1
            
            test_seq = seq
            mutations_applied = []
            
            # Apply all mutations in order (from end to start to avoid position shifts)
            sorted_codons = sorted(overlapping_codons, key=lambda x: x[0], reverse=True)
            for idx, (codon_start, codon_old, _) in enumerate(sorted_codons):
                # Find the index in original list
                orig_idx = next(i for i, (cs, _, _) in enumerate(overlapping_codons) if cs == codon_start)
                alt_codon = combo[orig_idx]
                test_seq = test_seq[:codon_start] + alt_codon + test_seq[codon_start + 3:]
                mutations_applied.append((codon_start, codon_old, alt_codon))
            
            # Check if the specific site at site_pos is gone
            check_start = max(0, site_pos - 3)
            check_end = min(len(test_seq), site_pos + site_len + 3)
            site_gone = site_seq not in test_seq[check_start:check_end].upper()
            
            if site_gone:
                other_site = ASCI_SITE if site_seq == PSTI_SITE else PSTI_SITE
                other_site_in_region = other_site in test_seq[check_start:check_end].upper()
                if not other_site_in_region:
                    # Return first mutation info
                    return test_seq, mutations_applied[-1]  # Return the last one (first position)
    
    # If no mutation worked, return original
    return seq, None


def remove_restriction_sites(seq):
    """
    Remove restriction sites by making silent mutations.
    Returns: (cleaned_sequence, mutations_list)
    mutations_list contains tuples: (site_name, site_position, site_sequence, mutation_position, old_codon, new_codon)
    """
    seq = seq.upper()
    sites = check_restriction_sites(seq)
    mutations = []
    
    max_iterations = 100
    iteration = 0
    
    while sites and iteration < max_iterations:
        iteration += 1
        # Flatten sites dictionary to list of (site_name, pos) tuples
        all_sites = []
        for site_name, positions in sites.items():
            if isinstance(positions, list):
                for pos in positions:
                    all_sites.append((site_name, pos))
            else:
                all_sites.append((site_name, positions))
        
        # Sort by position (reverse order to process from end to beginning)
        sorted_sites = sorted(all_sites, key=lambda x: x[1], reverse=True)
        
        prev_seq = seq
        for site_name, pos in sorted_sites:
            site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
            actual_site_seq = seq[pos:pos+len(site_seq)]
            print(f"    Attempting to remove {site_name} site at position {pos}: {actual_site_seq}")
            new_seq, mutation_info = make_silent_mutation(seq, pos, site_seq)
            if new_seq != seq and mutation_info:
                mut_pos, old_codon, new_codon = mutation_info
                mutations.append((site_name, pos, actual_site_seq, mut_pos, old_codon, new_codon))
                print(f"      ✓ Removed {site_name} site at position {pos}: {actual_site_seq}")
                print(f"      Mutation at position {mut_pos}: {old_codon} -> {new_codon} (silent)")
                print(f"      Amino acid preserved: {GENETIC_CODE.get(old_codon, '?')} = {GENETIC_CODE.get(new_codon, '?')}")
                seq = new_seq
                sites = check_restriction_sites(seq)
                if not sites:
                    break
            else:
                print(f"      ✗ Could not remove {site_name} site at position {pos} - no suitable silent mutation found")
        
        if seq == prev_seq:
            break
    
    return seq, mutations


def align_and_merge_with_c(chain_seq, c_seq):
    """
    Attach C sequence to chain sequence.
    C sequence should be in uppercase in the final assembly.
    No overlap detection - just append the full C sequence.
    """
    c_seq = c_seq.upper()
    
    # Simply append the full C sequence (uppercase) to the chain sequence
    merged = chain_seq + c_seq.upper()
    
    return merged, 0  # Return 0 overlap since we're not doing overlap detection


def get_template_path(gene_name, gene_type):
    """Get path to template FASTA file."""
    safe_name = gene_name.replace("/", "_")
    
    if gene_type == "TRBV":
        return os.path.join(TEMPLATES_DIR, "TRBV", f"{safe_name}.fasta")
    elif gene_type == "TRBJ":
        return os.path.join(TEMPLATES_DIR, "TRBJ", f"{safe_name}.fasta")
    elif gene_type == "TRBC":
        return os.path.join(TEMPLATES_DIR, "TRBC", f"{safe_name}.fasta")
    elif gene_type == "TRAV":
        return os.path.join(TEMPLATES_DIR, "TRAV", f"{safe_name}.fasta")
    elif gene_type == "TRAJ":
        return os.path.join(TEMPLATES_DIR, "TRAJ", f"{safe_name}.fasta")
    elif gene_type == "TRAC":
        return os.path.join(TEMPLATES_DIR, "TRAC", "TRAC_head.fasta")
    return None


def replace_cdr_in_sequence(seq, cdr_nn, cdr_aa, cdr_type="CDR1"):
    """Replace CDR in sequence if it differs from template."""
    if not cdr_nn or pd.isna(cdr_nn):
        return seq
    
    cdr_nn = str(cdr_nn).upper()
    seq_upper = seq.upper()
    
    # Try to find CDR in sequence (exact match)
    cdr_pos = seq_upper.find(cdr_nn)
    if cdr_pos != -1:
        # Found exact match, check if amino acid differs
        template_aa = translate(seq_upper[cdr_pos:cdr_pos+len(cdr_nn)])
        if template_aa != cdr_aa:
            # Replace with new CDR
            return seq[:cdr_pos] + cdr_nn + seq[cdr_pos+len(cdr_nn):]
        return seq
    
    return seq


def print_validation_info(label, nn_seq, aa_seq=None):
    """Print validation information for a sequence with full sequences."""
    print(f"\n  {label}:")
    print(f"    Nucleotide ({len(nn_seq)} bp):")
    # Print full sequence with line breaks every 60 characters for readability
    for i in range(0, len(nn_seq), 60):
        print(f"      {nn_seq[i:i+60]}")
    
    if aa_seq:
        print(f"    Amino acid ({len(aa_seq)} aa):")
        # Print full sequence with line breaks every 60 characters for readability
        for i in range(0, len(aa_seq), 60):
            print(f"      {aa_seq[i:i+60]}")
    else:
        aa = translate(nn_seq)
        print(f"    Amino acid ({len(aa)} aa):")
        # Print full sequence with line breaks every 60 characters for readability
        for i in range(0, len(aa), 60):
            print(f"      {aa[i:i+60]}")


def save_validation_file(content, filepath):
    """Save validation content to a text file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        f.write(content)


def assemble_trb_sequence(row, use_seq_cdrs=False, tcr_id="", output_dir=""):
    """Assemble TRB sequence with detailed validation output."""
    print(f"\n{'='*60}")
    print(f"Assembling TRB sequence ({'Seq_CDRs' if use_seq_cdrs else 'IMGT_CDRs'})")
    print(f"{'='*60}")
    
    # Get gene names
    trbv_name = row['TRB_v_gene']
    trbj_name = row['TRB_j_gene']
    trbc_name = row['TRB_c_gene']
    cdr3_nt = row['TRB_cdr3_nt']
    cdr3_aa = row['TRB_cdr3']
    
    # Load templates
    trbv_path = get_template_path(trbv_name, "TRBV")
    trbj_path = get_template_path(trbj_name, "TRBJ")
    trbc_path = get_template_path(trbc_name, "TRBC")
    
    trbv_seq = read_fasta(trbv_path)
    trbj_seq = read_fasta(trbj_path)
    trbc_seq = read_fasta(trbc_path)
    
    if not all([trbv_seq, trbj_seq, trbc_seq, cdr3_nt]):
        print("  Error: Missing templates or CDR3")
        return None, []
    
    # Print template sequences
    print("\n  Template Sequences:")
    print_validation_info(f"TRBV ({trbv_name})", trbv_seq)
    print_validation_info(f"TRBJ ({trbj_name})", trbj_seq)
    print_validation_info(f"TRBC ({trbc_name})", trbc_seq)
    print_validation_info("TRB CDR3", cdr3_nt, cdr3_aa)
    
    # Step 1: Align TRBV with CDR3
    print("\n  Step 1: Aligning TRBV with CDR3")
    trb_seq, overlap_len, pos = align_and_merge_v_cdr3(trbv_seq, cdr3_nt)
    print(f"    Overlap: {overlap_len} bp at position {pos}")
    print_validation_info("TRBV + CDR3", trb_seq)
    
    # Step 2: Align CDR3 with TRBJ
    print("\n  Step 2: Aligning CDR3 with TRBJ")
    trb_seq, j_overlap = align_and_merge_cdr3_j(trb_seq, trbj_seq, cdr3_nt)
    print(f"    J overlap: {j_overlap} bp")
    print_validation_info("TRBV + CDR3 + TRBJ", trb_seq)
    
    # Step 3: Handle CDR1 and CDR2 if available
    cdr1_nn = row.get('TRB_CDR1_nn', '')
    cdr1_aa = row.get('TRB_CDR1_aa', '')
    cdr2_nn = row.get('TRB_CDR2_nn', '')
    cdr2_aa = row.get('TRB_CDR2_aa', '')
    
    if use_seq_cdrs and cdr1_nn and cdr2_nn and pd.notna(cdr1_nn) and pd.notna(cdr2_nn):
        print("\n  Step 3: Replacing CDR1 and CDR2")
        if cdr1_nn and cdr1_aa:
            print_validation_info("TRB CDR1 (Seq)", cdr1_nn, str(cdr1_aa))
            trb_seq = replace_cdr_in_sequence(trb_seq, cdr1_nn, str(cdr1_aa), "CDR1")
        if cdr2_nn and cdr2_aa:
            print_validation_info("TRB CDR2 (Seq)", cdr2_nn, str(cdr2_aa))
            trb_seq = replace_cdr_in_sequence(trb_seq, cdr2_nn, str(cdr2_aa), "CDR2")
        print_validation_info("After CDR replacement", trb_seq)
    
    # Step 4: Remove restriction sites
    print("\n  Step 4: Removing restriction sites")
    sites_before = check_restriction_sites(trb_seq.upper())
    if sites_before:
        print(f"    Found restriction sites: {sites_before}")
    trb_seq, trb_mutations = remove_restriction_sites(trb_seq)
    trb_mutations = trb_mutations or []
    sites_after = check_restriction_sites(trb_seq.upper())
    if sites_after:
        print(f"    Warning: Still has restriction sites: {sites_after}")
    if not trb_mutations:
        print("    No mutations needed")
    print_validation_info("After restriction site removal", trb_seq)
    
    # Step 5: Add TRBC
    print("\n  Step 5: Adding TRBC")
    trb_seq, c_overlap = align_and_merge_with_c(trb_seq, trbc_seq)
    print(f"    TRBC attached (no overlap detection)")
    print_validation_info("Final TRB sequence", trb_seq)
    
    # Save validation file
    if output_dir:
        validation_content = f"""TRB Assembly Validation ({'Seq_CDRs' if use_seq_cdrs else 'IMGT_CDRs'})

Templates:
- TRBV ({trbv_name}): {len(trbv_seq)} bp
  {trbv_seq}
  Translation: {translate(trbv_seq)}

- TRBJ ({trbj_name}): {len(trbj_seq)} bp
  {trbj_seq}
  Translation: {translate(trbj_seq)}

- TRBC ({trbc_name}): {len(trbc_seq)} bp
  {trbc_seq}
  Translation: {translate(trbc_seq)}

- TRB CDR3: {len(cdr3_nt)} bp
  {cdr3_nt}
  Translation: {cdr3_aa}

Assembly Steps:
1. TRBV + CDR3 (overlap: {overlap_len} bp)
   {trb_seq[:len(trbv_seq) + len(cdr3_nt)]}
   Translation: {translate(trb_seq[:len(trbv_seq) + len(cdr3_nt)])}

2. + TRBJ (overlap: {j_overlap} bp)
   {trb_seq}
   Translation: {translate(trb_seq)}

3. + TRBC (overlap: {c_overlap} bp)
   Final: {trb_seq}
   Translation: {translate(trb_seq)}
"""
        save_validation_file(validation_content, os.path.join(output_dir, f"TRB_{'Seq' if use_seq_cdrs else 'IMGT'}_CDRs_validation.txt"))
    
    return trb_seq, trb_mutations


def assemble_tra_sequence(row, use_seq_cdrs=False, tcr_id="", output_dir=""):
    """Assemble TRA sequence with detailed validation output."""
    print(f"\n{'='*60}")
    print(f"Assembling TRA sequence ({'Seq_CDRs' if use_seq_cdrs else 'IMGT_CDRs'})")
    print(f"{'='*60}")
    
    # Get gene names
    trav_name = row['TRA_v_gene']
    traj_name = row['TRA_j_gene']
    cdr3_nt = row['TRA_cdr3_nt']
    cdr3_aa = row['TRA_cdr3']
    
    # Load templates
    trav_path = get_template_path(trav_name, "TRAV")
    traj_path = get_template_path(traj_name, "TRAJ")
    trac_path = get_template_path("TRAC", "TRAC")
    
    trav_seq = read_fasta(trav_path)
    traj_seq = read_fasta(traj_path)
    trac_seq = read_fasta(trac_path)
    
    if not all([trav_seq, traj_seq, trac_seq, cdr3_nt]):
        print("  Error: Missing templates or CDR3")
        return None, []
    
    # Print template sequences
    print("\n  Template Sequences:")
    print_validation_info(f"TRAV ({trav_name})", trav_seq)
    print_validation_info(f"TRAJ ({traj_name})", traj_seq)
    print_validation_info("TRAC", trac_seq)
    print_validation_info("TRA CDR3", cdr3_nt, cdr3_aa)
    
    # Step 1: Align TRAV with CDR3
    print("\n  Step 1: Aligning TRAV with CDR3")
    tra_seq, overlap_len, pos = align_and_merge_v_cdr3(trav_seq, cdr3_nt)
    print(f"    Overlap: {overlap_len} bp at position {pos}")
    print_validation_info("TRAV + CDR3", tra_seq)
    
    # Step 2: Align CDR3 with TRAJ
    print("\n  Step 2: Aligning CDR3 with TRAJ")
    tra_seq, j_overlap = align_and_merge_cdr3_j(tra_seq, traj_seq, cdr3_nt)
    print(f"    J overlap: {j_overlap} bp")
    print_validation_info("TRAV + CDR3 + TRAJ", tra_seq)
    
    # Step 3: Handle CDR1 and CDR2 if available
    cdr1_nn = row.get('TRA_CDR1_nn', '')
    cdr1_aa = row.get('TRA_CDR1_aa', '')
    cdr2_nn = row.get('TRA_CDR2_nn', '')
    cdr2_aa = row.get('TRA_CDR2_aa', '')
    
    if use_seq_cdrs and cdr1_nn and cdr2_nn and pd.notna(cdr1_nn) and pd.notna(cdr2_nn):
        print("\n  Step 3: Replacing CDR1 and CDR2")
        if cdr1_nn and cdr1_aa:
            print_validation_info("TRA CDR1 (Seq)", cdr1_nn, str(cdr1_aa))
            tra_seq = replace_cdr_in_sequence(tra_seq, cdr1_nn, str(cdr1_aa), "CDR1")
        if cdr2_nn and cdr2_aa:
            print_validation_info("TRA CDR2 (Seq)", cdr2_nn, str(cdr2_aa))
            tra_seq = replace_cdr_in_sequence(tra_seq, cdr2_nn, str(cdr2_aa), "CDR2")
        print_validation_info("After CDR replacement", tra_seq)
    
    # Step 4: Remove restriction sites
    print("\n  Step 4: Removing restriction sites")
    sites_before = check_restriction_sites(tra_seq.upper())
    if sites_before:
        print(f"    Found restriction sites: {sites_before}")
    tra_seq, tra_mutations = remove_restriction_sites(tra_seq)
    tra_mutations = tra_mutations or []
    sites_after = check_restriction_sites(tra_seq.upper())
    if sites_after:
        print(f"    Warning: Still has restriction sites: {sites_after}")
    if not tra_mutations:
        print("    No mutations needed")
    print_validation_info("After restriction site removal", tra_seq)
    
    # Step 7: Handle TRAC attachment
    print("\n  Step 5: Adding TRAC")
    ends_with_at = tra_seq.lower().endswith("at")
    print(f"    TRA ends with 'at': {ends_with_at}")
    print(f"    Keeping 'at' ending and attaching TRAC directly")
    
    # Initialize variables for validation file
    best_removal = None
    
    # Do NOT remove "at" ending - keep it and attach TRAC directly
    # Test different TRAC_head removal scenarios (0, 1, 2 positions)
    best_tra_seq = None
    
    for removal in [0, 1, 2]:
        trac_to_test = trac_seq[removal:]
        test_tra_seq = tra_seq + trac_to_test.upper()
        test_aa = translate(test_tra_seq)
        
        print(f"\n    Testing removal of {removal} positions from TRAC_head:")
        print(f"      TRAC_head (removed {removal}): {trac_to_test}")
        print(f"      Full TRA sequence length: {len(test_tra_seq)} bp")
        print(f"      Translation length: {len(test_aa)} aa")
        print(f"      Translation ends with: {test_aa[-5:] if len(test_aa) >= 5 else test_aa}")
        
        # Check if translation ends with "IQNPD"
        if test_aa.endswith("IQNPD"):
            print(f"      ✓ Translation ends with IQNPD - SELECTED")
            best_tra_seq = test_tra_seq
            best_removal = removal
            break
        else:
            print(f"      ✗ Translation does not end with IQNPD")
    
    if best_tra_seq:
        tra_seq = best_tra_seq
        print(f"\n    Selected: Removal of {best_removal} positions from TRAC_head")
    else:
        # If none match, use the original logic (remove 2 positions)
        print(f"\n    Warning: No removal scenario produced IQNPD ending, using default (remove 2 positions)")
        best_removal = 2
        tra_seq = tra_seq + trac_seq[2:].upper()
    
    print_validation_info("Final TRA sequence", tra_seq)
    
    # Save validation file
    if output_dir:
        validation_content = f"""TRA Assembly Validation ({'Seq_CDRs' if use_seq_cdrs else 'IMGT_CDRs'})

Templates:
- TRAV ({trav_name}): {len(trav_seq)} bp
  {trav_seq}
  Translation: {translate(trav_seq)}

- TRAJ ({traj_name}): {len(traj_seq)} bp
  {traj_seq}
  Translation: {translate(traj_seq)}

- TRAC: {len(trac_seq)} bp
  {trac_seq}
  Translation: {translate(trac_seq)}

- TRA CDR3: {len(cdr3_nt)} bp
  {cdr3_nt}
  Translation: {cdr3_aa}

Assembly Steps:
1. TRAV + CDR3 (overlap: {overlap_len} bp)
   {tra_seq[:len(trav_seq) + len(cdr3_nt)]}
   Translation: {translate(tra_seq[:len(trav_seq) + len(cdr3_nt)])}

2. + TRAJ (overlap: {j_overlap} bp)
   {tra_seq_before_trac if ends_with_at else (tra_seq[:len(tra_seq) - len(trac_seq)] if not ends_with_at else '')}
   Translation: {translate(tra_seq_before_trac if ends_with_at else (tra_seq[:len(tra_seq) - len(trac_seq)] if not ends_with_at else ''))}

3. + TRAC (ends with 'at': {ends_with_at}, removal: {best_removal if ends_with_at else 'N/A'})
   Final: {tra_seq}
   Translation: {translate(tra_seq)}
"""
        save_validation_file(validation_content, os.path.join(output_dir, f"TRA_{'Seq' if use_seq_cdrs else 'IMGT'}_CDRs_validation.txt"))
    
    return tra_seq, tra_mutations


def format_trb_tra_sequence(seq, cdr3_nt, mutations, v_seq_length=None):
    """
    Format TRB or TRA sequence with proper casing:
    - Uppercase: mutated restriction sites, CDR3
    - Lowercase: everything else
    
    Args:
        seq: The sequence to format
        cdr3_nt: The CDR3 nucleotide sequence
        mutations: List of mutations (site_name, site_pos, site_seq, mut_pos, old_codon, new_codon)
        v_seq_length: Length of V gene sequence (to help locate CDR3)
    
    Returns: Formatted sequence
    """
    if not seq:
        return seq
    
    seq_list = list(seq.lower())  # Start with all lowercase
    
    # Mark CDR3 as uppercase
    if cdr3_nt:
        cdr3_nt_upper = cdr3_nt.upper()
        seq_upper = ''.join(seq_list).upper()
        cdr3_pos = seq_upper.find(cdr3_nt_upper)
        if cdr3_pos != -1:
            for i in range(len(cdr3_nt)):
                if cdr3_pos + i < len(seq_list):
                    seq_list[cdr3_pos + i] = cdr3_nt_upper[i]
    
    # Mark mutated codons as uppercase
    if mutations:
        for site_name, site_pos, site_seq, mut_pos, old_codon, new_codon in mutations:
            # Mark the mutated codon in uppercase
            for i in range(len(new_codon)):
                if mut_pos + i < len(seq_list):
                    seq_list[mut_pos + i] = new_codon[i].upper()
    
    return ''.join(seq_list)


def format_final_assembly(seq, trb_cdr3_nt, tra_cdr3_nt, trb_mutations=None, tra_mutations=None, component_info=None):
    """
    Format final TRB-P2A-TRA assembly with proper casing:
    - Uppercase: ASCI_SITE, KOZAK_SEQ, P2A, PSTI_SITE, mutated restriction sites, CDR3s
    - Lowercase: everything else
    
    Args:
        seq: The full assembled sequence (all lowercase)
        trb_cdr3_nt: TRB CDR3 nucleotide sequence
        tra_cdr3_nt: TRA CDR3 nucleotide sequence
        trb_mutations: List of TRB mutations
        tra_mutations: List of TRA mutations
        component_info: Dict with 'trb_length', 'p2a_length', 'tra_length'
    
    Returns: Formatted sequence
    """
    if not seq:
        return seq
    
    seq_list = list(seq.lower())  # Start with all lowercase
    
    # We know the structure: ASCI_SITE + KOZAK_SEQ + TRB + P2A + TRA + PSTI_SITE
    asci_pos = 0
    kozak_pos = len(ASCI_SITE)
    
    # Load P2A sequence
    p2a_path = os.path.join(TEMPLATES_DIR, "codon_optimized_P2A.fasta")
    p2a_seq = read_fasta(p2a_path)
    
    # Calculate positions based on known structure
    trb_start = kozak_pos + len(KOZAK_SEQ)
    if component_info:
        trb_length = component_info['trb_length']
        p2a_length = component_info['p2a_length']
        tra_length = component_info['tra_length']
    else:
        # Estimate if not provided
        trb_length = 950  # Typical TRB length
        p2a_length = len(p2a_seq) if p2a_seq else 60
        tra_length = 400  # Typical TRA length
    
    p2a_pos = trb_start + trb_length
    tra_start = p2a_pos + p2a_length
    psti_pos = len(seq) - len(PSTI_SITE)
    
    # Mark ASCI_SITE as uppercase
    for i in range(len(ASCI_SITE)):
        if asci_pos + i < len(seq_list):
            seq_list[asci_pos + i] = ASCI_SITE[i]
    
    # Mark KOZAK_SEQ as uppercase
    for i in range(len(KOZAK_SEQ)):
        if kozak_pos + i < len(seq_list):
            seq_list[kozak_pos + i] = KOZAK_SEQ[i]
    
    # Mark P2A as uppercase
    if p2a_seq and p2a_pos < len(seq_list):
        for i in range(len(p2a_seq)):
            if p2a_pos + i < len(seq_list):
                seq_list[p2a_pos + i] = p2a_seq[i].upper()
    
    # Mark PSTI_SITE as uppercase
    if psti_pos >= 0:
        for i in range(len(PSTI_SITE)):
            if psti_pos + i < len(seq_list):
                seq_list[psti_pos + i] = PSTI_SITE[i]
    
    # Mark TRB CDR3 as uppercase
    if trb_cdr3_nt:
        trb_cdr3_upper = trb_cdr3_nt.upper()
        seq_upper = ''.join(seq_list).upper()
        trb_cdr3_pos = seq_upper.find(trb_cdr3_upper, trb_start)
        if trb_cdr3_pos != -1:
            for i in range(len(trb_cdr3_nt)):
                if trb_cdr3_pos + i < len(seq_list):
                    seq_list[trb_cdr3_pos + i] = trb_cdr3_upper[i]
    
    # Mark TRA CDR3 as uppercase
    if tra_cdr3_nt:
        tra_cdr3_upper = tra_cdr3_nt.upper()
        seq_upper = ''.join(seq_list).upper()
        tra_cdr3_pos = seq_upper.find(tra_cdr3_upper, tra_start)
        if tra_cdr3_pos != -1:
            for i in range(len(tra_cdr3_nt)):
                if tra_cdr3_pos + i < len(seq_list):
                    seq_list[tra_cdr3_pos + i] = tra_cdr3_upper[i]
    
    # Mark TRB mutated codons as uppercase
    if trb_mutations:
        for site_name, site_pos, site_seq, mut_pos, old_codon, new_codon in trb_mutations:
            # Adjust position for ASCI + KOZAK in final assembly
            adjusted_mut_pos = trb_start + mut_pos
            for i in range(len(new_codon)):
                if adjusted_mut_pos + i < len(seq_list):
                    seq_list[adjusted_mut_pos + i] = new_codon[i].upper()
    
    # Mark TRA mutated codons as uppercase
    if tra_mutations:
        for site_name, site_pos, site_seq, mut_pos, old_codon, new_codon in tra_mutations:
            adjusted_mut_pos = tra_start + mut_pos
            for i in range(len(new_codon)):
                if adjusted_mut_pos + i < len(seq_list):
                    seq_list[adjusted_mut_pos + i] = new_codon[i].upper()
    
    return ''.join(seq_list)


def assemble_trb_p2a_tra(trb_seq, tra_seq, trb_mutations=None, tra_mutations=None, trb_cdr3_nt=None, tra_cdr3_nt=None):
    """
    Assemble final TRB-P2A-TRA sequence.
    Returns unformatted sequence (formatting will be applied before saving).
    Also returns component lengths for formatting.
    """
    if not trb_seq or not tra_seq:
        return None, None
    
    # Load P2A sequence
    p2a_path = os.path.join(TEMPLATES_DIR, "codon_optimized_P2A.fasta")
    p2a_seq = read_fasta(p2a_path)
    if not p2a_seq:
        return None, None
    
    # Assemble: AscI + Kozak + TRB + P2A + TRA + PstI (all lowercase for now, will format later)
    final_seq = ASCI_SITE.lower() + KOZAK_SEQ.lower() + trb_seq.lower() + p2a_seq.lower() + tra_seq.lower() + PSTI_SITE.lower()
    
    # Return sequence and component info for formatting
    component_info = {
        'trb_length': len(trb_seq),
        'p2a_length': len(p2a_seq),
        'tra_length': len(tra_seq)
    }
    
    return final_seq, component_info


def save_sequence(seq, filepath):
    """Save sequence to FASTA file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        f.write(f">{os.path.basename(filepath)}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")


def process_tcr(row):
    """Process a single TCR entry with detailed validation."""
    tcr_id = row['tcr_id']
    print(f"\n{'#'*60}")
    print(f"Processing TCR: {tcr_id}")
    print(f"{'#'*60}")
    
    # Create output directory for this TCR
    tcr_output_dir = os.path.join(OUTPUT_DIR, tcr_id)
    os.makedirs(tcr_output_dir, exist_ok=True)
    
    # Assemble TRB sequences
    trb_imgt_result = assemble_trb_sequence(row, use_seq_cdrs=False, tcr_id=tcr_id, output_dir=tcr_output_dir)
    trb_seq_result = assemble_trb_sequence(row, use_seq_cdrs=True, tcr_id=tcr_id, output_dir=tcr_output_dir)
    
    if trb_imgt_result:
        trb_imgt, trb_imgt_mutations = trb_imgt_result
    else:
        trb_imgt, trb_imgt_mutations = None, []
    
    if trb_seq_result:
        trb_seq, trb_seq_mutations = trb_seq_result
    else:
        trb_seq, trb_seq_mutations = None, []
    
    # Assemble TRA sequences
    tra_imgt_result = assemble_tra_sequence(row, use_seq_cdrs=False, tcr_id=tcr_id, output_dir=tcr_output_dir)
    tra_seq_result = assemble_tra_sequence(row, use_seq_cdrs=True, tcr_id=tcr_id, output_dir=tcr_output_dir)
    
    if tra_imgt_result:
        tra_imgt, tra_imgt_mutations = tra_imgt_result
    else:
        tra_imgt, tra_imgt_mutations = None, []
    
    if tra_seq_result:
        tra_seq, tra_seq_mutations = tra_seq_result
    else:
        tra_seq, tra_seq_mutations = None, []
    
    # Get CDR3 sequences for formatting
    trb_cdr3_nt = row.get('TRB_cdr3_nt', '')
    tra_cdr3_nt = row.get('TRA_cdr3_nt', '')
    
    # Save TRB sequences (formatted)
    if trb_imgt:
        formatted_trb_imgt = format_trb_tra_sequence(trb_imgt, trb_cdr3_nt, trb_imgt_mutations)
        save_sequence(formatted_trb_imgt, os.path.join(tcr_output_dir, "TRB_IMGT_CDRs.fasta"))
    
    if trb_seq:
        formatted_trb_seq = format_trb_tra_sequence(trb_seq, trb_cdr3_nt, trb_seq_mutations)
        save_sequence(formatted_trb_seq, os.path.join(tcr_output_dir, "TRB_Seq_CDRs.fasta"))
    
    # Save TRA sequences (formatted)
    if tra_imgt:
        formatted_tra_imgt = format_trb_tra_sequence(tra_imgt, tra_cdr3_nt, tra_imgt_mutations)
        save_sequence(formatted_tra_imgt, os.path.join(tcr_output_dir, "TRA_IMGT_CDRs.fasta"))
    
    if tra_seq:
        formatted_tra_seq = format_trb_tra_sequence(tra_seq, tra_cdr3_nt, tra_seq_mutations)
        save_sequence(formatted_tra_seq, os.path.join(tcr_output_dir, "TRA_Seq_CDRs.fasta"))
    
    # Assemble final TRB-P2A-TRA sequences
    print(f"\n{'='*60}")
    print("Assembling Final TRB-P2A-TRA Sequences")
    print(f"{'='*60}")
    
    # TRB_IMGT_CDRs with TRA_IMGT_CDRs
    if trb_imgt and tra_imgt:
        print("\n  TRB_P2A_TRA_IMGT:")
        final_imgt, comp_info_imgt = assemble_trb_p2a_tra(trb_imgt, tra_imgt, trb_imgt_mutations, tra_imgt_mutations, trb_cdr3_nt, tra_cdr3_nt)
        if final_imgt:
            formatted_final_imgt = format_final_assembly(final_imgt, trb_cdr3_nt, tra_cdr3_nt, trb_imgt_mutations, tra_imgt_mutations, comp_info_imgt)
            print_validation_info("Final TRB-P2A-TRA (IMGT)", formatted_final_imgt)
            save_sequence(formatted_final_imgt, os.path.join(tcr_output_dir, "TRB_P2A_TRA_IMGT.fasta"))
    
    # TRB_Seq_CDRs with TRA_Seq_CDRs
    if trb_seq and tra_seq:
        print("\n  TRB_P2A_TRA_Seq:")
        final_seq, comp_info_seq = assemble_trb_p2a_tra(trb_seq, tra_seq, trb_seq_mutations, tra_seq_mutations, trb_cdr3_nt, tra_cdr3_nt)
        if final_seq:
            formatted_final_seq = format_final_assembly(final_seq, trb_cdr3_nt, tra_cdr3_nt, trb_seq_mutations, tra_seq_mutations, comp_info_seq)
            print_validation_info("Final TRB-P2A-TRA (Seq)", formatted_final_seq)
            save_sequence(formatted_final_seq, os.path.join(tcr_output_dir, "TRB_P2A_TRA_Seq.fasta"))
    
    return True


class Tee:
    """Class to write to both file and stdout."""
    def __init__(self, file_path):
        self.file = open(file_path, 'w', encoding='utf-8')
        self.stdout = sys.stdout
    
    def write(self, text):
        self.file.write(text)
        self.stdout.write(text)
        self.file.flush()
    
    def flush(self):
        self.file.flush()
        self.stdout.flush()
    
    def close(self):
        self.file.close()


def main(input_file=None, output_dir=None, max_tcrs=None):
    """
    Main function.
    
    Args:
        input_file: Path to input CSV file (default: uses global INPUT_FILE)
        output_dir: Path to output directory (default: uses global OUTPUT_DIR)
        max_tcrs: Maximum number of TCRs to process (default: uses global MAX_TCRS)
    
    Returns:
        tuple: (successful, failed) counts
    """
    # Declare global variables first
    global OUTPUT_DIR
    
    # Use provided parameters or fall back to global variables
    input_path = input_file if input_file is not None else INPUT_FILE
    output_path = output_dir if output_dir is not None else OUTPUT_DIR
    max_tcrs_limit = max_tcrs if max_tcrs is not None else MAX_TCRS
    
    # Update global OUTPUT_DIR for use in process_tcr function
    OUTPUT_DIR = output_path
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Set up log file
    log_file = os.path.join(OUTPUT_DIR, f"test_assembly_log_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt")
    tee = Tee(log_file)
    sys.stdout = tee
    
    try:
        print("=" * 60)
        if max_tcrs_limit is not None:
            print(f"TCR Sequence Assembly - TEST MODE (First {max_tcrs_limit} TCRs)")
        else:
            print("TCR Sequence Assembly - TEST MODE (All TCRs)")
        print("=" * 60)
        print(f"Log file: {log_file}\n")
        
        # Load input file
        print(f"\nLoading input file: {input_path}")
        df = pd.read_csv(input_path)
        print(f"Loaded {len(df)} TCR entries")
        
        # Limit to MAX_TCRS if specified
        if max_tcrs_limit is not None and max_tcrs_limit > 0:
            df = df.head(max_tcrs_limit)
            print(f"\nProcessing first {len(df)} TCR entries (MAX_TCRS = {max_tcrs_limit})...")
        else:
            print(f"\nProcessing all {len(df)} TCR entries...")
        
        # Process each TCR
        successful = 0
        failed = 0
        
        for idx, row in df.iterrows():
            try:
                if process_tcr(row):
                    successful += 1
                else:
                    failed += 1
            except Exception as e:
                print(f"Error processing {row.get('tcr_id', 'unknown')}: {e}")
                import traceback
                traceback.print_exc()
                failed += 1
        
        print(f"\n{'='*60}")
        print("Test complete!")
        print(f"Summary: {successful} successful, {failed} failed")
        print(f"Output directory: {OUTPUT_DIR}")
        print(f"Log file: {log_file}")
        print("=" * 60)
        
        return successful, failed
    finally:
        # Restore stdout and close log file
        sys.stdout = tee.stdout
        tee.close()
        print(f"\nLog saved to: {log_file}")


if __name__ == "__main__":
    main()

