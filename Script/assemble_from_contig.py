#!/usr/bin/env python3
"""
Script to assemble TRB-P2A-TRA sequences from full contig sequences.

This script:
1. Aligns full contig sequences to V gene templates to find assembly start
2. Uses full contig sequences directly (no V-CDR3-J assembly)
3. Removes restriction sites by silent mutations
4. Completes TRBC by finding J-C junction and replacing with template
5. Keeps TRAJ-TRAC junction from contig, ensures TRAC length matches template
6. Assembles final TRB-P2A-TRA sequences
"""

import os
import re
import pandas as pd
from pathlib import Path

# Input and output paths
INPUT_FILE = "Input/hERV/hERV_TCR_synthesize_input_010226_with_CDR12.csv"  # Using 010226 for testing; change to 010216 when available
TEMPLATES_DIR = "Templates"
OUTPUT_DIR = "Output_original_contig"

# Number of TCRs to process (None = process all)
MAX_TCRS = None  # Set to a number (e.g., 10) to process only the first N TCRs

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


def check_stop_codons_in_sequence(aa_seq, allow_terminal_stop=True):
    """
    Check for stop codons in amino acid sequence.
    
    Args:
        aa_seq: Amino acid sequence
        allow_terminal_stop: If True, allows stop codon at the very end
    
    Returns:
        (has_stop_codons, stop_positions) where stop_positions is a list of positions
    """
    stop_positions = []
    for i, aa in enumerate(aa_seq):
        if aa == '*':
            # If allow_terminal_stop and this is the last position, skip
            if allow_terminal_stop and i == len(aa_seq) - 1:
                continue
            stop_positions.append(i)
    
    return len(stop_positions) > 0, stop_positions


def validate_and_log_translation(trb_seq, tra_seq, tcr_id, output_dir, trb_mutations=None, tra_mutations=None, trb_failed_sites=None, tra_failed_sites=None):
    """
    Validate assembled sequences by translating to amino acids and checking for stop codons.
    Also logs restriction site removal information.
    Writes validation log to a text file.
    
    Args:
        trb_seq: TRB nucleotide sequence
        tra_seq: TRA nucleotide sequence
        tcr_id: TCR identifier
        output_dir: Output directory for this TCR
        trb_mutations: List of TRB mutations (site_name, site_pos, site_seq, mut_pos, old_codon, new_codon)
        tra_mutations: List of TRA mutations (site_name, site_pos, site_seq, mut_pos, old_codon, new_codon)
    
    Returns:
        (trb_valid, tra_valid) - Boolean tuples indicating if sequences are valid
    """
    log_file = os.path.join(output_dir, f"{tcr_id}_validation_log.txt")
    trb_mutations = trb_mutations or []
    tra_mutations = tra_mutations or []
    trb_failed_sites = trb_failed_sites or []
    tra_failed_sites = tra_failed_sites or []
    
    with open(log_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write(f"Translation Validation Log for {tcr_id}\n")
        f.write("=" * 80 + "\n\n")
        
        # Check for remaining restriction sites
        trb_remaining_sites = check_restriction_sites(trb_seq)
        tra_remaining_sites = check_restriction_sites(tra_seq)
        
        # Validate TRB sequence
        f.write("TRB Sequence Validation\n")
        f.write("-" * 80 + "\n")
        f.write(f"Nucleotide length: {len(trb_seq)}\n")
        f.write(f"In frame: {len(trb_seq) % 3 == 0}\n\n")
        
        # Log restriction site information for TRB
        f.write("TRB Restriction Site Removal\n")
        f.write("-" * 80 + "\n")
        if trb_mutations:
            f.write(f"Total restriction sites found and removed: {len(trb_mutations)}\n\n")
            for site_name, site_pos, site_seq, mut_pos, old_codon, new_codon in trb_mutations:
                old_aa = GENETIC_CODE.get(old_codon, '?')
                new_aa = GENETIC_CODE.get(new_codon, '?')
                f.write(f"  ✓ Removed {site_name} site at position {site_pos}: {site_seq}\n")
                f.write(f"    Mutation at position {mut_pos}: {old_codon} -> {new_codon}\n")
                f.write(f"    Amino acid preserved: {old_aa} = {new_aa}\n\n")
        else:
            f.write("No restriction sites found in TRB sequence.\n\n")
        
        # Check for remaining restriction sites in TRB
        if trb_failed_sites:
            f.write("WARNING: Restriction sites in TRB sequence that could not be removed:\n")
            for site_name, pos, site_seq in trb_failed_sites:
                f.write(f"  ✗ {site_name} site at position {pos}: {site_seq}\n")
            f.write("\n")
        elif trb_remaining_sites:
            f.write("WARNING: Remaining restriction sites in TRB sequence (could not be removed):\n")
            for site_name, positions in trb_remaining_sites.items():
                site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
                if isinstance(positions, list):
                    for pos in positions:
                        actual_seq = trb_seq[pos:pos+len(site_seq)]
                        f.write(f"  ✗ {site_name} site at position {pos}: {actual_seq}\n")
                else:
                    pos = positions
                    actual_seq = trb_seq[pos:pos+len(site_seq)]
                    f.write(f"  ✗ {site_name} site at position {pos}: {actual_seq}\n")
            f.write("\n")
        else:
            f.write("✓ No remaining restriction sites in TRB sequence.\n\n")
        
        f.write("\n")
        
        # Always translate and check for stop codons, even if not in frame
        trb_aa = translate(trb_seq)
        f.write(f"Amino acid length: {len(trb_aa)}\n")
        f.write(f"Amino acid sequence:\n{trb_aa}\n\n")
        
        if len(trb_seq) % 3 != 0:
            f.write("WARNING: TRB sequence is not in frame (length not divisible by 3)\n")
            f.write("(Translation may be incomplete or shifted)\n\n")
        
        has_stops, stop_positions = check_stop_codons_in_sequence(trb_aa, allow_terminal_stop=True)
        if has_stops:
            f.write(f"WARNING: TRB sequence contains stop codons at positions: {stop_positions}\n")
            f.write("(Stop codons in the middle of the sequence indicate potential issues)\n\n")
            trb_valid = False
        else:
            f.write("✓ TRB sequence: No stop codons found (valid)\n\n")
            trb_valid = True
        
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Validate TRA sequence
        f.write("TRA Sequence Validation\n")
        f.write("-" * 80 + "\n")
        f.write(f"Nucleotide length: {len(tra_seq)}\n")
        f.write(f"In frame: {len(tra_seq) % 3 == 0}\n\n")
        
        # Log restriction site information for TRA
        f.write("TRA Restriction Site Removal\n")
        f.write("-" * 80 + "\n")
        if tra_mutations:
            f.write(f"Total restriction sites found and removed: {len(tra_mutations)}\n\n")
            for site_name, site_pos, site_seq, mut_pos, old_codon, new_codon in tra_mutations:
                old_aa = GENETIC_CODE.get(old_codon, '?')
                new_aa = GENETIC_CODE.get(new_codon, '?')
                f.write(f"  ✓ Removed {site_name} site at position {site_pos}: {site_seq}\n")
                f.write(f"    Mutation at position {mut_pos}: {old_codon} -> {new_codon}\n")
                f.write(f"    Amino acid preserved: {old_aa} = {new_aa}\n\n")
        else:
            f.write("No restriction sites found in TRA sequence.\n\n")
        
        # Check for remaining restriction sites in TRA
        if tra_failed_sites:
            f.write("WARNING: Restriction sites in TRA sequence that could not be removed:\n")
            for site_name, pos, site_seq in tra_failed_sites:
                f.write(f"  ✗ {site_name} site at position {pos}: {site_seq}\n")
            f.write("\n")
        elif tra_remaining_sites:
            f.write("WARNING: Remaining restriction sites in TRA sequence (could not be removed):\n")
            for site_name, positions in tra_remaining_sites.items():
                site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
                if isinstance(positions, list):
                    for pos in positions:
                        actual_seq = tra_seq[pos:pos+len(site_seq)]
                        f.write(f"  ✗ {site_name} site at position {pos}: {actual_seq}\n")
                else:
                    pos = positions
                    actual_seq = tra_seq[pos:pos+len(site_seq)]
                    f.write(f"  ✗ {site_name} site at position {pos}: {actual_seq}\n")
            f.write("\n")
        else:
            f.write("✓ No remaining restriction sites in TRA sequence.\n\n")
        
        f.write("\n")
        
        # Always translate and check for stop codons, even if not in frame
        tra_aa = translate(tra_seq)
        f.write(f"Amino acid length: {len(tra_aa)}\n")
        f.write(f"Amino acid sequence:\n{tra_aa}\n\n")
        
        if len(tra_seq) % 3 != 0:
            f.write("NOTE: TRA sequence is not in frame (length not divisible by 3)\n")
            f.write("(Translation may be incomplete or shifted, but this is OK for contig sequences)\n\n")
        
        has_stops, stop_positions = check_stop_codons_in_sequence(tra_aa, allow_terminal_stop=True)
        if has_stops:
            f.write(f"WARNING: TRA sequence contains stop codons at positions: {stop_positions}\n")
            f.write("(Stop codons in the middle of the sequence indicate potential issues)\n\n")
            tra_valid = False
        else:
            f.write("✓ TRA sequence: No stop codons found (valid)\n\n")
            tra_valid = True
        
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Summary
        f.write("Validation Summary\n")
        f.write("-" * 80 + "\n")
        f.write(f"TRB valid: {'Yes' if trb_valid else 'No'}\n")
        f.write(f"TRA valid: {'Yes' if tra_valid else 'No'}\n")
        f.write(f"TRB restriction sites removed: {len(trb_mutations)}\n")
        f.write(f"TRA restriction sites removed: {len(tra_mutations)}\n")
        f.write(f"TRB restriction sites failed to remove: {len(trb_failed_sites)}\n")
        f.write(f"TRA restriction sites failed to remove: {len(tra_failed_sites)}\n")
        # Count remaining sites (handle both list and single value formats)
        trb_remaining_count = sum(len(v) if isinstance(v, list) else 1 for v in trb_remaining_sites.values()) if trb_remaining_sites else 0
        tra_remaining_count = sum(len(v) if isinstance(v, list) else 1 for v in tra_remaining_sites.values()) if tra_remaining_sites else 0
        f.write(f"TRB remaining restriction sites: {trb_remaining_count}\n")
        f.write(f"TRA remaining restriction sites: {tra_remaining_count}\n")
        f.write(f"Overall: {'VALID' if (trb_valid and tra_valid and len(trb_remaining_sites) == 0 and len(tra_remaining_sites) == 0) else 'INVALID - Contains issues'}\n")
        f.write("=" * 80 + "\n")
    
    return trb_valid, tra_valid


def find_v_gene_start_in_contig(contig_seq, v_template):
    """
    Find where the V gene template starts in the contig sequence.
    Uses alignment to find the best match position.
    
    Returns: (start_position, match_score) or (None, 0) if not found
    """
    contig_seq = contig_seq.upper()
    v_template = v_template.upper()
    
    # Try to find the V gene start by looking for a significant match
    # We'll look for a substring match of at least 30 nucleotides from the V template start
    v_start_region = v_template[:min(60, len(v_template))]  # First 60 nt of V
    
    best_match_pos = None
    best_match_len = 0
    
    # Try different positions in contig
    for i in range(len(contig_seq) - 30):
        # Try matching different lengths
        for match_len in range(30, min(len(v_start_region), len(contig_seq) - i) + 1):
            contig_region = contig_seq[i:i+match_len]
            v_region = v_start_region[:match_len]
            
            # Count matches
            matches = sum(1 for a, b in zip(contig_region, v_region) if a == b)
            match_ratio = matches / match_len
            
            # If we have a good match (>= 80% identity), use it
            if match_ratio >= 0.80 and match_len > best_match_len:
                best_match_pos = i
                best_match_len = match_len
    
    if best_match_pos is not None:
        return best_match_pos, best_match_len / len(v_start_region)
    
    # Fallback: try exact substring match
    for i in range(len(contig_seq) - len(v_start_region) + 1):
        if contig_seq[i:i+len(v_start_region)] == v_start_region:
            return i, 1.0
    
    return None, 0.0


def find_j_c_junction(contig_seq, j_template, c_template):
    """
    Find the J-C junction in the contig sequence.
    Looks for where J ends and C begins.
    
    Returns: (junction_position, j_end_pos, c_start_pos) or (None, None, None)
    """
    contig_seq = contig_seq.upper()
    j_template = j_template.upper()
    c_template = c_template.upper()
    
    # Find J gene end in contig (look for end of J template)
    j_end_region = j_template[-min(30, len(j_template)):]  # Last 30 nt of J
    
    j_end_pos = None
    best_j_match = 0
    
    # Find where J ends in contig
    for i in range(len(contig_seq) - len(j_end_region) + 1):
        contig_region = contig_seq[i:i+len(j_end_region)]
        matches = sum(1 for a, b in zip(contig_region, j_end_region) if a == b)
        match_ratio = matches / len(j_end_region)
        
        if match_ratio >= 0.80 and match_ratio > best_j_match:
            j_end_pos = i + len(j_end_region)
            best_j_match = match_ratio
    
    if j_end_pos is None:
        return None, None, None
    
    # Find where C starts in contig (look for start of C template)
    c_start_region = c_template[:min(30, len(c_template))]  # First 30 nt of C
    
    c_start_pos = None
    best_c_match = 0
    
    # Search after J end position
    for i in range(j_end_pos, len(contig_seq) - len(c_start_region) + 1):
        contig_region = contig_seq[i:i+len(c_start_region)]
        matches = sum(1 for a, b in zip(contig_region, c_start_region) if a == b)
        match_ratio = matches / len(c_start_region)
        
        if match_ratio >= 0.80 and match_ratio > best_c_match:
            c_start_pos = i
            best_c_match = match_ratio
    
    if c_start_pos is None:
        return None, j_end_pos, None
    
    # Junction is between j_end_pos and c_start_pos
    return (j_end_pos, c_start_pos), j_end_pos, c_start_pos


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
    Same as in original script.
    """
    site_len = len(site_seq)
    
    # Find all codons that actually contain any part of the restriction site
    overlapping_codons = []
    
    first_codon_start = (site_pos // 3) * 3
    last_codon_start = ((site_pos + site_len - 1) // 3) * 3
    
    for codon_start in range(first_codon_start, last_codon_start + 1, 3):
        codon_end = codon_start + 3
        if codon_end > len(seq):
            break
        
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
            site_gone = True
            # Check a window around the original position (accounting for potential shifts)
            check_start = max(0, site_pos - 3)
            check_end = min(len(test_seq), site_pos + site_len + 3)
            if site_seq in test_seq[check_start:check_end].upper():
                site_gone = False
            
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
                if i >= j:
                    continue
                for alt1 in codon1_alts:
                    for alt2 in codon2_alts:
                        test_seq = seq
                        if codon1_start > codon2_start:
                            test_seq = test_seq[:codon1_start] + alt1 + test_seq[codon1_start + 3:]
                            test_seq = test_seq[:codon2_start] + alt2 + test_seq[codon2_start + 3:]
                        else:
                            test_seq = test_seq[:codon2_start] + alt2 + test_seq[codon2_start + 3:]
                            test_seq = test_seq[:codon1_start] + alt1 + test_seq[codon1_start + 3:]
                        
                        # Check if the specific site at site_pos is gone
                        check_start = max(0, site_pos - 3)
                        check_end = min(len(test_seq), site_pos + site_len + 3)
                        site_gone = site_seq not in test_seq[check_start:check_end].upper()
                        
                        if site_gone:
                            other_site = ASCI_SITE if site_seq == PSTI_SITE else PSTI_SITE
                            other_site_in_region = other_site in test_seq[check_start:check_end].upper()
                            if not other_site_in_region:
                                return test_seq, (codon1_start, codon1_old, alt1)
    
    # Strategy 3: Try mutating all overlapping codons simultaneously
    if len(overlapping_codons) >= 2:
        import itertools
        max_combinations = 100
        combinations_tried = 0
        
        for combo in itertools.product(*[alts for _, _, alts in overlapping_codons]):
            if combinations_tried >= max_combinations:
                break
            combinations_tried += 1
            
            test_seq = seq
            mutations_applied = []
            
            sorted_codons = sorted(overlapping_codons, key=lambda x: x[0], reverse=True)
            for idx, (codon_start, codon_old, _) in enumerate(sorted_codons):
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
                    return test_seq, mutations_applied[-1]
    
    return seq, None


def remove_restriction_sites(seq):
    """
    Remove restriction sites by making silent mutations.
    Returns: (cleaned_sequence, mutations_list, failed_sites_list)
    failed_sites_list contains tuples: (site_name, site_position, site_sequence)
    """
    seq = seq.upper()
    sites = check_restriction_sites(seq)
    mutations = []
    failed_sites = []
    processed_sites = set()  # Track sites we've attempted to remove
    
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
            # Skip if we've already tried this site
            site_key = (site_name, pos)
            if site_key in processed_sites:
                continue
            processed_sites.add(site_key)
            
            site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
            actual_site_seq = seq[pos:pos+len(site_seq)]
            new_seq, mutation_info = make_silent_mutation(seq, pos, site_seq)
            if new_seq != seq and mutation_info:
                mut_pos, old_codon, new_codon = mutation_info
                mutations.append((site_name, pos, actual_site_seq, mut_pos, old_codon, new_codon))
                print(f"    Removed {site_name} site at position {pos}: {actual_site_seq}")
                print(f"      Mutation at position {mut_pos}: {old_codon} -> {new_codon} (silent)")
                seq = new_seq
                sites = check_restriction_sites(seq)
                if not sites:
                    break
            else:
                failed_sites.append((site_name, pos, actual_site_seq))
                print(f"    Warning: Could not remove {site_name} site at position {pos}")
        
        if seq == prev_seq:
            print(f"    Warning: Could not remove all restriction sites after {iteration} iterations")
            # Add any remaining sites to failed_sites
            for site_name, positions in sites.items():
                site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
                if isinstance(positions, list):
                    for pos in positions:
                        site_key = (site_name, pos)
                        if site_key not in processed_sites:
                            actual_site_seq = seq[pos:pos+len(site_seq)]
                            failed_sites.append((site_name, pos, actual_site_seq))
                else:
                    pos = positions
                    site_key = (site_name, pos)
                    if site_key not in processed_sites:
                        actual_site_seq = seq[pos:pos+len(site_seq)]
                        failed_sites.append((site_name, pos, actual_site_seq))
            break
    
    return seq, mutations, failed_sites


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


def process_trb_contig(row):
    """
    Process TRB full contig:
    1. Align to TRBV template to find start
    2. Use full contig from V start
    3. Find J-C junction and replace incomplete TRBC with complete template
    4. Remove restriction sites
    """
    trb_full_contig = row.get('TRB_full_contig', '')
    if not trb_full_contig or pd.isna(trb_full_contig) or trb_full_contig == '':
        print(f"    Warning: TRB_full_contig not available")
        return None, [], []
    
    trb_full_contig = str(trb_full_contig).upper().replace('-', '').replace(' ', '')
    
    # Get gene names and templates
    trbv_name = row['TRB_v_gene']
    trbj_name = row['TRB_j_gene']
    trbc_name = row['TRB_c_gene']
    
    trbv_path = get_template_path(trbv_name, "TRBV")
    trbj_path = get_template_path(trbj_name, "TRBJ")
    trbc_path = get_template_path(trbc_name, "TRBC")
    
    trbv_template = read_fasta(trbv_path)
    trbj_template = read_fasta(trbj_path)
    trbc_template = read_fasta(trbc_path)
    
    if not trbv_template:
        print(f"    Warning: TRBV template not found: {trbv_path}")
        return None, [], []
    if not trbj_template:
        print(f"    Warning: TRBJ template not found: {trbj_path}")
        return None, [], []
    if not trbc_template:
        print(f"    Warning: TRBC template not found: {trbc_path}")
        return None, [], []
    
    # Step 1: Find V gene start in contig
    v_start_pos, match_score = find_v_gene_start_in_contig(trb_full_contig, trbv_template)
    if v_start_pos is None:
        print(f"    Warning: Could not find TRBV start in contig")
        return None, [], []
    
    print(f"    Found TRBV start at position {v_start_pos} (match score: {match_score:.2f})")
    
    # Step 2: Extract contig from V start
    trb_seq = trb_full_contig[v_start_pos:]
    
    # Step 3: Find J-C junction and replace incomplete TRBC with complete template
    junction_info, j_end_pos, c_start_pos = find_j_c_junction(trb_seq, trbj_template, trbc_template)
    
    if junction_info is not None:
        j_end_pos, c_start_pos = junction_info
        print(f"    Found J-C junction: J ends at {j_end_pos}, C starts at {c_start_pos}")
        
        # Replace incomplete C with complete C template
        trb_seq = trb_seq[:c_start_pos] + trbc_template.upper()
        print(f"    Replaced incomplete TRBC with complete template (length: {len(trbc_template)})")
    else:
        print(f"    Warning: Could not find J-C junction, using contig as-is")
        # Try to append TRBC if not found
        if trbc_template:
            trb_seq = trb_seq + trbc_template.upper()
    
    # Step 4: Remove restriction sites
    trb_seq, trb_mutations, trb_failed_sites = remove_restriction_sites(trb_seq)
    trb_mutations = trb_mutations or []
    trb_failed_sites = trb_failed_sites or []
    
    return trb_seq, trb_mutations, trb_failed_sites


def process_tra_contig(row):
    """
    Process TRA full contig:
    1. Align to TRAV template to find start
    2. Use full contig from V start
    3. Keep TRAJ-TRAC junction from contig, ensure TRAC length matches template
    4. Remove restriction sites
    """
    tra_full_contig = row.get('TRA_full_contig', '')
    if not tra_full_contig or pd.isna(tra_full_contig) or tra_full_contig == '':
        print(f"    Warning: TRA_full_contig not available")
        return None, [], []
    
    tra_full_contig = str(tra_full_contig).upper().replace('-', '').replace(' ', '')
    
    # Get gene names and templates
    trav_name = row['TRA_v_gene']
    traj_name = row['TRA_j_gene']
    
    trav_path = get_template_path(trav_name, "TRAV")
    traj_path = get_template_path(traj_name, "TRAJ")
    trac_path = get_template_path("TRAC", "TRAC")
    
    trav_template = read_fasta(trav_path)
    traj_template = read_fasta(traj_path)
    trac_template = read_fasta(trac_path)
    
    if not trav_template:
        print(f"    Warning: TRAV template not found: {trav_path}")
        return None, [], []
    if not traj_template:
        print(f"    Warning: TRAJ template not found: {traj_path}")
        return None, [], []
    if not trac_template:
        print(f"    Warning: TRAC template not found: {trac_path}")
        return None, [], []
    
    # Step 1: Find V gene start in contig
    v_start_pos, match_score = find_v_gene_start_in_contig(tra_full_contig, trav_template)
    if v_start_pos is None:
        print(f"    Warning: Could not find TRAV start in contig")
        return None, [], []
    
    print(f"    Found TRAV start at position {v_start_pos} (match score: {match_score:.2f})")
    
    # Step 2: Extract contig from V start
    tra_seq = tra_full_contig[v_start_pos:]
    
    # Step 3: Find TRAJ end and ensure TRAC length matches template
    # Look for TRAJ end in the contig
    traj_end_region = traj_template[-min(30, len(traj_template)):]  # Last 30 nt of TRAJ
    
    traj_end_pos = None
    best_match = 0
    
    for i in range(len(tra_seq) - len(traj_end_region) + 1):
        contig_region = tra_seq[i:i+len(traj_end_region)]
        matches = sum(1 for a, b in zip(contig_region, traj_end_region) if a == b)
        match_ratio = matches / len(traj_end_region)
        
        if match_ratio >= 0.80 and match_ratio > best_match:
            traj_end_pos = i + len(traj_end_region)
            best_match = match_ratio
    
    if traj_end_pos is not None:
        print(f"    Found TRAJ end at position {traj_end_pos}")
        # Keep the TRAJ-TRAC junction from contig
        traj_part = tra_seq[:traj_end_pos]
        trac_part_from_contig = tra_seq[traj_end_pos:]
        
        # Ensure TRAC length matches template (use same number of nucleotides)
        trac_template_length = len(trac_template)
        
        if len(trac_part_from_contig) < trac_template_length:
            # If contig TRAC is shorter, pad with template
            trac_part = trac_part_from_contig + trac_template[len(trac_part_from_contig):].upper()
        elif len(trac_part_from_contig) > trac_template_length:
            # If contig TRAC is longer, truncate to template length
            trac_part = trac_part_from_contig[:trac_template_length]
        else:
            # Same length, use as-is
            trac_part = trac_part_from_contig
        
        tra_seq = traj_part + trac_part
        print(f"    Kept TRAJ-TRAC junction from contig, TRAC length: {len(trac_part)} (template: {trac_template_length})")
    else:
        print(f"    Warning: Could not find TRAJ end, using contig as-is")
        # Ensure TRAC length matches template if we can find it
        if len(tra_seq) > len(trac_template):
            # Try to find where TRAC might start (look for TRAC start)
            trac_start_region = trac_template[:min(30, len(trac_template))]
            for i in range(len(tra_seq) - len(trac_start_region), max(0, len(tra_seq) - 200), -1):
                contig_region = tra_seq[i:i+len(trac_start_region)]
                matches = sum(1 for a, b in zip(contig_region, trac_start_region) if a == b)
                if matches / len(trac_start_region) >= 0.80:
                    # Found TRAC start, adjust length
                    trac_part = tra_seq[i:i+trac_template_length] if i+trac_template_length <= len(tra_seq) else tra_seq[i:]
                    if len(trac_part) < trac_template_length:
                        trac_part = trac_part + trac_template[len(trac_part):].upper()
                    tra_seq = tra_seq[:i] + trac_part
                    break
    
    # Step 4: Remove restriction sites
    tra_seq, tra_mutations, tra_failed_sites = remove_restriction_sites(tra_seq)
    tra_mutations = tra_mutations or []
    tra_failed_sites = tra_failed_sites or []
    
    return tra_seq, tra_mutations, tra_failed_sites


def format_final_assembly(seq, trb_cdr3_nt, tra_cdr3_nt, trb_mutations=None, tra_mutations=None, component_info=None):
    """
    Format final TRB-P2A-TRA assembly with proper casing.
    Same as original script.
    """
    if not seq:
        return seq
    
    seq_list = list(seq.lower())
    
    asci_pos = 0
    kozak_pos = len(ASCI_SITE)
    
    p2a_path = os.path.join(TEMPLATES_DIR, "codon_optimized_P2A.fasta")
    p2a_seq = read_fasta(p2a_path)
    
    trb_start = kozak_pos + len(KOZAK_SEQ)
    if component_info:
        trb_length = component_info['trb_length']
        p2a_length = component_info['p2a_length']
        tra_length = component_info['tra_length']
    else:
        trb_length = 950
        p2a_length = len(p2a_seq) if p2a_seq else 60
        tra_length = 400
    
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
    Note: We don't remove the heading 0, 1, 2 positions of TRAC since we're using the full contig junction.
    """
    if not trb_seq or not tra_seq:
        print(f"      Error: Missing TRB or TRA sequence for assembly")
        return None, None
    
    p2a_path = os.path.join(TEMPLATES_DIR, "codon_optimized_P2A.fasta")
    p2a_seq = read_fasta(p2a_path)
    if not p2a_seq:
        print(f"      Error: P2A sequence not found at {p2a_path}")
        return None, None
    
    # Assemble: AscI + Kozak + TRB + P2A + TRA + PstI (all lowercase for now)
    final_seq = ASCI_SITE.lower() + KOZAK_SEQ.lower() + trb_seq.lower() + p2a_seq.lower() + tra_seq.lower() + PSTI_SITE.lower()
    
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


def generate_final_summary_log(output_dir, processed_tcrs):
    """
    Generate a final summary log file with all failed restriction sites and invalid sequences.
    
    Args:
        output_dir: Output directory
        processed_tcrs: List of tuples (tcr_id, trb_valid, tra_valid, trb_failed_sites, tra_failed_sites, trb_remaining_sites, tra_remaining_sites)
    """
    summary_log_file = os.path.join(output_dir, "final_validation_summary.txt")
    
    with open(summary_log_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("Final Validation Summary - All TCR Assemblies\n")
        f.write("=" * 80 + "\n\n")
        
        # Collect invalid sequences
        invalid_tcrs = []
        tcrs_with_failed_sites = []
        
        for tcr_id, trb_valid, tra_valid, trb_failed, tra_failed, trb_remaining, tra_remaining in processed_tcrs:
            is_invalid = not (trb_valid and tra_valid)
            has_failed_sites = len(trb_failed) > 0 or len(tra_failed) > 0 or len(trb_remaining) > 0 or len(tra_remaining) > 0
            
            if is_invalid:
                invalid_tcrs.append((tcr_id, trb_valid, tra_valid))
            if has_failed_sites:
                tcrs_with_failed_sites.append((tcr_id, trb_failed, tra_failed, trb_remaining, tra_remaining))
        
        # Section 1: Invalid TCR Sequences
        f.write("INVALID TCR ASSEMBLY SEQUENCES\n")
        f.write("=" * 80 + "\n")
        if invalid_tcrs:
            f.write(f"Total invalid sequences: {len(invalid_tcrs)}\n\n")
            for tcr_id, trb_valid, tra_valid in invalid_tcrs:
                f.write(f"TCR ID: {tcr_id}\n")
                f.write(f"  TRB valid: {'Yes' if trb_valid else 'No'}\n")
                f.write(f"  TRA valid: {'Yes' if tra_valid else 'No'}\n")
                f.write(f"  Validation log: {output_dir}/{tcr_id}/{tcr_id}_validation_log.txt\n")
                f.write("\n")
        else:
            f.write("✓ All TCR sequences are valid (no stop codons found)\n\n")
        
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Section 2: Failed Restriction Sites
        f.write("FAILED RESTRICTION SITE REMOVAL\n")
        f.write("=" * 80 + "\n")
        if tcrs_with_failed_sites:
            f.write(f"Total TCRs with failed restriction sites: {len(tcrs_with_failed_sites)}\n\n")
            for tcr_id, trb_failed, tra_failed, trb_remaining, tra_remaining in tcrs_with_failed_sites:
                f.write(f"TCR ID: {tcr_id}\n")
                
                if trb_failed:
                    f.write(f"  TRB - Failed to remove ({len(trb_failed)} sites):\n")
                    for site_name, pos, site_seq in trb_failed:
                        f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                
                if tra_failed:
                    f.write(f"  TRA - Failed to remove ({len(tra_failed)} sites):\n")
                    for site_name, pos, site_seq in tra_failed:
                        f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                
                if trb_remaining:
                    # Count sites (handle both list and single value formats)
                    trb_count = sum(len(v) if isinstance(v, list) else 1 for v in trb_remaining.values())
                    f.write(f"  TRB - Remaining sites ({trb_count} sites):\n")
                    for site_name, positions in trb_remaining.items():
                        site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
                        if isinstance(positions, list):
                            for pos in positions:
                                f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                        else:
                            pos = positions
                            f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                
                if tra_remaining:
                    # Count sites (handle both list and single value formats)
                    tra_count = sum(len(v) if isinstance(v, list) else 1 for v in tra_remaining.values())
                    f.write(f"  TRA - Remaining sites ({tra_count} sites):\n")
                    for site_name, positions in tra_remaining.items():
                        site_seq = ASCI_SITE if site_name == 'AscI' else PSTI_SITE
                        if isinstance(positions, list):
                            for pos in positions:
                                f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                        else:
                            pos = positions
                            f.write(f"    ✗ {site_name} at position {pos}: {site_seq}\n")
                
                f.write(f"  Validation log: {output_dir}/{tcr_id}/{tcr_id}_validation_log.txt\n")
                f.write("\n")
        else:
            f.write("✓ All restriction sites successfully removed from all TCRs\n\n")
        
        f.write("\n" + "=" * 80 + "\n\n")
        
        # Section 3: Summary Statistics
        f.write("SUMMARY STATISTICS\n")
        f.write("=" * 80 + "\n")
        f.write(f"Total TCRs processed: {len(processed_tcrs)}\n")
        f.write(f"Invalid sequences: {len(invalid_tcrs)}\n")
        f.write(f"TCRs with failed restriction sites: {len(tcrs_with_failed_sites)}\n")
        f.write(f"Valid sequences: {len(processed_tcrs) - len(invalid_tcrs)}\n")
        f.write(f"Success rate: {100 * (len(processed_tcrs) - len(invalid_tcrs)) / len(processed_tcrs):.1f}%\n")
        f.write("=" * 80 + "\n")


def process_tcr(row):
    """Process a single TCR entry from full contig sequences."""
    tcr_id = row['tcr_id']
    print(f"\nProcessing {tcr_id}...")
    
    # Create output directory for this TCR
    tcr_output_dir = os.path.join(OUTPUT_DIR, tcr_id)
    os.makedirs(tcr_output_dir, exist_ok=True)
    
    # Process TRB contig
    print(f"  Processing TRB contig...")
    trb_result = process_trb_contig(row)
    if not trb_result or trb_result[0] is None:
        print(f"    Error: Failed to process TRB contig")
        return False, False, False, [], [], {}, {}
    trb_seq, trb_mutations, trb_failed_sites = trb_result
    
    # Process TRA contig
    print(f"  Processing TRA contig...")
    tra_result = process_tra_contig(row)
    if not tra_result or tra_result[0] is None:
        print(f"    Error: Failed to process TRA contig")
        return False, False, False, [], [], {}, {}
    tra_seq, tra_mutations, tra_failed_sites = tra_result
    
    # Get CDR3 sequences for formatting
    trb_cdr3_nt = row.get('TRB_cdr3_nt', '')
    tra_cdr3_nt = row.get('TRA_cdr3_nt', '')
    
    # Check for remaining restriction sites after processing
    trb_remaining_sites = check_restriction_sites(trb_seq)
    tra_remaining_sites = check_restriction_sites(tra_seq)
    
    # Validate sequences by translation
    print(f"  Validating sequences (translation and stop codon check)...")
    trb_valid, tra_valid = validate_and_log_translation(trb_seq, tra_seq, tcr_id, tcr_output_dir, trb_mutations, tra_mutations, trb_failed_sites, tra_failed_sites)
    
    if not trb_valid:
        print(f"    Warning: TRB sequence validation failed - check log file")
    if not tra_valid:
        print(f"    Warning: TRA sequence validation failed - check log file")
    
    if trb_valid and tra_valid:
        print(f"    ✓ Both sequences validated successfully")
    else:
        print(f"    ⚠ Validation issues detected - see {tcr_id}_validation_log.txt")
    
    # Assemble final TRB-P2A-TRA sequence
    print(f"  Assembling TRB-P2A-TRA...")
    final_seq, comp_info = assemble_trb_p2a_tra(trb_seq, tra_seq, trb_mutations, tra_mutations, trb_cdr3_nt, tra_cdr3_nt)
    
    if final_seq:
        formatted_final = format_final_assembly(final_seq, trb_cdr3_nt, tra_cdr3_nt, trb_mutations, tra_mutations, comp_info)
        output_file = os.path.join(tcr_output_dir, "TRB_P2A_TRA_original_contig.fasta")
        save_sequence(formatted_final, output_file)
        print(f"    Saved: {tcr_id}/TRB_P2A_TRA_original_contig.fasta")
        print(f"    Validation log: {tcr_id}_validation_log.txt")
        return True, trb_valid, tra_valid, trb_failed_sites, tra_failed_sites, trb_remaining_sites, tra_remaining_sites
    else:
        print(f"    Error: Failed to assemble final sequence")
        return False, trb_valid, tra_valid, trb_failed_sites, tra_failed_sites, trb_remaining_sites, tra_remaining_sites


def main(input_file=None, output_dir=None, max_tcrs=None):
    """
    Main function.
    
    Args:
        input_file: Path to input CSV file (default: uses global INPUT_FILE)
        output_dir: Path to output directory (default: uses global OUTPUT_DIR)
        max_tcrs: Maximum number of TCRs to process (default: uses global MAX_TCRS)
    """
    # Declare global variables first
    global OUTPUT_DIR
    
    # Use provided parameters or fall back to global variables
    input_path = input_file if input_file is not None else INPUT_FILE
    output_path = output_dir if output_dir is not None else OUTPUT_DIR
    max_tcrs_limit = max_tcrs if max_tcrs is not None else MAX_TCRS
    
    # Update global OUTPUT_DIR for use in process_tcr function
    OUTPUT_DIR = output_path
    
    print("=" * 60)
    print("TCR Sequence Assembly from Full Contig Sequences")
    print("=" * 60)
    
    # Load input file
    print(f"\nLoading input file: {input_path}")
    df = pd.read_csv(input_path)
    print(f"Loaded {len(df)} TCR entries")
    
    # Filter entries that have full contig sequences
    df_with_contigs = df[
        df['TRA_full_contig'].notna() & (df['TRA_full_contig'] != '') &
        df['TRB_full_contig'].notna() & (df['TRB_full_contig'] != '')
    ]
    print(f"Entries with full contig sequences: {len(df_with_contigs)}")
    
    # Limit to MAX_TCRS if specified
    if max_tcrs_limit is not None and max_tcrs_limit > 0:
        df_with_contigs = df_with_contigs.head(max_tcrs_limit)
        print(f"Processing first {len(df_with_contigs)} TCR entries (MAX_TCRS = {max_tcrs_limit})")
    else:
        print(f"Processing all {len(df_with_contigs)} TCR entries")
    
    # Process each TCR
    successful = 0
    failed = 0
    processed_tcrs = []  # Store validation info for summary log
    
    for idx, row in df_with_contigs.iterrows():
        try:
            result = process_tcr(row)
            if isinstance(result, tuple):
                success, trb_valid, tra_valid, trb_failed, tra_failed, trb_remaining, tra_remaining = result
                if success:
                    successful += 1
                else:
                    failed += 1
                # Store validation info for summary log
                tcr_id = row['tcr_id']
                processed_tcrs.append((tcr_id, trb_valid, tra_valid, trb_failed, tra_failed, trb_remaining, tra_remaining))
            else:
                # Old return format (shouldn't happen, but handle gracefully)
                if result:
                    successful += 1
                else:
                    failed += 1
        except Exception as e:
            print(f"Error processing {row.get('tcr_id', 'unknown')}: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    # Generate final summary log
    if processed_tcrs:
        print(f"\nGenerating final validation summary log...")
        generate_final_summary_log(output_path, processed_tcrs)
        print(f"  Saved: {output_path}/final_validation_summary.txt")
    
    print(f"\n{'='*60}")
    print(f"Summary: {successful} successful, {failed} failed")
    print(f"Output directory: {output_path}")
    print(f"Final summary log: {output_path}/final_validation_summary.txt")
    print("=" * 60)
    
    return successful, failed


if __name__ == "__main__":
    main()

