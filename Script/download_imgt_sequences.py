#!/usr/bin/env python3
"""
Script to download TRBV, TRBJ, TRAV, TRAJ fragments from IMGT and organize into folders.
"""

import os
import re
import time
import requests
from bs4 import BeautifulSoup
from urllib.parse import urljoin, urlparse
import sys

# Base IMGT URL
IMGT_BASE_URL = "https://www.imgt.org"

# Gene group URLs
GENE_GROUPS = {
    "TRAV": "https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?species=Homo%20sapiens&group=TRAV",
    "TRBV": "https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?species=Homo%20sapiens&group=TRBV",
    "TRAJ": "https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?species=Homo%20sapiens&group=TRAJ",
    "TRBJ": "https://www.imgt.org/IMGTrepertoire/LocusGenes/genetable/autotable.php?species=Homo%20sapiens&group=TRBJ"
}

# Output format for each gene group
# TRAV and TRBV use fastaNConstructs, TRAJ and TRBJ use fastaNAllP
OUTPUT_FORMATS = {
    "TRAV": "fastaNConstructs",
    "TRBV": "fastaNConstructs",
    "TRAJ": "fastaNAllP",
    "TRBJ": "fastaNAllP"
}

# Output base directory
OUTPUT_BASE = "Templates"

# Session for maintaining cookies
session = requests.Session()
session.headers.update({
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
})


def create_output_folders():
    """Create output folders for each gene group."""
    for group in GENE_GROUPS.keys():
        folder_path = os.path.join(OUTPUT_BASE, group)
        os.makedirs(folder_path, exist_ok=True)
        print(f"Created/verified folder: {folder_path}")


def get_gene_names_from_table(url):
    """
    Extract all gene names from the gene table page.
    Returns a list of gene names (strings).
    """
    print(f"\nFetching gene table from: {url}")
    try:
        response = session.get(url, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        gene_names = []
        
        # IMGT gene tables typically have gene names as links in tables
        # Look for all tables first
        tables = soup.find_all('table')
        
        for table in tables:
            # Find all rows in the table
            rows = table.find_all('tr')
            for row in rows:
                # Find all links in this row
                links = row.find_all('a', href=True)
                for link in links:
                    gene_name = link.get_text(strip=True)
                    
                    # Filter for valid gene names (typically start with TRAV, TRBV, TRAJ, TRBJ)
                    # and look like gene identifiers (e.g., TRAV1-1, TRBV2-1, TRAV38-2/DV8, etc.)
                    # Pattern: TRAV1, TRAV1-1, TRBV2-1, TRAV38-2/DV8, TRAV14/DV4, etc.
                    # Includes support for /DV variants
                    if gene_name and re.match(r'^(TRAV|TRBV|TRAJ|TRBJ)[\d]+(-[\d]+)?(/DV[\d]+)?(\*[\d]+)?$', gene_name):
                        gene_names.append(gene_name)
        
        # If no links found in tables, try broader search
        if not gene_names:
            for link in soup.find_all('a', href=True):
                gene_name = link.get_text(strip=True)
                # Include support for /DV variants
                if gene_name and re.match(r'^(TRAV|TRBV|TRAJ|TRBJ)[\d]+(-[\d]+)?(/DV[\d]+)?(\*[\d]+)?$', gene_name):
                    gene_names.append(gene_name)
        
        # Remove duplicates while preserving order
        seen = set()
        unique_names = []
        dv_genes = []
        for name in gene_names:
            if name not in seen:
                seen.add(name)
                unique_names.append(name)
                # Track genes with /DV
                if '/DV' in name:
                    dv_genes.append(name)
        
        print(f"Found {len(unique_names)} unique gene names")
        if dv_genes:
            print(f"  Found {len(dv_genes)} genes with /DV: {', '.join(dv_genes[:10])}" + 
                  (f" ... and {len(dv_genes) - 10} more" if len(dv_genes) > 10 else ""))
        return unique_names
        
    except Exception as e:
        print(f"Error fetching gene table: {e}")
        import traceback
        traceback.print_exc()
        return []


def download_fasta_from_api(gene_name, output_format="fastaNConstructs"):
    """
    Download FASTA sequence directly from IMGT FASTA API.
    URL format: https://www.imgt.org/genedb/fasta?outputFormat={output_format}&selectGenes={gene_name}&selectSpecies=Homo%20sapiens
    Returns the FASTA text or None.
    
    Args:
        gene_name: The gene name (e.g., TRAV1, TRBV1, TRAJ1, TRBJ1, TRAV38-2/DV8)
        output_format: The output format - "fastaNConstructs" for TRAV/TRBV, "fastaNAllP" for TRAJ/TRBJ
    """
    # URL encode the gene name to handle special characters like "/"
    from urllib.parse import quote
    encoded_gene_name = quote(gene_name, safe='')
    
    # Construct the FASTA API URL
    fasta_url = f"https://www.imgt.org/genedb/fasta?outputFormat={output_format}&selectGenes={encoded_gene_name}&selectSpecies=Homo%20sapiens"
    
    print(f"  Downloading FASTA from API: {fasta_url}")
    try:
        response = session.get(fasta_url, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Look for FASTA text in <pre> tags (most common format)
        pre_tags = soup.find_all('pre')
        for pre in pre_tags:
            text = pre.get_text()
            if text.strip().startswith('>'):
                return text.strip()
        
        # Alternative: Look for FASTA text in <textarea> tags
        textareas = soup.find_all('textarea')
        for textarea in textareas:
            text = textarea.get_text()
            if text.strip().startswith('>'):
                return text.strip()
        
        # Alternative: Extract FASTA from page text
        # FASTA sequences start with '>' and contain nucleotide sequences
        body_text = soup.get_text()
        lines = body_text.split('\n')
        fasta_lines = []
        in_fasta = False
        
        for line in lines:
            stripped = line.strip()
            if stripped.startswith('>'):
                in_fasta = True
                fasta_lines = [stripped]
            elif in_fasta:
                # Continue collecting FASTA lines (nucleotide sequences)
                if stripped and not stripped.startswith('<') and not stripped.startswith('|'):
                    # Check if it looks like a nucleotide sequence (A, T, G, C, N, and whitespace)
                    if re.match(r'^[ATGCNatgcn\s]+$', stripped):
                        fasta_lines.append(stripped)
                    elif stripped.startswith('>'):
                        # New FASTA entry, save previous and start new
                        if len(fasta_lines) > 1:
                            return '\n'.join(fasta_lines)
                        fasta_lines = [stripped]
                elif stripped.startswith('Number of results') or stripped.startswith('IMGT'):
                    # End of FASTA section
                    if len(fasta_lines) > 1:
                        break
        
        if fasta_lines and len(fasta_lines) > 1:
            return '\n'.join(fasta_lines)
        
        print(f"    Warning: Could not extract FASTA text from API response")
        return None
        
    except Exception as e:
        print(f"    Error downloading FASTA from API: {e}")
        import traceback
        traceback.print_exc()
        return None


def save_fasta(fasta_text, output_path):
    """Save FASTA text to a file."""
    try:
        if not fasta_text or not fasta_text.strip():
            print(f"    Error: Empty FASTA text")
            return False
        
        # Ensure it starts with > (FASTA format)
        if not fasta_text.strip().startswith('>'):
            print(f"    Warning: Text doesn't appear to be in FASTA format")
        
        # Save the file
        with open(output_path, 'w') as f:
            f.write(fasta_text.strip())
            f.write('\n')  # Ensure file ends with newline
        
        print(f"    Saved: {output_path}")
        return True
        
    except Exception as e:
        print(f"    Error saving FASTA: {e}")
        return False


def process_gene_group(group_name, group_url, skip_if_complete=False):
    """
    Process a single gene group: get all genes and download their FASTA files.
    
    Args:
        group_name: Name of the gene group (TRAV, TRBV, TRAJ, TRBJ)
        group_url: URL to the gene table page
        skip_if_complete: If True, skip this group if output folder already has files
    """
    print(f"\n{'='*60}")
    print(f"Processing {group_name}")
    print(f"{'='*60}")
    
    # Create output folder
    output_folder = os.path.join(OUTPUT_BASE, group_name)
    os.makedirs(output_folder, exist_ok=True)
    
    # Get all gene names from the table
    gene_names = get_gene_names_from_table(group_url)
    
    if not gene_names:
        print(f"No gene names found for {group_name}")
        return
    
    # Check which files already exist (handle both / and _ in filenames)
    existing_files = set()
    existing_safe_names = set()
    if os.path.exists(output_folder):
        for f in os.listdir(output_folder):
            if f.endswith('.fasta'):
                safe_name = f[:-6]  # Remove .fasta
                existing_safe_names.add(safe_name)
                # Convert back to original format for comparison
                original_name = safe_name.replace('_', '/')
                existing_files.add(original_name)
    
    # Always check for /DV genes and download if missing
    dv_genes = [g for g in gene_names if '/DV' in g]
    if dv_genes:
        # Get unique /DV gene names (without allele info)
        unique_dv_genes = []
        seen_dv = set()
        for g in dv_genes:
            # Remove allele info (*01, *02, etc.) for comparison
            base_name = g.split('*')[0]
            if base_name not in seen_dv:
                seen_dv.add(base_name)
                unique_dv_genes.append(base_name)
        
        # Check which /DV genes are missing
        missing_dv = [g for g in unique_dv_genes if g not in existing_files]
        
        if missing_dv:
            print(f"\nFound {len(unique_dv_genes)} unique /DV genes, {len(missing_dv)} are missing")
            print(f"Missing /DV genes: {missing_dv}")
            # Download missing /DV genes (use the full name with allele if available)
            genes_to_download = []
            for missing in missing_dv:
                # Find the first occurrence in gene_names (prefer without allele)
                for g in gene_names:
                    if g.split('*')[0] == missing:
                        genes_to_download.append(g)
                        break
        else:
            print(f"\nAll {len(unique_dv_genes)} /DV genes already exist")
            if skip_if_complete:
                # Check if all genes are present
                missing_genes = [g for g in gene_names if g.split('*')[0] not in existing_files]
                if not missing_genes:
                    print(f"Skipping {group_name} - all files already exist")
                    return
                else:
                    genes_to_download = missing_genes
            else:
                genes_to_download = gene_names
    else:
        # No /DV genes, proceed normally
        if skip_if_complete and existing_files:
            missing_genes = [g for g in gene_names if g.split('*')[0] not in existing_files]
            if not missing_genes:
                print(f"Skipping {group_name} - all {len(existing_files)} files already exist")
                return
            else:
                genes_to_download = missing_genes
        else:
            genes_to_download = gene_names
    
    # Get the appropriate output format for this gene group
    output_format = OUTPUT_FORMATS.get(group_name, "fastaNConstructs")
    print(f"Using output format: {output_format}")
    print(f"Genes to download: {len(genes_to_download)}")
    
    downloaded = 0
    failed = 0
    
    for gene_name in genes_to_download:
        print(f"\nProcessing {gene_name}...")
        
        # Download FASTA directly from the API with the correct output format
        fasta_text = download_fasta_from_api(gene_name, output_format)
        
        if fasta_text:
            # Create output filename - replace "/" with "_" for filesystem compatibility
            # (e.g., "TRAV38-2/DV8" becomes "TRAV38-2_DV8.fasta")
            safe_gene_name = gene_name.replace("/", "_")
            output_filename = f"{safe_gene_name}.fasta"
            output_path = os.path.join(output_folder, output_filename)
            
            # Save the FASTA file
            if save_fasta(fasta_text, output_path):
                downloaded += 1
            else:
                failed += 1
        else:
            print(f"  Failed to download FASTA for {gene_name}")
            failed += 1
        
        # Be polite - add a small delay
        time.sleep(1)
    
    print(f"\n{group_name} Summary: {downloaded} downloaded, {failed} failed")


def main(skip_groups=None):
    """
    Main function to download all gene groups.
    
    Args:
        skip_groups: List of group names to skip (e.g., ['TRAV', 'TRBV'])
    """
    print("IMGT FASTA Downloader")
    print("=" * 60)
    
    if skip_groups:
        print(f"Skipping groups: {', '.join(skip_groups)}")
    
    # Create output folders
    create_output_folders()
    
    # Process each gene group
    for group_name, group_url in GENE_GROUPS.items():
        # Skip if in skip list
        if skip_groups and group_name in skip_groups:
            print(f"\nSkipping {group_name} (in skip list)")
            continue
        
        try:
            # Skip if folder already has files (only for groups not explicitly skipped)
            skip_if_complete = True  # Always check if complete to avoid re-downloading
            process_gene_group(group_name, group_url, skip_if_complete=skip_if_complete)
        except Exception as e:
            print(f"\nError processing {group_name}: {e}")
            continue
    
    print("\n" + "=" * 60)
    print("Download complete!")
    print("=" * 60)


if __name__ == "__main__":
    # For this run, skip TRAV and TRBV since they're already completed
    # Set skip_groups to None to process all groups, or pass a list to skip specific ones
    skip_groups = ['TRAV', 'TRBV']  # Skip completed groups
    main(skip_groups=skip_groups)


if __name__ == "__main__":
    main()

