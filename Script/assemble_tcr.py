#!/usr/bin/env python3
"""
Master script for TCR sequence assembly.

This script can automatically determine which assembly method to use based on the input file,
or the user can explicitly specify which method(s) to use:
- If TRA_full_contig and TRB_full_contig columns exist: uses assemble_from_contig.py
- Otherwise: uses assemble_IMGT_tcr.py

When running IMGT template-based assembly, the script also automatically runs
test_assemble_IMGT_tcr.py after the main assembly.

Outputs are organized in subdirectories:
- Contig_tcr/ for contig-based assembly
- IMGT_tcr/ for IMGT template-based assembly
"""

import os
import sys
import argparse
import pandas as pd
from pathlib import Path

# Import the assembly scripts
# Add Script directory to path for imports
script_dir = os.path.dirname(os.path.abspath(__file__))
if script_dir not in sys.path:
    sys.path.insert(0, script_dir)

from assemble_from_contig import main as assemble_from_contig_main
from assemble_IMGT_tcr import main as assemble_IMGT_tcr_main
from test_assemble_IMGT_tcr import main as test_assemble_IMGT_tcr_main


def check_contig_columns(input_file):
    """
    Check if the input file has TRA_full_contig and TRB_full_contig columns.
    
    Args:
        input_file: Path to input CSV file
    
    Returns:
        bool: True if both contig columns exist and have data, False otherwise
    """
    try:
        df = pd.read_csv(input_file, nrows=1)  # Read just the header
        has_tra_contig = 'TRA_full_contig' in df.columns
        has_trb_contig = 'TRB_full_contig' in df.columns
        
        if has_tra_contig and has_trb_contig:
            # Check if there's actual data in these columns
            df_full = pd.read_csv(input_file)
            tra_has_data = df_full['TRA_full_contig'].notna().any() and (df_full['TRA_full_contig'] != '').any()
            trb_has_data = df_full['TRB_full_contig'].notna().any() and (df_full['TRB_full_contig'] != '').any()
            return tra_has_data and trb_has_data
        
        return False
    except Exception as e:
        print(f"Error checking input file: {e}")
        return False


def get_output_directory(input_file, output_path=None):
    """
    Determine output directory based on input file name.
    
    Args:
        input_file: Path to input CSV file
        output_path: Optional custom output path
    
    Returns:
        str: Output directory path
    """
    if output_path:
        return output_path
    
    # Extract input filename without extension
    input_basename = Path(input_file).stem
    output_dir = os.path.join("Output", input_basename)
    return output_dir


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description='TCR Sequence Assembly Master Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Use default input file and output directory (auto-detect method)
  python Script/assemble_tcr.py
  
  # Specify input file
  python Script/assemble_tcr.py -i Input/my_data.csv
  
  # Specify assembly method explicitly
  python Script/assemble_tcr.py -i Input/my_data.csv -m imgt
  python Script/assemble_tcr.py -i Input/my_data.csv -m contig
  python Script/assemble_tcr.py -i Input/my_data.csv -m both
  
  # Specify input file and output directory
  python Script/assemble_tcr.py -i Input/my_data.csv -o Output/my_results
  
  # Process only first 10 TCRs
  python Script/assemble_tcr.py -i Input/my_data.csv -n 10
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        type=str,
        default='Input/sample_input.csv',
        help='Path to input CSV file (default: Input/sample_input.csv)'
    )
    
    parser.add_argument(
        '-o', '--output',
        type=str,
        default=None,
        help='Path to output directory (default: Output/{input_filename})'
    )
    
    parser.add_argument(
        '-n', '--num-tcrs',
        type=int,
        default=None,
        help='Number of TCRs to process (default: all TCRs)'
    )
    
    parser.add_argument(
        '-m', '--method',
        type=str,
        choices=['imgt', 'contig', 'both'],
        default=None,
        help='Assembly method to use: "imgt" for IMGT template-based, "contig" for contig-based, "both" for both methods. If not specified, automatically detects based on input file.'
    )
    
    args = parser.parse_args()
    
    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Determine output directory
    output_base = get_output_directory(args.input, args.output)
    
    # Determine which assembly method(s) to use
    print("=" * 80)
    print("TCR Sequence Assembly Master Script")
    print("=" * 80)
    print(f"\nInput file: {args.input}")
    print(f"Output base directory: {output_base}")
    print(f"Number of TCRs: {args.num_tcrs if args.num_tcrs else 'All'}")
    print()
    
    # Determine assembly method
    if args.method:
        # User specified method explicitly
        run_imgt = args.method in ['imgt', 'both']
        run_contig = args.method in ['contig', 'both']
        print(f"Assembly method specified: {args.method}")
    else:
        # Auto-detect based on input file
        print("Checking input file for contig sequences...")
        has_contigs = check_contig_columns(args.input)
        run_imgt = not has_contigs
        run_contig = has_contigs
        if has_contigs:
            print("✓ Found TRA_full_contig and TRB_full_contig columns")
        else:
            print("✗ TRA_full_contig and/or TRB_full_contig columns not found")
    
    # Run IMGT assembly if needed
    if run_imgt:
        print("  Using IMGT template-based assembly (assemble_IMGT_tcr.py)")
        output_dir = os.path.join(output_base, "IMGT_tcr")
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Output directory: {output_dir}\n")
        
        # Run IMGT template-based assembly
        try:
            successful, failed = assemble_IMGT_tcr_main(
                input_file=args.input,
                output_dir=output_dir,
                max_tcrs=args.num_tcrs
            )
            print(f"\n{'='*80}")
            print("IMGT template-based assembly completed!")
            print(f"  Successful: {successful}")
            print(f"  Failed: {failed}")
            print(f"  Output: {output_dir}")
            print("=" * 80)
            
            # Run test script after IMGT assembly
            print(f"\n{'='*80}")
            print("Running test assembly (test_assemble_IMGT_tcr.py)...")
            print("=" * 80)
            try:
                test_successful, test_failed = test_assemble_IMGT_tcr_main(
                    input_file=args.input,
                    output_dir=output_dir,
                    max_tcrs=args.num_tcrs
                )
                print(f"\n{'='*80}")
                print("Test assembly completed!")
                print(f"  Successful: {test_successful}")
                print(f"  Failed: {test_failed}")
                print("=" * 80)
            except Exception as e:
                print(f"\nError during test assembly: {e}")
                import traceback
                traceback.print_exc()
                # Don't exit, continue with contig assembly if needed
        except Exception as e:
            print(f"\nError during IMGT template-based assembly: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)
    
    # Run contig assembly if needed
    if run_contig:
        print("  Using contig-based assembly (assemble_from_contig.py)")
        output_dir = os.path.join(output_base, "Contig_tcr")
        os.makedirs(output_dir, exist_ok=True)
        print(f"  Output directory: {output_dir}\n")
        
        # Run contig-based assembly
        try:
            successful, failed = assemble_from_contig_main(
                input_file=args.input,
                output_dir=output_dir,
                max_tcrs=args.num_tcrs
            )
            print(f"\n{'='*80}")
            print("Contig-based assembly completed!")
            print(f"  Successful: {successful}")
            print(f"  Failed: {failed}")
            print(f"  Output: {output_dir}")
            print("=" * 80)
        except Exception as e:
            print(f"\nError during contig-based assembly: {e}")
            import traceback
            traceback.print_exc()
            sys.exit(1)


if __name__ == "__main__":
    main()

