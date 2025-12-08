#!/usr/bin/env python3
"""
Script to add InterPro domain counts to fungal reference proteomes data.

This script reads proteome data from a TSV file, queries the UniProt API to count
proteins containing a specific InterPro domain for each proteome, and outputs
the results to a new TSV file with the counts added.

Usage:
    python add_interpro_counts.py

Requirements:
    - requests library (pip install requests)
    - fungal_ref_proteomes.tsv file in the same directory

Output:
    - Creates fungal_ref_proteomes_with_interpro.tsv with added InterPro counts

Author: Generated script for fungal proteome InterPro analysis
"""

import csv
import time
import requests
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# === CONFIGURATION ===
INPUT_TSV = "fungal_ref_proteomes.tsv"
OUTPUT_TSV = "fungal_ref_proteomes_with_interpro.tsv"
INTERPRO_ID = "IPR013902"  # 👈 REPLACE with your InterPro domain ID
MAX_WORKERS = 10  # Number of concurrent requests (adjust based on your needs)
# 📊 MAX_WORKERS Tuning Guide:
#   - Conservative (5-10): Safe start, good for testing
#   - Moderate (15-20): Faster processing, balanced approach  
#   - Aggressive (25-50): Maximum speed, watch for rate limits
# 💡 For large datasets (1000+ proteomes): try 15-25 workers
# 🚨 Reduce if you see HTTP 429 errors or timeouts

# UniProt API base URL
BASE_URL = "https://rest.uniprot.org/uniprotkb/search"

# Global session for connection reuse
session = requests.Session()
session.headers.update({'User-Agent': 'InterPro-Counter/1.0'})

# Thread-safe lock for progress reporting
progress_lock = Lock()


def count_proteins_with_interpro(proteome_id, interpro_id):
    """
    Returns the number of proteins in a proteome that have the given InterPro domain.
    Uses UniProt's count endpoint for efficiency.
    
    Args:
        proteome_id (str): UniProt proteome ID (upid)
        interpro_id (str): InterPro domain ID to search for
        
    Returns:
        int: Number of proteins containing the InterPro domain
        
    Note:
        Uses the x-total-results header from UniProt API for efficient counting
        without downloading all results.
    """
    # Use the search endpoint with size=0 to get total number from headers
    url = f"https://rest.uniprot.org/uniprotkb/search?query=proteome:{proteome_id} AND database:interpro AND {interpro_id}&format=tsv&size=0"
    
    try:
        response = session.get(url)
        if response.status_code == 200:
            # The count is in the 'x-total-results' header
            return int(response.headers.get('x-total-results', 0))
        else:
            print(f"⚠️ Warning: Failed to fetch count for proteome {proteome_id} (status {response.status_code})")
            return 0
    except requests.RequestException as e:
        print(f"⚠️ Warning: Network error for proteome {proteome_id}: {e}")
        return 0
    except ValueError as e:
        print(f"⚠️ Warning: Invalid response for proteome {proteome_id}: {e}")
        return 0


def process_single_proteome(row_data, index, total_count):
    """
    Process a single proteome and return the updated row with InterPro count.
    
    Args:
        row_data (dict): Row data from the input TSV
        index (int): Current row index for progress reporting
        total_count (int): Total number of proteomes to process
        
    Returns:
        dict: Updated row with InterPro count added
    """
    proteome_id = row_data['Proteome Id']
    organism = row_data.get('Organism', 'N/A')
    
    # Thread-safe progress reporting
    with progress_lock:
        print(f"Processing {index}/{total_count}: {proteome_id} ({organism})")
    
    count = count_proteins_with_interpro(proteome_id, INTERPRO_ID)
    row_data[f'interpro_{INTERPRO_ID}_count'] = count
    
    # Small delay to be respectful to the server
    time.sleep(0.05)
    
    return row_data, count >= 0  # Return row and success status


def validate_input_file(filename):
    """
    Validate that the input TSV file exists and has the expected format.
    
    Args:
        filename (str): Path to the input TSV file
        
    Returns:
        bool: True if file is valid, False otherwise
    """
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='\t')
            if 'Proteome Id' not in reader.fieldnames:
                print(f"❌ Error: Input file '{filename}' must contain 'Proteome Id' column")
                return False
        return True
    except FileNotFoundError:
        print(f"❌ Error: Input file '{filename}' not found")
        return False
    except Exception as e:
        print(f"❌ Error: Failed to read input file '{filename}': {e}")
        return False


def main():
    """
    Main function to process proteomes and add InterPro domain counts using concurrent processing.
    
    Reads the input TSV file, queries UniProt API for proteomes concurrently,
    and writes results to output TSV file with added InterPro counts.
    """
    print("=== Fungal Proteome InterPro Domain Counter (Concurrent Version) ===")
    print(f"Input file: {INPUT_TSV}")
    print(f"Output file: {OUTPUT_TSV}")
    print(f"InterPro domain: {INTERPRO_ID}")
    print(f"Max concurrent workers: {MAX_WORKERS}")
    print()
    
    # Validate input file
    if not validate_input_file(INPUT_TSV):
        sys.exit(1)
    
    try:
        # First pass: read all data into memory
        print("📖 Reading input data...")
        rows_data = []
        with open(INPUT_TSV, mode='r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile, delimiter='\t')
            fieldnames = reader.fieldnames + [f'interpro_{INTERPRO_ID}_count']
            rows_data = list(reader)
        
        total_count = len(rows_data)
        print(f"📊 Found {total_count} proteomes to process")
        print("🚀 Starting concurrent processing...\n")
        
        # Process proteomes concurrently
        processed_rows = [None] * total_count  # Maintain order
        successful_queries = 0
        
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            # Submit all tasks
            future_to_index = {
                executor.submit(process_single_proteome, row, i+1, total_count): i 
                for i, row in enumerate(rows_data)
            }
            
            # Collect results as they complete
            for future in as_completed(future_to_index):
                index = future_to_index[future]
                try:
                    processed_row, success = future.result()
                    processed_rows[index] = processed_row
                    if success:
                        successful_queries += 1
                except Exception as e:
                    print(f"❌ Error processing proteome at index {index+1}: {e}")
                    # Keep original row without InterPro count
                    processed_rows[index] = rows_data[index]
                    processed_rows[index][f'interpro_{INTERPRO_ID}_count'] = -1
        
        # Write results to output file
        print(f"\n💾 Writing results to {OUTPUT_TSV}...")
        with open(OUTPUT_TSV, mode='w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in processed_rows:
                if row is not None:
                    writer.writerow(row)

        print(f"\n✅ Done! Output written to {OUTPUT_TSV}")
        print(f"📊 Summary:")
        print(f"   - Total proteomes processed: {total_count}")
        print(f"   - Successful queries: {successful_queries}")
        print(f"   - Failed queries: {total_count - successful_queries}")
        
    except Exception as e:
        print(f"❌ Error during processing: {e}")
        sys.exit(1)


if __name__ == "__main__":
    # Help message
    if len(sys.argv) > 1 and sys.argv[1] in ["-h", "--help", "help"]:
        print(__doc__)
        sys.exit(0)
    
    main()