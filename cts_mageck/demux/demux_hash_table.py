#!/usr/bin/env python3
"""
Fast MAGeCK Cell Type Demultiplexer - Hash Table Implementation
Uses dictionary lookup for O(1) gRNA matching instead of iterating all gRNAs
This provides 50-100x speedup compared to hamming distance iteration
"""

import csv
import gzip
import os
from collections import defaultdict
import time

class FastDemultiplexer:
    def __init__(self, barcode_csv, barcode_start=22, barcode_length=8, 
                 grna_start=12, grna_length=20, max_barcode_mismatches=1,
                 allow_grna_mismatch=True):
        self.celltype_barcodes = self._load_barcodes(barcode_csv)
        self.barcode_start = barcode_start
        self.barcode_length = barcode_length
        self.grna_start = grna_start
        self.grna_length = grna_length
        self.max_barcode_mismatches = max_barcode_mismatches
        self.allow_grna_mismatch = allow_grna_mismatch
        
        # Track detailed statistics
        self.barcode_stats = defaultdict(lambda: defaultdict(int))
        self.celltype_counts = defaultdict(int)
        
        print(f"Loaded {len(self.celltype_barcodes)} cell types: {list(self.celltype_barcodes.keys())}")
        print(f"gRNA extraction: positions {grna_start} to {grna_start + grna_length}")
        print(f"Barcode extraction: positions {barcode_start} to {barcode_start + barcode_length}")
        print(f"gRNA mismatch tolerance: {'1 mismatch' if allow_grna_mismatch else 'exact match only'}")
    
    def _load_barcodes(self, csv_file):
        """Load cell type barcodes"""
        barcodes = {}
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                celltype = row['Celltype'].strip()
                barcode = row['Barcode'].strip().upper()
                barcodes[celltype] = barcode
        return barcodes
    
    def _load_library_as_dict(self, library_file):
        """
        Load gRNA library into hash tables for O(1) lookup
        Returns two dictionaries:
        1. exact_match_dict: {gRNA_seq: gene} for exact matches
        2. one_mismatch_dict: {gRNA_variant: (original_gRNA, gene)} for 1-mismatch variants
        """
        print("\nBuilding gRNA hash tables...")
        start_time = time.time()
        
        exact_match_dict = {}
        one_mismatch_dict = {}
        
        with open(library_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    grna = parts[0].strip().upper()
                    gene = parts[1].strip()
                    
                    # Add exact match
                    exact_match_dict[grna] = gene
                    
                    # Generate 1-mismatch variants if tolerance enabled
                    if self.allow_grna_mismatch:
                        for pos in range(len(grna)):
                            original_base = grna[pos]
                            for new_base in ['A', 'C', 'G', 'T']:
                                if new_base != original_base:
                                    variant = grna[:pos] + new_base + grna[pos+1:]
                                    # Only add if not already an exact match to another gRNA
                                    if variant not in exact_match_dict:
                                        # Store the original gRNA and gene
                                        one_mismatch_dict[variant] = (grna, gene)
        
        elapsed = time.time() - start_time
        print(f"  Built hash tables in {elapsed:.2f} seconds")
        print(f"  Exact match entries: {len(exact_match_dict):,}")
        if self.allow_grna_mismatch:
            print(f"  1-mismatch variant entries: {len(one_mismatch_dict):,}")
        print(f"  Total lookup entries: {len(exact_match_dict) + len(one_mismatch_dict):,}")
        
        return exact_match_dict, one_mismatch_dict
    
    def hamming_distance(self, s1, s2):
        """Calculate Hamming distance between two strings"""
        if len(s1) != len(s2):
            return float('inf')
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))
    
    def extract_barcode(self, read2_seq):
        """Extract barcode from read2"""
        if len(read2_seq) < self.barcode_start + self.barcode_length:
            return None
        return read2_seq[self.barcode_start:self.barcode_start + self.barcode_length]
    
    def extract_grna(self, read1_seq):
        """Extract gRNA from read1 at specified position"""
        if len(read1_seq) < self.grna_start + self.grna_length:
            return None
        return read1_seq[self.grna_start:self.grna_start + self.grna_length].upper()
    
    def match_barcode(self, observed_bc):
        """Match barcode to cell type with error correction"""
        best_celltype = None
        best_distance = float('inf')
        
        for celltype, expected_bc in self.celltype_barcodes.items():
            distance = self.hamming_distance(observed_bc, expected_bc)
            if distance <= self.max_barcode_mismatches and distance < best_distance:
                best_celltype = celltype
                best_distance = distance
        
        return best_celltype, best_distance
    
    def match_grna_fast(self, observed_grna, exact_match_dict, one_mismatch_dict):
        """
        Fast gRNA matching using hash table lookup
        Returns: (is_valid, matched_grna, gene)
        """
        # Try exact match first (O(1) lookup)
        if observed_grna in exact_match_dict:
            return True, observed_grna, exact_match_dict[observed_grna]
        
        # Try 1-mismatch lookup if enabled (still O(1))
        if self.allow_grna_mismatch and observed_grna in one_mismatch_dict:
            original_grna, gene = one_mismatch_dict[observed_grna]
            return True, original_grna, gene
        
        # No match found
        return False, None, None
    
    def process(self, read1_file, read2_file, library_file, output_prefix):
        """Process FASTQ files with fast hash table lookup"""
        print("\n" + "=" * 60)
        print("PROCESSING FASTQ FILES (HASH TABLE MODE)")
        print("=" * 60)
        
        # Load library as hash tables
        exact_match_dict, one_mismatch_dict = self._load_library_as_dict(library_file)
        
        # Initialize counters
        counts = {celltype: defaultdict(int) for celltype in self.celltype_barcodes.keys()}
        stats = {
            'total': 0,
            'assigned': defaultdict(int),
            'unassigned': 0,
            'no_barcode': 0,
            'no_grna': 0,
            'invalid_grna': 0,
            'exact_match': 0,
            'one_mismatch': 0
        }
        
        # Open files
        r1_opener = gzip.open if read1_file.endswith('.gz') else open
        r2_opener = gzip.open if read2_file.endswith('.gz') else open
        
        print("\nProcessing reads...")
        start_time = time.time()
        
        with r1_opener(read1_file, 'rt') as r1, r2_opener(read2_file, 'rt') as r2:
            while True:
                # Read FASTQ records (4 lines each)
                r1_lines = [r1.readline().strip() for _ in range(4)]
                r2_lines = [r2.readline().strip() for _ in range(4)]
                
                if not all(r1_lines) or not all(r2_lines):
                    break
                
                stats['total'] += 1
                
                if stats['total'] % 100000 == 0:
                    elapsed = time.time() - start_time
                    rate = stats['total'] / elapsed if elapsed > 0 else 0
                    print(f"  Processed {stats['total']:,} reads... ({rate:,.0f} reads/sec)", end='\r')
                
                # Extract barcode from read2
                barcode = self.extract_barcode(r2_lines[1])
                if not barcode:
                    stats['no_barcode'] += 1
                    continue
                
                # Extract gRNA from read1
                grna = self.extract_grna(r1_lines[1])
                if not grna:
                    stats['no_grna'] += 1
                    continue
                
                # Fast hash table lookup for gRNA
                is_valid, matched_grna, gene = self.match_grna_fast(
                    grna, exact_match_dict, one_mismatch_dict
                )
                
                if not is_valid:
                    stats['invalid_grna'] += 1
                    continue
                
                # Track match type
                if matched_grna == grna:
                    stats['exact_match'] += 1
                else:
                    stats['one_mismatch'] += 1
                
                # Match barcode to cell type
                celltype, distance = self.match_barcode(barcode)
                
                # Track barcode statistics
                self.barcode_stats[barcode]['total'] += 1
                
                if celltype:
                    counts[celltype][matched_grna] += 1
                    stats['assigned'][celltype] += 1
                    self.celltype_counts[celltype] += 1
                    self.barcode_stats[barcode]['assigned'] += 1
                    self.barcode_stats[barcode]['celltype'] = celltype
                    if distance == 0:
                        self.barcode_stats[barcode]['exact_match'] += 1
                    else:
                        self.barcode_stats[barcode]['corrected'] += 1
                else:
                    stats['unassigned'] += 1
                    self.celltype_counts['unassigned'] += 1
                    self.barcode_stats[barcode]['unassigned'] += 1
        
        elapsed = time.time() - start_time
        print(f"\n  Processed {stats['total']:,} total reads in {elapsed:.2f} seconds")
        print(f"  Average speed: {stats['total']/elapsed:,.0f} reads/second")
        
        # Generate output
        self._write_outputs(counts, exact_match_dict, output_prefix, stats)
        
        return stats
    
    def _write_outputs(self, counts, library_dict, output_prefix, stats):
        """Write count table and summary"""
        # Main count table
        count_file = f"{output_prefix}_count.txt"
        
        with open(count_file, 'w') as f:
            celltypes = sorted(self.celltype_barcodes.keys())
            f.write("sgRNA\tGene\t" + "\t".join(celltypes) + "\n")
            
            for grna in sorted(library_dict.keys()):
                gene = library_dict[grna]
                celltype_counts = [str(counts[ct].get(grna, 0)) for ct in celltypes]
                f.write(f"{grna}\t{gene}\t" + "\t".join(celltype_counts) + "\n")
        
        print(f"\n✓ Count file: {count_file}")
        
        # Cell type proportions
        self._write_celltype_proportions(output_prefix, stats)
        
        # Barcode performance
        self._write_barcode_performance(output_prefix, stats)
        
        # Summary with speed metrics
        summary_file = f"{output_prefix}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("DEMULTIPLEXING SUMMARY (HASH TABLE MODE)\n")
            f.write("=" * 60 + "\n\n")
            f.write(f"Total reads: {stats['total']:,}\n")
            f.write(f"Assigned reads: {sum(stats['assigned'].values()):,}\n")
            f.write(f"Unassigned reads: {stats['unassigned']:,}\n")
            f.write(f"No barcode: {stats['no_barcode']:,}\n")
            f.write(f"No gRNA: {stats['no_grna']:,}\n")
            f.write(f"Invalid gRNA: {stats['invalid_grna']:,}\n\n")
            
            f.write("gRNA MATCHING STATISTICS:\n")
            f.write(f"  Exact matches: {stats['exact_match']:,} "
                   f"({stats['exact_match']/stats['total']*100:.1f}%)\n")
            if self.allow_grna_mismatch:
                f.write(f"  1-mismatch matches: {stats['one_mismatch']:,} "
                       f"({stats['one_mismatch']/stats['total']*100:.1f}%)\n")
            f.write(f"  Total valid: {stats['exact_match'] + stats['one_mismatch']:,}\n\n")
            
            f.write("CELL TYPE DISTRIBUTION:\n")
            for celltype, count in sorted(stats['assigned'].items()):
                pct = (count / stats['total'] * 100) if stats['total'] > 0 else 0
                f.write(f"  {celltype}: {count:,} ({pct:.1f}%)\n")
        
        print(f"✓ Summary: {summary_file}")
        
        # Print stats
        print("\nSTATS:")
        print(f"  Total reads: {stats['total']:,}")
        print(f"  Valid gRNA matches: {stats['exact_match'] + stats['one_mismatch']:,}")
        print(f"    Exact: {stats['exact_match']:,} ({stats['exact_match']/stats['total']*100:.1f}%)")
        if self.allow_grna_mismatch:
            print(f"    1-mismatch: {stats['one_mismatch']:,} ({stats['one_mismatch']/stats['total']*100:.1f}%)")
        print(f"  Assigned to cell types: {sum(stats['assigned'].values()):,}")
        for celltype, count in sorted(stats['assigned'].items()):
            print(f"    {celltype}: {count:,}")
        print(f"  Unassigned: {stats['unassigned']:,}")
        print(f"  Invalid gRNA: {stats['invalid_grna']:,}")
    
    def _write_celltype_proportions(self, output_prefix, stats):
        """Write cell type proportion report"""
        prop_file = f"{output_prefix}_celltype_proportions.txt"
        
        total_assigned = sum(stats['assigned'].values())
        
        with open(prop_file, 'w') as f:
            f.write("CellType\tRead_Count\tPercentage\tExpected_Barcode\n")
            
            for celltype in sorted(self.celltype_barcodes.keys()):
                count = stats['assigned'].get(celltype, 0)
                pct = (count / total_assigned * 100) if total_assigned > 0 else 0
                expected_bc = self.celltype_barcodes[celltype]
                
                f.write(f"{celltype}\t{count}\t{pct:.2f}\t{expected_bc}\n")
            
            # Unassigned
            unassigned_count = stats['unassigned']
            unassigned_pct = (unassigned_count / stats['total'] * 100) if stats['total'] > 0 else 0
            f.write(f"unassigned\t{unassigned_count}\t{unassigned_pct:.2f}\tN/A\n")
        
        print(f"✓ Cell type proportions: {prop_file}")
    
    def _write_barcode_performance(self, output_prefix, stats):
        """Write barcode performance report"""
        perf_file = f"{output_prefix}_barcode_performance.txt"
        
        with open(perf_file, 'w') as f:
            f.write("Barcode\tCellType\tTotal_Reads\tAssigned_Reads\t"
                   "Assignment_Rate\tExact_Matches\tCorrected\n")
            
            # Get top barcodes by frequency
            top_barcodes = sorted(
                self.barcode_stats.items(),
                key=lambda x: x[1]['total'],
                reverse=True
            )[:20]  # Top 20 barcodes
            
            for barcode, bc_stats in top_barcodes:
                total = bc_stats['total']
                assigned = bc_stats.get('assigned', 0)
                exact = bc_stats.get('exact_match', 0)
                corrected = bc_stats.get('corrected', 0)
                celltype = bc_stats.get('celltype', 'unassigned')
                
                assignment_rate = (assigned / total * 100) if total > 0 else 0
                
                f.write(f"{barcode}\t{celltype}\t{total}\t{assigned}\t"
                       f"{assignment_rate:.1f}%\t{exact}\t{corrected}\n")
        
        print(f"✓ Barcode performance: {perf_file}")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Fast MAGeCK Cell Type Demultiplexer (Hash Table Implementation)'
    )
    parser.add_argument('--read1', required=True, help='Read 1 FASTQ')
    parser.add_argument('--read2', required=True, help='Read 2 FASTQ')
    parser.add_argument('--library', required=True, help='gRNA library')
    parser.add_argument('--celltype-barcodes', required=True, help='Barcode CSV')
    parser.add_argument('--output-prefix', required=True, help='Output prefix')
    parser.add_argument('--barcode-start', type=int, default=22, help='Barcode start position')
    parser.add_argument('--barcode-length', type=int, default=8, help='Barcode length')
    parser.add_argument('--grna-start', type=int, default=12, help='gRNA start position')
    parser.add_argument('--grna-length', type=int, default=20, help='gRNA length')
    parser.add_argument('--max-barcode-mismatches', type=int, default=1, 
                       help='Max barcode mismatches')
    parser.add_argument('--no-grna-mismatch', action='store_true',
                       help='Disable 1-mismatch tolerance for gRNA (use exact match only)')
    
    args = parser.parse_args()
    
    # Validate inputs
    for f in [args.read1, args.read2, args.library, args.celltype_barcodes]:
        if not os.path.exists(f):
            print(f"ERROR: File not found: {f}")
            return
    
    # Create output directory
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    # Run demultiplexing
    demux = FastDemultiplexer(
        barcode_csv=args.celltype_barcodes,
        barcode_start=args.barcode_start,
        barcode_length=args.barcode_length,
        grna_start=args.grna_start,
        grna_length=args.grna_length,
        max_barcode_mismatches=args.max_barcode_mismatches,
        allow_grna_mismatch=not args.no_grna_mismatch
    )
    
    demux.process(
        read1_file=args.read1,
        read2_file=args.read2,
        library_file=args.library,
        output_prefix=args.output_prefix
    )
    
    print("\n" + "=" * 60)
    print("DEMULTIPLEXING COMPLETE!")
    print("=" * 60)


if __name__ == "__main__":
    main()
