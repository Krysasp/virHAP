#!/usr/bin/env python3
"""
Generate protein-level VCF and mutation summary files.
WORKFLOW:
1. Align each sample to reference (BOTH can have gaps)
2. Track coordinate mapping: each position in aligned sample maps to reference position
3. Extract CDS regions using ORIGINAL reference coordinates, but follow the alignment
4. Handle insertions (gaps in reference) and deletions (gaps in sample) properly
5. Translate extracted CDS regions to amino acids, preserving gap information
6. Compare amino acid sequences to call protein variants
7. Report ONE entry per amino acid change
"""

import argparse
import csv
from collections import defaultdict
from datetime import datetime
import os
import re
import sys
import time
import threading
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio.Align import PairwiseAligner


class Timer:
    """Simple timer for measuring duration of operations."""
    
    def __init__(self, name: str, verbose: bool = True):
        self.name = name
        self.verbose = verbose
        self.start_time = None
        self.end_time = None
        self.elapsed = None
    
    def __enter__(self):
        self.start_time = time.time()
        if self.verbose:
            print(f"\n[Timer] Starting: {self.name}")
        return self
    
    def __exit__(self, *args):
        self.end_time = time.time()
        self.elapsed = self.end_time - self.start_time
        if self.verbose:
            print(f"[Timer] Completed: {self.name} - Duration: {self.elapsed:.3f} seconds")


class ProgressTracker:
    """Track and display progress of operations with thread safety."""
    
    def __init__(self, total_items: int, description: str = "Processing", verbose: bool = True):
        self.total = total_items
        self.current = 0
        self.description = description
        self.verbose = verbose
        self.start_time = time.time()
        self.last_percent = 0
        self.lock = threading.Lock()
        
        if self.verbose and self.total > 0:
            print(f"\n{self.description}: 0% (0/{self.total})", end='', flush=True)
    
    def update(self, increment: int = 1):
        with self.lock:
            self.current += increment
            
            if self.verbose and self.total > 0:
                percent = int((self.current / self.total) * 100)
                if percent >= self.last_percent + 5 or self.current == self.total:
                    elapsed = time.time() - self.start_time
                    if self.current > 0:
                        rate = self.current / elapsed
                        remaining = (self.total - self.current) / rate if rate > 0 else 0
                        time_str = f", {elapsed:.1f}s elapsed, {remaining:.1f}s remaining"
                    else:
                        time_str = ""
                    
                    print(f"\r{self.description}: {percent}% ({self.current}/{self.total}){time_str}", 
                          end='', flush=True)
                    self.last_percent = percent
    
    def complete(self):
        if self.verbose and self.total > 0:
            elapsed = time.time() - self.start_time
            print(f"\r{self.description}: 100% ({self.total}/{self.total}) - {elapsed:.1f}s total")
        elif self.verbose:
            print(f"\n{self.description}: Complete")


@dataclass
class CDSFeature:
    """Store CDS feature information from GenBank."""
    name: str
    start: int
    end: int
    strand: int
    product: str
    locus_tag: str
    translation: str = ""
    
    @property
    def length(self) -> int:
        return self.end - self.start + 1
    
    @property
    def codon_length(self) -> int:
        return self.length // 3


@dataclass
class ProteinVariant:
    """
    Represents a SINGLE protein-level change.
    This is the atomic unit of variation - one per amino acid change.
    """
    cds_name: str
    codon_idx: int
    ref_aa: str
    alt_aa: str
    samples: Set[str]
    
    @property
    def protein_change(self) -> str:
        """Get protein change string (e.g., 'G:M163T')."""
        return f"{self.cds_name}:{self.ref_aa}{self.codon_idx}{self.alt_aa}"


class GenBankParser:
    """Parse GenBank file to extract CDS features with single-letter protein names."""
    
    PROTEIN_MAP = {
        'NUCLEOPROTEIN': 'N',
        'PHOSPHOPROTEIN': 'P',
        'MATRIX': 'M',
        'GLYCOPROTEIN': 'G',
        'LARGE': 'L',
        'POLYMERASE': 'L'
    }
    
    def __init__(self, genbank_file: str, verbose: bool = True):
        self.genbank_file = genbank_file
        self.cds_features: Dict[str, CDSFeature] = {}
        self.sequence = None
        self.verbose = verbose
        with Timer("GenBank parsing", verbose):
            self._parse()
    
    def _extract_single_letter_name(self, product: str, locus_tag: str = "") -> str:
        product = product.strip('"').strip().upper()
        
        if self.verbose:
            print(f"    Extracting name from product: '{product}'")
        
        direct_match = re.match(r'^([A-Z0-9]+)$', product)
        if direct_match:
            name = direct_match.group(1)
            if self.verbose:
                print(f"      → Found name: '{name}' (direct match)")
            return name
        
        for protein_name, letter in self.PROTEIN_MAP.items():
            pattern = rf'{protein_name}\s+([A-Z0-9]+)'
            match = re.search(pattern, product)
            if match:
                specific = match.group(1)
                if specific == letter or len(specific) == 1:
                    name = specific
                else:
                    name = f"{letter}{specific}"
                if self.verbose:
                    print(f"      → Found name: '{name}' (protein + specific)")
                return name
        
        for protein_name, letter in self.PROTEIN_MAP.items():
            if protein_name in product:
                if self.verbose:
                    print(f"      → Found name: '{letter}' (protein name only)")
                return letter
        
        match = re.search(r'([NPGML][0-9]*)$', product)
        if match:
            name = match.group(1)
            if self.verbose:
                print(f"      → Found name: '{name}' (letter at end)")
            return name
        
        if locus_tag:
            num_match = re.search(r'gp(\d+)', locus_tag.lower())
            if num_match:
                num = num_match.group(1)
                name = f"M{num}"
                if self.verbose:
                    print(f"      → Using locus tag derived name: '{name}'")
                return name
        
        name = product[:2] if len(product) > 1 else product[0]
        if self.verbose:
            print(f"      → Using fallback name: '{name}'")
        return name
    
    def _parse(self):
        print(f"\n[GenBank Parser] Reading file: {self.genbank_file}")
        
        with open(self.genbank_file, 'r') as f:
            records = list(SeqIO.parse(f, 'genbank'))
            if not records:
                raise ValueError(f"No valid GenBank records found in {self.genbank_file}")
            
            record = records[0]
            self.sequence = str(record.seq).upper()
            print(f"  Loaded sequence: {len(self.sequence)} bp")
            
            cds_count = sum(1 for feature in record.features if feature.type == 'CDS')
            print(f"  Found {cds_count} CDS features")
            
            progress = ProgressTracker(cds_count, "  Parsing CDS features", self.verbose)
            
            for feature in record.features:
                if feature.type == 'CDS':
                    start = int(feature.location.start) + 1
                    end = int(feature.location.end)
                    strand = feature.location.strand
                    
                    product = feature.qualifiers.get('product', [''])[0]
                    locus_tag = feature.qualifiers.get('locus_tag', [''])[0]
                    translation = feature.qualifiers.get('translation', [''])[0]
                    
                    if self.verbose:
                        print(f"\n    Processing CDS: {locus_tag}")
                        print(f"      Coordinates: {start}-{end} (strand={strand})")
                        print(f"      Product: '{product}'")
                    
                    protein_name = self._extract_single_letter_name(product, locus_tag)
                    
                    original_name = protein_name
                    suffix = 1
                    while protein_name in self.cds_features:
                        protein_name = f"{original_name}_{suffix}"
                        suffix += 1
                    
                    cds = CDSFeature(
                        name=protein_name,
                        start=start,
                        end=end,
                        strand=strand,
                        product=product,
                        locus_tag=locus_tag,
                        translation=translation
                    )
                    
                    self.cds_features[protein_name] = cds
                    
                    if self.verbose:
                        print(f"      → Registered as: '{protein_name}'")
                        print(f"      → Length: {cds.length} nt ({cds.codon_length} codons)")
                    
                    progress.update()
            
            progress.complete()
            
            print(f"\n[GenBank Parser] Extracted {len(self.cds_features)} CDS features:")
            for name, cds in sorted(self.cds_features.items()):
                print(f"  {name}: {cds.start}-{cds.end} (strand={cds.strand}, {cds.codon_length} codons)")


class SequenceAligner:
    """Handle sequence alignment - both sequences can have gaps."""
    
    def __init__(self):
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -5
        self.aligner.extend_gap_score = -2
    
    def align_pairwise(self, seq1: str, seq2: str, name1: str = "seq1", name2: str = "seq2") -> Tuple[str, str, List[Tuple[int, int]], float]:
        """
        Perform pairwise alignment between two sequences.
        Both sequences can have gaps in the alignment.
        Returns aligned sequences and a detailed mapping between positions.
        """
        # Remove any existing gaps for proper alignment
        seq1_clean = seq1.replace('-', '')
        seq2_clean = seq2.replace('-', '')
        
        # Perform alignment
        alignments = self.aligner.align(seq1_clean, seq2_clean)
        if not alignments:
            print(f"      Warning: No alignment found between {name1} and {name2}, using simple alignment")
            return self._simple_align_pair(seq1, seq2)
        
        # Get the best alignment
        alignment = alignments[0]
        aligned_seq1 = alignment[0]
        aligned_seq2 = alignment[1]
        score = alignment.score
        
        # Create detailed mapping between positions
        # Each entry is (ref_pos, sample_pos) where:
        # - ref_pos: position in reference sequence (1-based, 0 for gaps)
        # - sample_pos: position in sample sequence (1-based, 0 for gaps)
        pos_map = []
        ref_idx = 0
        sample_idx = 0
        
        for ref_char, sample_char in zip(aligned_seq1, aligned_seq2):
            ref_pos = ref_idx + 1 if ref_char != '-' else 0
            sample_pos = sample_idx + 1 if sample_char != '-' else 0
            
            pos_map.append((ref_pos, sample_pos))
            
            if ref_char != '-':
                ref_idx += 1
            if sample_char != '-':
                sample_idx += 1
        
        return aligned_seq1, aligned_seq2, pos_map, score
    
    def _simple_align_pair(self, seq1: str, seq2: str) -> Tuple[str, str, List[Tuple[int, int]], float]:
        """Simple alignment fallback."""
        if len(seq1) == len(seq2):
            pos_map = [(i+1, i+1) for i in range(len(seq1))]
            return seq1, seq2, pos_map, 0.0
        
        # Pad the shorter sequence with Ns
        max_len = max(len(seq1), len(seq2))
        seq1_padded = seq1.ljust(max_len, 'N')
        seq2_padded = seq2.ljust(max_len, 'N')
        pos_map = [(i+1, i+1) for i in range(max_len)]
        
        return seq1_padded, seq2_padded, pos_map, 0.0


class ProteinTranslator:
    """Translate nucleotide sequences to amino acids, handling gaps appropriately."""
    
    def __init__(self):
        self.codon_table = {
            'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
            'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
            'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
            'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
            'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
            'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
            'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
            'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
            'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
            'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
            'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
            'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
            'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
            'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
            'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
            'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
        }
    
    def extract_cds_with_mapping(self, aligned_ref: str, aligned_sample: str, 
                                 pos_map: List[Tuple[int, int]], 
                                 cds: CDSFeature, sample_name: str = "sample") -> Tuple[str, str, bool, Dict[str, Any]]:
        """
        Extract CDS region using detailed position mapping.
        This function properly handles insertions and deletions.
        
        Args:
            aligned_ref: Reference sequence after alignment (may have gaps)
            aligned_sample: Sample sequence after alignment (may have gaps)
            pos_map: List of (ref_pos, sample_pos) tuples for each aligned column
            cds: CDS feature with original coordinates
            sample_name: Name of sample for logging
        
        Returns:
            ref_cds: Reference CDS nucleotides (gaps preserved as '-')
            sample_cds: Sample CDS nucleotides (gaps preserved as '-')
            success: Whether extraction was successful
            stats: Dictionary with extraction statistics
        """
        stats = {
            'total_ref_positions': cds.end - cds.start + 1,
            'ref_positions_found': 0,
            'sample_positions_found': 0,
            'insertions': 0,
            'deletions': 0,
            'matches': 0,
            'mismatches': 0
        }
        
        # Build mapping from reference position to aligned columns
        ref_to_columns = defaultdict(list)
        sample_to_columns = defaultdict(list)
        
        for col_idx, (ref_pos, sample_pos) in enumerate(pos_map):
            if ref_pos > 0:
                ref_to_columns[ref_pos].append(col_idx)
            if sample_pos > 0:
                sample_to_columns[sample_pos].append(col_idx)
        
        # Extract CDS region
        ref_cds_bases = []
        sample_cds_bases = []
        
        # Track which reference positions we've processed
        ref_positions_processed = set()
        
        # Process each reference position in the CDS
        for ref_pos in range(cds.start, cds.end + 1):
            if ref_pos in ref_to_columns:
                # This reference position maps to one or more columns
                columns = ref_to_columns[ref_pos]
                
                for col_idx in columns:
                    ref_char = aligned_ref[col_idx]
                    sample_char = aligned_sample[col_idx]
                    _, sample_pos = pos_map[col_idx]
                    
                    ref_cds_bases.append(ref_char)
                    sample_cds_bases.append(sample_char)
                    
                    ref_positions_processed.add(ref_pos)
                    stats['ref_positions_found'] += 1
                    
                    if sample_pos > 0:
                        stats['sample_positions_found'] += 1
                    
                    # Classify the event
                    if ref_char != '-' and sample_char != '-':
                        if ref_char == sample_char:
                            stats['matches'] += 1
                        else:
                            stats['mismatches'] += 1
                    elif ref_char == '-' and sample_char != '-':
                        stats['insertions'] += 1
                    elif ref_char != '-' and sample_char == '-':
                        stats['deletions'] += 1
        
        # Check if we found all needed positions
        missing_positions = set(range(cds.start, cds.end + 1)) - ref_positions_processed
        if missing_positions:
            print(f"      Warning: Sample {sample_name} missing CDS positions: {sorted(missing_positions)}")
            # Fill missing positions with appropriate gaps
            for ref_pos in sorted(missing_positions):
                # Add a gap in both sequences for missing positions
                ref_cds_bases.append(self.ref_seq[ref_pos-1] if hasattr(self, 'ref_seq') else 'N')
                sample_cds_bases.append('-')
                stats['deletions'] += 1
        
        # Convert to strings
        ref_cds = ''.join(ref_cds_bases)
        sample_cds = ''.join(sample_cds_bases)
        
        # Handle reverse strand
        if cds.strand == -1:
            ref_cds = self._reverse_complement_with_gaps(ref_cds)
            sample_cds = self._reverse_complement_with_gaps(sample_cds)
        
        return ref_cds, sample_cds, True, stats
    
    def _reverse_complement_with_gaps(self, seq: str) -> str:
        """Reverse complement a sequence, preserving gaps as '-'."""
        comp_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N', '-': '-'}
        complemented = ''.join(comp_map.get(base.upper(), 'N') for base in seq)
        return complemented[::-1]
    
    def translate_to_aa(self, nucleotide_seq: str, preserve_gaps: bool = True) -> str:
        """
        Translate a nucleotide sequence to amino acids.
        
        Args:
            nucleotide_seq: Nucleotide sequence (may contain gaps)
            preserve_gaps: If True, gaps in codons are represented as '-' in output
        
        Returns:
            Amino acid sequence with '-' for gaps and 'X' for unknown
        """
        amino_acids = []
        
        i = 0
        while i < len(nucleotide_seq):
            # Get next codon (may contain gaps)
            codon = nucleotide_seq[i:i+3]
            
            if len(codon) < 3:
                # Incomplete codon at end - pad with gaps
                codon = codon.ljust(3, '-')
            
            # Check if codon contains only gaps
            if codon == '---':
                amino_acids.append('-')
            # Check if codon has any gaps
            elif '-' in codon:
                if preserve_gaps:
                    # For mixed gaps, represent as unknown
                    amino_acids.append('X')
                else:
                    # Remove gaps and try to translate what's left
                    clean_codon = codon.replace('-', '')
                    if len(clean_codon) == 3 and 'N' not in clean_codon:
                        amino_acids.append(self.codon_table.get(clean_codon, 'X'))
                    else:
                        amino_acids.append('X')
            else:
                # No gaps, translate normally
                if 'N' not in codon:
                    amino_acids.append(self.codon_table.get(codon, 'X'))
                else:
                    amino_acids.append('X')
            
            i += 3
        
        return ''.join(amino_acids)


class ProteinVariantCaller:
    """
    Call variants at the protein level.
    WORKFLOW:
    1. Align each sample to reference (BOTH can have gaps)
    2. Create detailed position mapping
    3. Extract CDS regions using position mapping (handling indels)
    4. Translate to amino acids
    5. Compare amino acid sequences to call variants
    """
    
    def __init__(self, genbank_parser: GenBankParser, reference_file: str, sample_file: str, 
                 verbose: bool = True, num_threads: int = None, output_aligned: str = None,
                 output_aa: str = None):
        self.gb = genbank_parser
        self.verbose = verbose
        self.num_threads = num_threads or multiprocessing.cpu_count()
        self.output_aligned = output_aligned
        self.output_aa = output_aa
        self.aligner = SequenceAligner()
        self.translator = ProteinTranslator()
        
        print(f"\n[ProteinVariantCaller] Initializing with {self.num_threads} threads")
        print(f"  WORKFLOW:")
        print(f"    1. Align samples to reference (BOTH can have gaps)")
        print(f"    2. Create detailed position mapping")
        print(f"    3. Extract CDS regions using position mapping (handling indels)")
        print(f"    4. Translate to amino acids")
        print(f"    5. Call protein variants")
        
        with Timer("Loading reference and sample sequences", verbose):
            # Load reference
            self.ref_record = SeqIO.read(reference_file, 'fasta')
            self.ref_seq = str(self.ref_record.seq).upper()
            self.ref_name = self.ref_record.id
            print(f"  Reference: {self.ref_name}, length: {len(self.ref_seq)} bp")
            
            # Store reference sequence in translator for gap filling
            self.translator.ref_seq = self.ref_seq
            
            # Load all samples
            self.samples = list(SeqIO.parse(sample_file, 'fasta'))
            print(f"  Loaded {len(self.samples)} samples")
        
        with Timer("Extracting reference CDS sequences", verbose):
            # Store reference CDS sequences
            print(f"  Extracting reference CDS regions...")
            self.ref_cds_seqs = {}
            self.ref_aa_seqs = {}
            
            for cds_name, cds in self.gb.cds_features.items():
                # For reference, we can extract directly (no alignment needed)
                ref_cds = self.ref_seq[cds.start-1:cds.end]
                
                # Handle reverse strand
                if cds.strand == -1:
                    ref_cds = self.translator._reverse_complement_with_gaps(ref_cds)
                
                # Store nucleotide sequence (no gaps in reference)
                self.ref_cds_seqs[cds_name] = ref_cds
                
                # Translate to amino acids
                ref_aa = self.translator.translate_to_aa(ref_cds, preserve_gaps=False)
                self.ref_aa_seqs[cds_name] = ref_aa
                
                if self.verbose:
                    print(f"    {cds_name}: {len(ref_aa)} aa ({len(ref_cds)} nt)")
        
        # Store data for output
        self.aligned_pairs = []  # List of (ref_aligned, sample_aligned, sample_name, score, pos_map)
        self.aa_sequences = defaultdict(list)  # cds_name -> list of (sample_name, aa_seq)
        self.cds_nt_sequences = defaultdict(list)  # cds_name -> list of (sample_name, nt_seq_with_gaps)
        self.extraction_stats = defaultdict(list)  # sample_name -> list of (cds_name, stats)
        
        # Add reference to AA sequences (for output)
        for cds_name, ref_aa in self.ref_aa_seqs.items():
            self.aa_sequences[cds_name].append(("REFERENCE", ref_aa))
            self.cds_nt_sequences[cds_name].append(("REFERENCE", self.ref_cds_seqs[cds_name]))
    
    def process_single_sample(self, sample: SeqRecord) -> Tuple[str, Dict[str, str], Dict[str, str], str, str, List[Tuple[int, int]], float]:
        """
        Process a single sample with proper handling of indels.
        """
        sample_name = sample.id
        sample_seq = str(sample.seq).upper()
        
        if self.verbose and sample_name == self.samples[0].id:
            print(f"\n    Processing sample: {sample_name}")
            print(f"      Original length: {len(sample_seq)} bp")
            if '-' in sample_seq:
                gap_count = sample_seq.count('-')
                print(f"      Contains {gap_count} gaps in input sequence")
        
        # STEP 1: Align sample to reference
        aligned_ref, aligned_sample, pos_map, score = self.aligner.align_pairwise(
            self.ref_seq, sample_seq, "reference", sample_name
        )
        
        if self.verbose and sample_name == self.samples[0].id:
            print(f"      Aligned lengths: ref={len(aligned_ref)} bp, sample={len(aligned_sample)} bp")
            print(f"      Alignment score: {score:.2f}")
            
            # Count gaps in each sequence
            ref_gaps = aligned_ref.count('-')
            sample_gaps = aligned_sample.count('-')
            if ref_gaps > 0 or sample_gaps > 0:
                print(f"      Gaps after alignment - ref: {ref_gaps}, sample: {sample_gaps}")
        
        # STEP 2-3: Extract CDS regions using detailed mapping
        sample_proteins = {}
        sample_cds_nt = {}
        failed_cds = []
        
        for cds_name, cds in self.gb.cds_features.items():
            # Extract CDS using detailed position mapping
            ref_cds, sample_cds, success, stats = self.translator.extract_cds_with_mapping(
                aligned_ref, aligned_sample, pos_map, cds, sample_name
            )
            
            # Store extraction statistics
            self.extraction_stats[sample_name].append((cds_name, stats))
            
            if not success:
                failed_cds.append(cds_name)
                # Create a CDS sequence with appropriate gaps for failed extraction
                expected_length = cds.length
                sample_cds = '-' * expected_length
                ref_cds = self.ref_cds_seqs[cds_name]
                
                if self.verbose:
                    print(f"      Warning: Failed to extract {cds_name} for {sample_name}, using all gaps")
            
            # Store nucleotide sequence with gaps preserved
            sample_cds_nt[cds_name] = sample_cds
            
            # Translate to amino acids
            sample_aa = self.translator.translate_to_aa(sample_cds, preserve_gaps=True)
            sample_proteins[cds_name] = sample_aa
            
            # Store for output
            self.aa_sequences[cds_name].append((sample_name, sample_aa))
            self.cds_nt_sequences[cds_name].append((sample_name, sample_cds))
        
        if failed_cds and self.verbose:
            print(f"      Warning: Failed to extract CDS regions: {', '.join(failed_cds)}")
        
        return sample_name, sample_proteins, sample_cds_nt, aligned_ref, aligned_sample, pos_map, score
    
    def call_protein_variants(self, parallel: bool = True) -> Tuple[Dict[str, List[ProteinVariant]], Dict[str, Dict[str, str]]]:
        """
        Main method to call protein-level variants.
        """
        print(f"\n[ProteinVariantCaller] Processing samples...")
        
        # Store protein sequences for each sample
        sample_proteins = {}
        sample_cds_nt = {}
        alignment_scores = {}
        self.aligned_pairs = []
        
        with Timer(f"Aligning and translating {len(self.samples)} samples", self.verbose):
            if parallel and self.num_threads > 1:
                # Parallel processing
                with ThreadPoolExecutor(max_workers=self.num_threads) as executor:
                    future_to_sample = {
                        executor.submit(self.process_single_sample, sample): sample
                        for sample in self.samples
                    }
                    
                    progress = ProgressTracker(len(self.samples), "  Aligning and translating", self.verbose)
                    
                    for future in as_completed(future_to_sample):
                        sample = future_to_sample[future]
                        sample_name, proteins, cds_nt, aligned_ref, aligned_sample, pos_map, score = future.result()
                        sample_proteins[sample_name] = proteins
                        sample_cds_nt[sample_name] = cds_nt
                        alignment_scores[sample_name] = score
                        
                        self.aligned_pairs.append((aligned_ref, aligned_sample, sample_name, score, pos_map))
                        
                        progress.update()
                    
                    progress.complete()
            else:
                # Sequential processing
                progress = ProgressTracker(len(self.samples), "  Aligning and translating", self.verbose)
                
                for sample in self.samples:
                    sample_name, proteins, cds_nt, aligned_ref, aligned_sample, pos_map, score = self.process_single_sample(sample)
                    sample_proteins[sample_name] = proteins
                    sample_cds_nt[sample_name] = cds_nt
                    alignment_scores[sample_name] = score
                    
                    self.aligned_pairs.append((aligned_ref, aligned_sample, sample_name, score, pos_map))
                    
                    progress.update()
                
                progress.complete()
        
        # Report alignment and extraction statistics
        if alignment_scores:
            avg_score = sum(alignment_scores.values()) / len(alignment_scores)
            print(f"\n  Alignment statistics:")
            print(f"    Average alignment score: {avg_score:.2f}")
            print(f"    Min score: {min(alignment_scores.values()):.2f}")
            print(f"    Max score: {max(alignment_scores.values()):.2f}")
        
        # Summarize extraction statistics
        if self.verbose and self.extraction_stats:
            print(f"\n  CDS extraction statistics:")
            cds_totals = defaultdict(lambda: {'insertions': 0, 'deletions': 0, 'mismatches': 0})
            for sample_name, stats_list in self.extraction_stats.items():
                for cds_name, stats in stats_list:
                    cds_totals[cds_name]['insertions'] += stats['insertions']
                    cds_totals[cds_name]['deletions'] += stats['deletions']
                    cds_totals[cds_name]['mismatches'] += stats['mismatches']
            
            for cds_name, totals in sorted(cds_totals.items()):
                if totals['insertions'] > 0 or totals['deletions'] > 0 or totals['mismatches'] > 0:
                    print(f"    {cds_name}: ins={totals['insertions']}, del={totals['deletions']}, mis={totals['mismatches']}")
        
        # Write aligned sequences if requested
        if self.output_aligned:
            with Timer(f"Writing aligned sequences to {self.output_aligned}", self.verbose):
                self.write_aligned_sequences(self.output_aligned)
        
        # Write amino acid sequences if requested
        if self.output_aa:
            with Timer(f"Writing amino acid sequences to {self.output_aa}", self.verbose):
                self.write_aa_sequences(self.output_aa)
        
        # STEP 5: Call variants by comparing amino acid sequences
        print(f"\n  Calling protein-level variants...")
        
        with Timer("Variant calling", self.verbose):
            # Group variants by (cds_name, position, alt_aa)
            variant_dict = {}
            
            for cds_name, ref_aa_seq in self.ref_aa_seqs.items():
                cds_variant_count = 0
                for sample_name, proteins in sample_proteins.items():
                    if cds_name not in proteins:
                        continue
                    
                    sample_aa_seq = proteins[cds_name]
                    
                    # Compare amino acid by amino acid
                    min_len = min(len(ref_aa_seq), len(sample_aa_seq))
                    for pos in range(min_len):
                        ref_aa = ref_aa_seq[pos]
                        sample_aa = sample_aa_seq[pos]
                        
                        # Skip comparison if reference is gap (shouldn't happen) or sample is gap
                        if ref_aa == '-' or sample_aa == '-':
                            continue
                        
                        if ref_aa != sample_aa and sample_aa != 'X':
                            # This is a protein-level variant
                            key = (cds_name, pos + 1, ref_aa, sample_aa)
                            
                            if key not in variant_dict:
                                variant_dict[key] = set()
                            
                            variant_dict[key].add(sample_name)
                            cds_variant_count += 1
                
                if cds_variant_count > 0 and self.verbose:
                    print(f"    {cds_name}: {cds_variant_count} variant positions")
            
            # Create ProteinVariant objects
            protein_variants = []
            for (cds_name, pos, ref_aa, alt_aa), samples in variant_dict.items():
                variant = ProteinVariant(
                    cds_name=cds_name,
                    codon_idx=pos,
                    ref_aa=ref_aa,
                    alt_aa=alt_aa,
                    samples=samples
                )
                protein_variants.append(variant)
            
            # Sort by CDS and position
            protein_variants.sort(key=lambda v: (v.cds_name, v.codon_idx))
        
        print(f"  Found {len(protein_variants)} unique protein-level changes")
        
        # Create per-sample variant lists for summary
        sample_variants = defaultdict(list)
        for variant in protein_variants:
            for sample in variant.samples:
                sample_variants[sample].append(variant)
        
        return sample_variants, sample_proteins
    
    def write_aligned_sequences(self, output_file: str):
        """Write aligned nucleotide sequences to FASTA file."""
        if not self.aligned_pairs:
            print(f"  No aligned sequences to write")
            return
        
        print(f"\n[Aligned Nucleotide Sequences] Writing to: {output_file}")
        
        with open(output_file, 'w') as f:
            # Write reference first
            if self.aligned_pairs:
                first_ref, _, _, first_score, _ = self.aligned_pairs[0]
                ref_record = SeqRecord(
                    Seq(first_ref),
                    id=self.ref_name,
                    description=f"Reference sequence (aligned)"
                )
                SeqIO.write(ref_record, f, 'fasta')
            
            # Write each sample with its alignment
            for aligned_ref, aligned_sample, sample_name, score, pos_map in self.aligned_pairs:
                sample_record = SeqRecord(
                    Seq(aligned_sample),
                    id=sample_name,
                    description=f"Aligned to reference (score: {score:.2f})"
                )
                SeqIO.write(sample_record, f, 'fasta')
        
        print(f"  Wrote 1 reference + {len(self.aligned_pairs)} sample sequences")
    
    def write_aa_sequences(self, output_prefix: str):
        """Write amino acid and nucleotide sequences for each CDS."""
        if not self.aa_sequences:
            print(f"  No amino acid sequences to write")
            return
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_prefix)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"  Created output directory: {output_dir}")
        
        print(f"\n[Amino Acid Sequences] Writing to {output_prefix}_*.fasta")
        
        # Write amino acid sequences
        for cds_name, sequences in self.aa_sequences.items():
            output_file = f"{output_prefix}_{cds_name}_aa.fasta"
            records = []
            
            for name, aa_seq in sequences:
                gap_count = aa_seq.count('-')
                description = f"{name} {cds_name} protein"
                if gap_count > 0:
                    description += f" (contains {gap_count} gaps)"
                
                record = SeqRecord(
                    Seq(aa_seq),
                    id=name,
                    description=description
                )
                records.append(record)
            
            with open(output_file, 'w') as f:
                SeqIO.write(records, f, 'fasta')
            print(f"  Wrote {len(records)} sequences to {output_file}")
        
        # Write nucleotide CDS sequences
        print(f"\n[Nucleotide CDS Sequences] Writing to {output_prefix}_*_nt.fasta")
        for cds_name, sequences in self.cds_nt_sequences.items():
            output_file = f"{output_prefix}_{cds_name}_nt.fasta"
            records = []
            
            for name, nt_seq in sequences:
                gap_count = nt_seq.count('-')
                description = f"{name} {cds_name} nucleotide CDS"
                if gap_count > 0:
                    description += f" (contains {gap_count} gaps)"
                
                record = SeqRecord(
                    Seq(nt_seq),
                    id=name,
                    description=description
                )
                records.append(record)
            
            with open(output_file, 'w') as f:
                SeqIO.write(records, f, 'fasta')
            print(f"  Wrote {len(records)} sequences to {output_file}")


class VCFGenerator:
    """Generate VCF file - ONE ENTRY PER PROTEIN CHANGE."""
    
    @staticmethod
    def write_vcf(protein_variants: List[ProteinVariant],
                  output_file: str,
                  ref_name: str,
                  verbose: bool = True):
        """
        Write VCF with ONE ROW per protein change.
        """
        print(f"\n[VCFGenerator] Writing protein-level VCF to: {output_file}")
        
        with Timer(f"Writing VCF with {len(protein_variants)} entries", verbose):
            # Get all unique samples
            all_samples = set()
            for v in protein_variants:
                all_samples.update(v.samples)
            samples = sorted(all_samples)
            
            with open(output_file, 'w') as f:
                # Header
                f.write('##fileformat=VCFv4.2\n')
                f.write(f'##fileDate={datetime.now().strftime("%Y%m%d")}\n')
                f.write(f'##reference={ref_name}\n')
                f.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with this change">\n')
                f.write('##INFO=<ID=CDS,Number=1,Type=String,Description="CDS/protein name">\n')
                f.write('##INFO=<ID=AA_POS,Number=1,Type=Integer,Description="Amino acid position">\n')
                f.write('##INFO=<ID=REF_AA,Number=1,Type=String,Description="Reference amino acid">\n')
                f.write('##INFO=<ID=ALT_AA,Number=1,Type=String,Description="Alternate amino acid">\n')
                f.write('##INFO=<ID=PROTEIN_CHANGE,Number=1,Type=String,Description="Complete protein change">\n')
                f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (1/1 if present)">\n')
                f.write('##FORMAT=<ID=PROT,Number=1,Type=String,Description="Protein change">\n')
                
                # Column headers
                f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + 
                       '\t'.join(samples) + '\n')
                
                # Write variants
                progress = ProgressTracker(len(protein_variants), "  Writing VCF entries", verbose)
                
                for variant in protein_variants:
                    var_id = f"{variant.cds_name}_{variant.codon_idx}_{variant.ref_aa}{variant.alt_aa}"
                    
                    info = (f"NS={len(variant.samples)};CDS={variant.cds_name};AA_POS={variant.codon_idx};"
                           f"REF_AA={variant.ref_aa};ALT_AA={variant.alt_aa};"
                           f"PROTEIN_CHANGE={variant.protein_change}")
                    
                    # Protein-level VCF uses placeholder nucleotide values
                    pos = 1
                    ref_base = "N"
                    alt = "N"
                    
                    genotypes = []
                    for sample in samples:
                        if sample in variant.samples:
                            genotypes.append(f'1/1:{variant.protein_change}')
                        else:
                            genotypes.append('./.:.')
                    
                    f.write(f'{ref_name}\t{pos}\t{var_id}\t{ref_base}\t{alt}\t.\t.\t'
                           f'{info}\tGT:PROT\t' + '\t'.join(genotypes) + '\n')
                    
                    progress.update()
                
                progress.complete()
        
        print(f"  VCF written with {len(protein_variants)} entries")


class MutationSummary:
    """Generate mutation summary CSV files."""
    
    @staticmethod
    def write_aa_mutations(protein_variants: List[ProteinVariant],
                          samples: List[SeqRecord],
                          output_file: str,
                          verbose: bool = True):
        """
        Write amino acid mutations summary.
        """
        print(f"\n[MutationSummary] Writing protein-level mutations to: {output_file}")
        
        with Timer(f"Writing protein mutations summary", verbose):
            sample_to_variants = defaultdict(list)
            for variant in protein_variants:
                for sample in variant.samples:
                    sample_to_variants[sample].append(variant)
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Sample', 'Protein_Changes', 'Total_Changes'])
                
                progress = ProgressTracker(len(samples), "  Processing samples", verbose)
                
                for sample in samples:
                    sample_name = sample.id
                    variants = sample_to_variants.get(sample_name, [])
                    
                    def sort_key(v):
                        return (v.cds_name, v.codon_idx)
                    
                    variants.sort(key=sort_key)
                    
                    protein_changes = [v.protein_change for v in variants]
                    
                    writer.writerow([
                        sample_name,
                        ','.join(protein_changes) if protein_changes else '',
                        len(variants)
                    ])
                    
                    progress.update()
                
                progress.complete()
    
    @staticmethod
    def write_cds_summary(protein_variants: List[ProteinVariant],
                         gb: GenBankParser,
                         output_prefix: str,
                         verbose: bool = True):
        """
        Write CDS-level summary.
        """
        print(f"\n[MutationSummary] Writing CDS summary to: {output_prefix}_cds_summary.csv")
        
        with Timer(f"Writing CDS summary", verbose):
            cds_stats = defaultdict(lambda: {
                'total_changes': 0,
                'samples_affected': set(),
                'aa_positions': set(),
                'protein_changes': defaultdict(int)
            })
            
            for variant in protein_variants:
                stats = cds_stats[variant.cds_name]
                stats['total_changes'] += 1
                stats['samples_affected'].update(variant.samples)
                stats['aa_positions'].add(variant.codon_idx)
                stats['protein_changes'][variant.protein_change] = len(variant.samples)
            
            with open(f"{output_prefix}_cds_summary.csv", 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'Protein', 'Genomic_Coordinates', 'Length_aa', 'Total_Changes',
                    'Samples_Affected', 'Unique_Positions', 'Protein_Changes', 'Sample_Counts'
                ])
                
                progress = ProgressTracker(len(cds_stats), "  Writing CDS summaries", verbose)
                
                for cds_name in sorted(cds_stats.keys()):
                    cds = gb.cds_features[cds_name]
                    stats = cds_stats[cds_name]
                    
                    changes_str = ';'.join([
                        f"{chg}:{cnt}" for chg, cnt in sorted(stats['protein_changes'].items())
                    ])
                    
                    writer.writerow([
                        cds_name,
                        f"{cds.start}-{cds.end} (strand={cds.strand})",
                        cds.codon_length,
                        stats['total_changes'],
                        len(stats['samples_affected']),
                        len(stats['aa_positions']),
                        len(stats['protein_changes']),
                        changes_str
                    ])
                    
                    progress.update()
                
                progress.complete()


def main():
    parser = argparse.ArgumentParser(
        description='Generate protein-level VCF and mutation summaries.\n'
                   'Handles insertions and deletions properly.\n\n'
                   'EXAMPLE:\n'
                   '  python protein_variant_caller.py -g reference.gb -r reference.fasta -s samples.fasta -o output\n'
                   '  python protein_variant_caller.py -g reference.gb -r reference.fasta -s samples.fasta -o output --output-aa proteins/aa\n',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-g', '--genbank', required=True, help='GenBank annotation file')
    parser.add_argument('-r', '--reference', required=True, help='Reference FASTA file')
    parser.add_argument('-s', '--samples', required=True, help='Sample sequences FASTA file')
    parser.add_argument('-o', '--output', default='variants', help='Output prefix for VCF and summary files')
    parser.add_argument('--output-aligned', metavar='FILE', help='Output aligned nucleotide sequences')
    parser.add_argument('--output-aa', metavar='PREFIX', help='Output amino acid sequences (PREFIX_CDS_aa.fasta and PREFIX_CDS_nt.fasta)')
    parser.add_argument('--ref-name', default='reference', help='Reference name for VCF')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    parser.add_argument('-t', '--threads', type=int, default=None, help='Number of threads')
    parser.add_argument('--no-parallel', action='store_true', help='Disable parallel processing')
    
    args = parser.parse_args()
    verbose = not args.quiet
    
    print("\n" + "="*80)
    print("PROTEIN-LEVEL VARIANT CALLER (with indel support)")
    print("="*80)
    print("  WORKFLOW:")
    print("    1. Align samples to reference (BOTH can have gaps)")
    print("    2. Create detailed position mapping")
    print("    3. Extract CDS regions using mapping (handling indels)")
    print("    4. Translate to amino acids")
    print("    5. Call protein variants")
    print("-"*80)
    
    overall_timer = Timer("Complete pipeline", verbose)
    overall_timer.__enter__()
    
    # Parse GenBank
    gb = GenBankParser(args.genbank, verbose)
    
    # Initialize protein variant caller
    caller = ProteinVariantCaller(
        gb, args.reference, args.samples,
        verbose=verbose,
        num_threads=args.threads,
        output_aligned=args.output_aligned,
        output_aa=args.output_aa
    )
    
    # Call protein variants
    sample_variants, sample_proteins = caller.call_protein_variants(parallel=not args.no_parallel)
    
    # Get all unique protein variants
    all_variants = []
    seen = set()
    for variants in sample_variants.values():
        for v in variants:
            key = f"{v.cds_name}:{v.codon_idx}:{v.alt_aa}"
            if key not in seen:
                seen.add(key)
                all_variants.append(v)
    
    # Write VCF
    VCFGenerator.write_vcf(
        all_variants, f"{args.output}.vcf", args.ref_name,
        verbose=verbose
    )
    
    # Write protein mutations CSV
    samples = list(SeqIO.parse(args.samples, 'fasta'))
    MutationSummary.write_aa_mutations(
        all_variants, samples, f"{args.output}_protein_changes.csv",
        verbose=verbose
    )
    
    # Write CDS summary
    MutationSummary.write_cds_summary(all_variants, gb, args.output, verbose=verbose)
    
    overall_timer.__exit__()
    
    print("\n" + "="*80)
    print("PROCESSING COMPLETE")
    print("="*80)
    print(f"\nOutput files:")
    print(f"  VCF: {args.output}.vcf")
    print(f"  Protein changes: {args.output}_protein_changes.csv")
    print(f"  CDS Summary: {args.output}_cds_summary.csv")
    if args.output_aligned:
        print(f"  Aligned nucleotides: {args.output_aligned}")
    if args.output_aa:
        print(f"  Amino acid sequences: {args.output_aa}_*_aa.fasta")
        print(f"  Nucleotide CDS sequences: {args.output_aa}_*_nt.fasta")
    print("\nDone!")


if __name__ == '__main__':
    main()
