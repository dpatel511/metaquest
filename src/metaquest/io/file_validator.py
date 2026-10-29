"""
MetaQuest Enhanced File Validator v5.0.0 (FASTQ-only)
======================================================

Beautiful validation reports with spinner animations and table formatting.
Now exclusively validates FASTQ files.
"""

import os
from pathlib import Path
from Bio import SeqIO
import gzip
import hashlib
from itertools import zip_longest
import numpy as np
from collections import Counter
from datetime import datetime

class FileValidator:
    """Enhanced FASTQ file validation with beautiful output"""
    
    def __init__(self):
        self.min_sequences = 100
        self.quality_threshold = 20
        self.overrep_threshold = 0.1   # fraction (10%) — matches CLI --overrep-threshold default
        self.adapter_threshold = 5.0
        self.q30_threshold = 80.0
        self.sample_size = 10000
        self.MAX_READS_TO_ANALYZE = 1000_000
        self.strict = False
        self.common_adapters = [
            "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "CTGTCTCTTATACACATCT",
        ]
        
    def validate_and_analyze(self, input_file, file_type='fastq'):
        """
        Main validation with enhanced formatting.
        
        Args:
            input_file: Path to FASTQ file
            file_type: Must be 'fastq' (parameter kept for compatibility)
            
        Returns:
            tuple: (is_valid: bool, stats: dict)
        """
        from .output_formatter import get_formatter
        fmt = get_formatter()
        
        # Force FASTQ type
        if file_type != 'fastq':
            fmt.warning("Only FASTQ format is supported. Treating as FASTQ.")
            file_type = 'fastq'
        
        # Header
        _h = '\u2500' * 71
        _h73 = '\u2500' * 73
        fmt._print("\n\u250c" + _h + "\u2510")
        fmt._print("\u2502 \U0001f4cb FILE VALIDATION REPORT" + (' ' * 44) + "\u2502")
        fmt._print("\u251c" + _h + "\u2524")
        fmt._print(f"\u2502  File: {Path(input_file).name:<61}\u2502")
        fmt._print(f"\u2502  Path: {str(input_file)[:61]:<61}\u2502")
        fmt._print("\u2502  Type: FASTQ" + (' ' * 56) + "\u2502")
        fmt._print(f"\u2502  Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S'):<61}\u2502")
        fmt._print("\u2514" + _h + "\u2518\n")
        
        if not self._check_file_exists(input_file): 
            return False, None
        if not self._validate_format(input_file, 'fastq'): 
            return False, None
        
        # Analyze with spinner
        with fmt.spinner("Analyzing file contents"):
            stats = self._generate_statistics(input_file)
        
        self._display_statistics(stats)
        
        error_messages, warning_messages = self._get_validation_status(stats)
        has_errors = bool(error_messages)
        
        # Validation result box
        fmt._print("\n\u250c" + _h + "\u2510")
        if has_errors:
            fmt._print("\u2502 \u274c VALIDATION STATUS: FAILED" + (' ' * 42) + "\u2502")
            fmt._print("\u2514" + _h + "\u2518")
            fmt._print("\n\U0001f6d1 CRITICAL ISSUES DETECTED:")
            for msg in error_messages:
                fmt.error(msg)
            fmt._print("\n\U0001f4a1 RECOMMENDATIONS:")
            self._print_recommendations(stats, error_messages)
            fmt._print("\n" + _h73 + "\n")
            return False, stats
        else:
            for message in warning_messages:
                fmt.warning(message)
            fmt._print("\u2502 \u2705 VALIDATION STATUS: PASSED" + (' ' * 42) + "\u2502")
            fmt._print("\u2514" + _h + "\u2518")
            fmt._print("\n\u2713 All quality criteria met - file is ready for analysis")
            fmt._print(_h73 + "\n")
            return True, stats

    def _generate_statistics(self, input_file):
        """Generate comprehensive FASTQ file statistics"""
        stats = {
            'file_path': input_file,
            'file_size_mb': os.path.getsize(input_file) / (1024 * 1024),
            'is_compressed': str(input_file).endswith('.gz'),
            'md5_checksum': self._calculate_md5(input_file)
        }
        
        stats.update(self._analyze_fastq(input_file))
            
        return stats

    def _analyze_fastq(self, fastq_file):
        """
        Optimized FASTQ analysis with sampling.
        
        Analyzes up to MAX_READS_TO_ANALYZE (1M) reads for performance.
        Collects comprehensive quality and contamination statistics.
        
        Args:
            fastq_file: Path to FASTQ file (can be gzipped)
            
        Returns:
            Dict with statistics
        """
        from .output_formatter import get_formatter
        fmt = get_formatter()
        
        stats = {
            'total_sequences': 0,
            'total_bases': np.int64(0),
            'min_length': float('inf'),
            'max_length': 0,
            'mean_length': 0,
            'median_length': 0,
            'n50_length': 0,
            'gc_content': 0,
            'mean_quality': 0,
            'min_quality': float('inf'),
            'max_quality': 0,
            'q20_bases': 0,
            'q30_bases': 0,
            'quality_encoding': None,
            'adapter_content_percent': 0,
            'overrepresented_sequences': [],
            'top_5_sequences': [],
            'sampling_limited': False,
        }
        
        # Open file handle
        handle = gzip.open(fastq_file, 'rt') if str(fastq_file).endswith('.gz') else open(fastq_file, 'r')
        
        # Initialize collectors
        lengths = []
        gc_count = 0
        total_quality_sum = 0
        sequence_counts = Counter()
        adapter_hits = 0
        sampled_sequences = 0
        hit_sample_limit = False
        
        try:
            for record in SeqIO.parse(handle, 'fastq'):
                stats['total_sequences'] += 1
                
                # Performance optimization: Stop after 1M reads
                if stats['total_sequences'] > self.MAX_READS_TO_ANALYZE:
                    stats['total_sequences'] -= 1
                    hit_sample_limit = True
                    break
                
                # Sequence analysis
                seq_str = str(record.seq)
                seq_len = len(seq_str)
                lengths.append(seq_len)
                
                stats['min_length'] = min(stats['min_length'], seq_len)
                stats['max_length'] = max(stats['max_length'], seq_len)
                stats['total_bases'] += np.int64(seq_len)
                
                # GC content
                gc_count += seq_str.count('G') + seq_str.count('C')
                
                # Quality scores
                qual_scores = record.letter_annotations.get("phred_quality", [])
                if qual_scores:
                    total_quality_sum += sum(qual_scores)
                    stats['min_quality'] = min(stats['min_quality'], min(qual_scores))
                    stats['max_quality'] = max(stats['max_quality'], max(qual_scores))
                    stats['q20_bases'] += sum(1 for q in qual_scores if q >= 20)
                    stats['q30_bases'] += sum(1 for q in qual_scores if q >= 30)
                
                # Sample sequences for contamination analysis
                if sampled_sequences < self.sample_size:
                    sequence_counts[seq_str] += 1
                    
                    # Adapter detection
                    for adapter in self.common_adapters:
                        if adapter in seq_str:
                            adapter_hits += 1
                            break
                    
                    sampled_sequences += 1
        
        finally:
            handle.close()
        
        # Report sampling status
        if hit_sample_limit:
            stats['sampling_limited'] = True
            fmt.info(f"Performance optimization: Analyzed first {self.MAX_READS_TO_ANALYZE:,} reads (sampling mode)")
            fmt.info("Reported quality statistics describe the sampled reads, not the full file")
        
        # Calculate derived statistics
        if lengths:
            stats['mean_length'] = np.mean(lengths)
            stats['median_length'] = np.median(lengths)
            stats['n50_length'] = self._calculate_n50(lengths)
        
        if stats['total_bases'] > 0:
            stats['gc_content'] = (gc_count / stats['total_bases']) * 100
            stats['mean_quality'] = total_quality_sum / stats['total_bases']
        
        # Bio.SeqIO's "fastq" parser interprets qualities as Phred+33.
        stats['quality_encoding'] = "Phred+33 parser"
        
        # Adapter contamination percentage
        if sampled_sequences > 0:
            stats['adapter_content_percent'] = (adapter_hits / sampled_sequences) * 100
        
        # Overrepresented sequences
        for seq, count in sequence_counts.most_common(5):
            percentage = (count / sampled_sequences) * 100
            stats['top_5_sequences'].append((seq, count, percentage))
            if percentage >= self.overrep_threshold * 100:
                stats['overrepresented_sequences'].append((seq, count, percentage))
        
        return stats

       
    def _display_statistics(self, stats):
        """Display statistics using beautiful tables"""
        self._display_fastq_stats_table(stats)
    
    def _display_fastq_stats_table(self, stats):
        """Display FASTQ stats in beautiful table format"""
        from .output_formatter import TableFormatter, get_formatter
        fmt = get_formatter()
        
        # Build sections for table
        sections = [
            {
                'header': 'File Properties',
                'rows': {
                    'Size': f"{stats['file_size_mb']:.2f} MB ({'gzip compressed' if stats['is_compressed'] else 'uncompressed'})",
                    'MD5 Checksum': stats['md5_checksum'],
                    'Encoding': stats.get('quality_encoding', 'Unknown')
                }
            },
            {
                'header': 'Sequence Quality',
                'rows': {
                    'Reads analyzed': f"{stats['total_sequences']:,} reads",
                    'Total Bases': f"{stats['total_bases'] / 1_000_000:.1f} Mbp (mean: {stats['mean_length']:.0f} bp)",
                    'Length Range': f"{stats['min_length']}-{stats['max_length']} bp (median: {stats['median_length']:.0f} bp)",
                    'N50': f"{stats['n50_length']:.0f} bp",
                    'GC Content': f"{stats['gc_content']:.1f}%"
                }
            }
        ]
        
        # Display main table
        table = TableFormatter.format_table("\U0001f4ca SUMMARY STATISTICS", sections, width=73)
        fmt._print(f"\n{table}")
        
        # Display quality metrics with visual bars
        q20_pct = (stats['q20_bases'] / stats['total_bases'] * 100) if stats['total_bases'] > 0 else 0
        q30_pct = (stats['q30_bases'] / stats['total_bases'] * 100) if stats['total_bases'] > 0 else 0
        mean_q = stats['mean_quality']
        
        fmt._print("\n     \U0001f4c8 Quality Visualization:")
        fmt._print(f"        {'Mean quality score':.<30} Q{mean_q:.1f}")
        quality_metrics = {
            'Bases \u2265 Q20': (q20_pct, ""),
            'Bases \u2265 Q30': (q30_pct, "")
        }
        
        for metric, (pct, extra) in quality_metrics.items():
            bar = TableFormatter.format_visual_bar(pct, width=10, label=extra)
            fmt._print(f"        {metric:.<30} {bar}")
        
        # Contamination info
        fmt._print("\n     \U0001f52c Contamination Analysis:")
        fmt._print(f"        Adapter Content: {stats['adapter_content_percent']:.2f}% (sampled {min(self.sample_size, stats['total_sequences'])} reads)")
        
        if stats['overrepresented_sequences']:
            fmt._print(f"        \u26a0\ufe0f  Overrepresented Sequences: {len(stats['overrepresented_sequences'])} detected (\u2265{self.overrep_threshold:.1%})")
        else:
            fmt._print("        \u2713 No overrepresented sequences detected")
    
    def _get_validation_status(self, stats):
        """Return fatal structural errors and nonfatal quality warnings."""
        errors = []
        warnings = []

        if stats['total_sequences'] < self.min_sequences:
            errors.append(f"Insufficient sequences: {stats['total_sequences']} < {self.min_sequences} required")
            
        if stats['mean_quality'] < self.quality_threshold:
            warnings.append(
                f"Low quality score: Q{stats['mean_quality']:.1f} "
                f"< Q{self.quality_threshold} recommended"
            )
            
        q30_pct = (stats.get('q30_bases', 0) / stats.get('total_bases', 1) * 100)
        if q30_pct < self.q30_threshold:
            warnings.append(
                f"Low Q30 content: {q30_pct:.1f}% < {self.q30_threshold}% recommended"
            )
            
        if stats.get('adapter_content_percent', 0) > self.adapter_threshold:
            warnings.append(
                f"High adapter contamination: {stats['adapter_content_percent']:.2f}% "
                f"> {self.adapter_threshold}% threshold"
            )
            
        if stats.get('overrepresented_sequences'):
            max_overrep = max(pct for _, _, pct in stats['overrepresented_sequences'])
            warnings.append(
                f"Overrepresented sequences: {len(stats['overrepresented_sequences'])} "
                f"detected (max {max_overrep:.1f}%)"
            )
        
        if stats['total_sequences'] == 0:
            errors.append("No valid sequences found")

        if self.strict:
            errors.extend(warnings)
            warnings = []
        return errors, warnings

    @staticmethod
    def validate_pair(forward, reverse) -> tuple[bool, str]:
        """Fully verify paired FASTQ counts, ordering, and read identifiers."""
        def records(path):
            opener = gzip.open if str(path).endswith(".gz") else open
            with opener(path, "rt", encoding="utf-8") as handle:
                yield from SeqIO.parse(handle, "fastq")

        def normalized(identifier):
            return identifier.split()[0].removesuffix("/1").removesuffix("/2")

        for index, (left, right) in enumerate(
            zip_longest(records(forward), records(reverse)), 1
        ):
            if left is None or right is None:
                return False, f"Paired FASTQ files have different read counts at pair {index}"
            if normalized(left.id) != normalized(right.id):
                return False, (
                    f"Paired FASTQ identifiers differ at pair {index}: "
                    f"{left.id!r} != {right.id!r}"
                )
        return True, ""
    
    def _print_recommendations(self, stats, error_messages):
        """Print recommendations for quality improvement"""
        from .output_formatter import get_formatter
        fmt = get_formatter()
        fmt._print("\n     Quality preprocessing recommended:")
        fmt._print("\n     Option 1: fastp (recommended)")
        fmt._print("        fastp -i input.fq -o output.fq --qualified_quality_phred 20 \\")
        fmt._print("              --detect_adapter_for_pe --cut_right")
        fmt._print("\n     Option 2: Trimmomatic")
        fmt._print("        trimmomatic SE input.fq output.fq LEADING:20 TRAILING:20 \\")
        fmt._print("                    SLIDINGWINDOW:4:20 MINLEN:50")
    
    def _calculate_n50(self, lengths):
        """Calculate N50 metric"""
        if not lengths:
            return 0
        sorted_lengths = sorted(lengths, reverse=True)
        cumsum = np.cumsum(sorted_lengths)
        total = cumsum[-1]
        for i, cs in enumerate(cumsum):
            if cs >= total / 2:
                return sorted_lengths[i]
        return 0
    
    def _calculate_md5(self, filepath):
        """Calculate MD5 checksum"""
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def _check_file_exists(self, filepath):
        """Check file existence"""
        from .output_formatter import get_formatter
        if not Path(filepath).exists():
            get_formatter().error(f"File not found: {filepath}")
            return False
        return True
    
    def _validate_format(self, filepath, expected_type):
        """Validate every FASTQ record, not only the sampled QC prefix."""
        from .output_formatter import get_formatter
        try:
            opener = gzip.open if str(filepath).endswith('.gz') else open
            with opener(filepath, 'rt') as handle:
                records = sum(1 for _ in SeqIO.parse(handle, 'fastq'))
            if records == 0:
                raise StopIteration
            return True
            
        except StopIteration:
            get_formatter().error("File is empty")
            return False
        except Exception as e:
            get_formatter().error(f"Invalid FASTQ format: {e}")
            return False
