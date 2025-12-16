#!/usr/bin/env python3
"""
ATAC-seq preprocessing for Hi-Compass training data preparation.

This module handles:
1. Stratified subsampling from bulk ATAC-seq BAM files
2. Conversion to BigWig format with proper chromosome filtering
3. Hierarchical directory organization for training data
"""

import os
import subprocess
from pathlib import Path
from typing import List, Optional, Union
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ATACPreprocessor:
    """Preprocesses ATAC-seq data for Hi-Compass training."""
    
    def __init__(
        self,
        input_bam: Union[str, Path],
        cell_type: str,
        output_dir: Union[str, Path],
        chrom_sizes: Union[str, Path],
        genome: str = 'hg38',
        depth_mode: str = 'range',
        depth_list: Optional[List[Union[int, float]]] = None,
        min_depth: float = 2e5,
        max_depth: float = 2e7,
        step: float = 2e4,
        seed: int = 114514,
        keep_intermediate: bool = False
    ):
        """
        Initialize ATAC preprocessor.
        
        Args:
            input_bam: Path to input BAM/SAM file (bulk ATAC-seq)
            cell_type: Cell line name (e.g., 'GM12878', 'K562')
            output_dir: Output directory for processed files
            chrom_sizes: Chromosome sizes file (e.g., hg38.chrom.sizes)
            genome: Genome build ('hg38', 'mm10', etc.) for chromosome filtering
            depth_mode: 'range' or 'list' for depth specification
            depth_list: Custom depths (for 'list' mode)
            min_depth: Minimum depth (for 'range' mode, default: 2e5)
            max_depth: Maximum depth (for 'range' mode, default: 2e7)
            step: Step size (for 'range' mode, default: 2e4)
            seed: Random seed for reproducibility
            keep_intermediate: Keep intermediate files (BAM, bedGraph)
        """
        self.input_bam = Path(input_bam)
        self.cell_type = self._sanitize_cell_type(cell_type)
        self.output_dir = Path(output_dir)
        self.chrom_sizes = Path(chrom_sizes)
        self.genome = genome
        
        self.depth_mode = depth_mode
        self.depth_list = [int(d) for d in depth_list] if depth_list else None
        self.min_depth = int(min_depth)
        self.max_depth = int(max_depth)
        self.step = int(step)
        self.seed = seed
        self.keep_intermediate = keep_intermediate
        
        # Define valid chromosomes based on genome
        self._set_valid_chromosomes()
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._validate_inputs()
        self._total_reads = None
    
    def _set_valid_chromosomes(self):
        """Set valid chromosome names based on genome build."""
        if self.genome == 'hg38':
            # Human: chr1-22, chrX, chrY
            self.valid_chroms = set([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'])
        elif self.genome == 'mm10':
            # Mouse: chr1-19, chrX, chrY
            self.valid_chroms = set([f'chr{i}' for i in range(1, 20)] + ['chrX', 'chrY'])
        else:
            # For custom genomes, read from chrom_sizes file
            logger.warning(f"Custom genome '{self.genome}' - will use all chromosomes from chrom_sizes")
            self.valid_chroms = None  # Will be set from chrom_sizes file
    
    @staticmethod
    def _sanitize_cell_type(cell_type: str) -> str:
        """Sanitize cell type to valid filename characters."""
        import re
        return re.sub(r'[^a-zA-Z0-9_-]', '_', cell_type)
    
    def _validate_inputs(self):
        """Validate inputs and check required tools."""
        if not self.input_bam.exists():
            raise FileNotFoundError(f"Input BAM not found: {self.input_bam}")
        if not self.chrom_sizes.exists():
            raise FileNotFoundError(f"Chrom sizes not found: {self.chrom_sizes}")
        
        # Load valid chromosomes from chrom_sizes if not set
        if self.valid_chroms is None:
            self.valid_chroms = set()
            with open(self.chrom_sizes, 'r') as f:
                for line in f:
                    chrom = line.strip().split('\t')[0]
                    self.valid_chroms.add(chrom)
            logger.info(f"Loaded {len(self.valid_chroms)} chromosomes from {self.chrom_sizes.name}")
        
        # Validate depth configuration
        if self.depth_mode == 'list':
            if not self.depth_list:
                raise ValueError("depth_list required for mode 'list'")
            if any(d <= 0 for d in self.depth_list):
                raise ValueError("All depths must be positive")
        elif self.depth_mode == 'range':
            if self.min_depth >= self.max_depth:
                raise ValueError("min_depth must be < max_depth")
        else:
            raise ValueError(f"Invalid depth_mode: {self.depth_mode}")
        
        # Check required tools
        required_tools = ['samtools', 'bedtools', 'bedGraphToBigWig']
        missing_tools = []
        for tool in required_tools:
            if subprocess.run(['which', tool], capture_output=True).returncode != 0:
                missing_tools.append(tool)
        
        if missing_tools:
            raise EnvironmentError(
                f"Required tools not found: {', '.join(missing_tools)}\n"
                f"Please install: samtools, bedtools, bedGraphToBigWig (UCSC tools)"
            )
    
    def _get_total_reads(self) -> int:
        """Count total reads in BAM file."""
        if self._total_reads is None:
            logger.info("Counting reads in input BAM...")
            result = subprocess.run(
                ['samtools', 'view', '-c', str(self.input_bam)],
                capture_output=True, text=True, check=True
            )
            self._total_reads = int(result.stdout.strip())
            logger.info(f"Total reads in {self.input_bam.name}: {self._total_reads:,}")
        return self._total_reads
    
    @staticmethod
    def _format_depth(depth: int) -> str:
        """
        Format depth as scientific notation (e.g., 1e5, 2e6, 2.1e5).
        
        Examples:
            100000 -> '1e5'
            210000 -> '2.1e5'
            2000000 -> '2e6'
        """
        exp = 0
        value = float(depth)
        while value >= 10:
            value /= 10
            exp += 1
        
        if value == int(value):
            return f"{int(value)}e{exp}"
        return f"{value:.1f}e{exp}".rstrip('0').rstrip('.')
    
    def _generate_filename(self, depth: Union[int, str], extension: str) -> str:
        """
        Generate standardized filename: {CellType}~ATAC~{Depth}.{ext}
        
        Examples:
            GM12878, 210000, 'bw' -> 'GM12878~ATAC~2.1e5.bw'
            K562, 'bulk', 'bam' -> 'K562~ATAC~bulk.bam'
        """
        if isinstance(depth, int):
            depth_str = self._format_depth(depth)
        else:
            depth_str = depth
        return f"{self.cell_type}~ATAC~{depth_str}.{extension}"
    
    def _get_depth_dir(self, depth: Union[int, str]) -> Path:
        """
        Get subdirectory for a specific depth.
        
        Creates hierarchical structure:
        output_dir/
        └── CellType~ATAC~Depth/
            ├── CellType~ATAC~Depth.bam
            ├── CellType~ATAC~Depth.sorted.bam
            ├── CellType~ATAC~Depth.bedgraph
            ├── CellType~ATAC~Depth.sorted.bedgraph
            └── CellType~ATAC~Depth.bw
        """
        if isinstance(depth, int):
            depth_str = self._format_depth(depth)
        else:
            depth_str = depth
        
        depth_dir = self.output_dir / f"{self.cell_type}~ATAC~{depth_str}"
        depth_dir.mkdir(parents=True, exist_ok=True)
        return depth_dir
    
    def get_target_depths(self) -> List[int]:
        """Get list of target depths based on configuration."""
        if self.depth_mode == 'list':
            depths = sorted(self.depth_list)
        else:
            depths = list(range(self.min_depth, self.max_depth + 1, self.step))
        
        logger.info(f"Target depths: {len(depths)} levels")
        return depths
    
    def subsample_bam(self, target_depth: int) -> Path:
        """
        Subsample BAM to target depth using samtools.
        
        Args:
            target_depth: Target number of reads
            
        Returns:
            Path to subsampled BAM in depth-specific subdirectory
        """
        total_reads = self._get_total_reads()
        
        if target_depth >= total_reads:
            logger.warning(
                f"Target {target_depth:,} >= total {total_reads:,}, using full BAM"
            )
            return self.input_bam
        
        fraction = target_depth / total_reads
        
        # Get depth-specific directory
        depth_dir = self._get_depth_dir(target_depth)
        output_bam = depth_dir / self._generate_filename(target_depth, 'bam')
        
        logger.info(f"Subsampling to {target_depth:,} reads (fraction={fraction:.6f})")
        
        # Format for samtools -s: seed.fraction
        fraction_str = f"{fraction:.10f}"
        if fraction_str.startswith('0.'):
            fraction_decimal = fraction_str[2:]  # Remove "0."
        else:
            fraction_decimal = fraction_str
        
        fraction_decimal = fraction_decimal.rstrip('0') or '0'
        samtools_seed = f"{self.seed}.{fraction_decimal}"
        
        try:
            subprocess.run([
                'samtools', 'view',
                '-s', samtools_seed,
                '-b', '-o', str(output_bam),
                str(self.input_bam)
            ], check=True, capture_output=True, text=True)
            
            logger.info(f"✓ Created: {output_bam.relative_to(self.output_dir)}")
            return output_bam
            
        except subprocess.CalledProcessError as e:
            logger.error(
                f"samtools failed: -s {samtools_seed}\n"
                f"stderr: {e.stderr}"
            )
            raise
    
    def _filter_bedgraph(self, input_bedgraph: Path, output_bedgraph: Path):
        """
        Filter bedGraph to keep only valid chromosomes.
        
        This is critical to prevent bedGraphToBigWig errors from non-standard
        chromosomes (e.g., chrM, random contigs, alt haplotypes).
        
        Args:
            input_bedgraph: Input bedGraph file (sorted)
            output_bedgraph: Output filtered bedGraph file
        """
        logger.info(f"Filtering chromosomes (keeping {len(self.valid_chroms)} standard chroms)...")
        
        filtered_lines = 0
        total_lines = 0
        
        with open(input_bedgraph, 'r') as infile, open(output_bedgraph, 'w') as outfile:
            for line in infile:
                total_lines += 1
                chrom = line.split('\t')[0]
                if chrom in self.valid_chroms:
                    outfile.write(line)
                    filtered_lines += 1
        
        kept_pct = (filtered_lines / total_lines * 100) if total_lines > 0 else 0
        logger.info(f"  Kept {filtered_lines:,}/{total_lines:,} lines ({kept_pct:.1f}%)")
    
    def bam_to_bigwig(self, bam_file: Path, depth: Union[int, str]) -> Path:
        """
        Convert BAM to BigWig with chromosome filtering.
        
        Pipeline:
        1. BAM -> sorted BAM (coordinate sorted)
        2. sorted BAM -> bedGraph (coverage)
        3. bedGraph -> sorted bedGraph (lexicographic sort)
        4. sorted bedGraph -> filtered bedGraph (standard chroms only)
        5. filtered bedGraph -> BigWig
        
        Args:
            bam_file: Input BAM file
            depth: Depth label (for directory structure)
            
        Returns:
            Path to output BigWig
        """
        depth_dir = self._get_depth_dir(depth)
        
        base = bam_file.stem
        sorted_bam = depth_dir / f"{base}.sorted.bam"
        bedgraph = depth_dir / f"{base}.bedgraph"
        sorted_bedgraph_temp = depth_dir / f"{base}.sorted.temp.bedgraph"
        sorted_bedgraph = depth_dir / f"{base}.sorted.bedgraph"
        bigwig = depth_dir / f"{base}.bw"
        
        try:
            # Step 1: Sort BAM by coordinates
            logger.info(f"[1/5] Sorting BAM...")
            subprocess.run([
                'samtools', 'sort', '-o', str(sorted_bam), str(bam_file)
            ], check=True, capture_output=True)
            
            # Step 2: BAM to bedGraph (coverage)
            logger.info("[2/5] Converting to bedGraph...")
            with open(bedgraph, 'w') as f:
                subprocess.run([
                    'bedtools', 'genomecov', '-ibam', str(sorted_bam), '-bg'
                ], stdout=f, check=True)
            
            # Step 3: Sort bedGraph
            logger.info("[3/5] Sorting bedGraph...")
            with open(sorted_bedgraph_temp, 'w') as f:
                subprocess.run([
                    'sort', '-k1,1', '-k2,2n', str(bedgraph)
                ], stdout=f, check=True)
            
            # Step 4: Filter to standard chromosomes
            logger.info("[4/5] Filtering standard chromosomes...")
            self._filter_bedgraph(sorted_bedgraph_temp, sorted_bedgraph)
            
            # Step 5: bedGraph to BigWig
            logger.info("[5/5] Converting to BigWig...")
            subprocess.run([
                'bedGraphToBigWig',
                str(sorted_bedgraph),
                str(self.chrom_sizes),
                str(bigwig)
            ], check=True, capture_output=True)
            
            logger.info(f"✓ Created: {bigwig.relative_to(self.output_dir)}")
            
            # Cleanup intermediate files
            if not self.keep_intermediate:
                for f in [sorted_bam, bedgraph, sorted_bedgraph_temp, sorted_bedgraph]:
                    f.unlink(missing_ok=True)
                if bam_file != self.input_bam and bam_file.parent == depth_dir:
                    bam_file.unlink(missing_ok=True)
            
            return bigwig
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Conversion failed: {e}")
            if hasattr(e, 'stderr') and e.stderr:
                logger.error(f"stderr: {e.stderr}")
            raise
    
    def process(self, include_bulk: bool = True) -> List[Path]:
        """
        Process all target depths and optionally generate bulk BigWig.
        
        Directory structure created:
        output_dir/
        ├── CellType~ATAC~2.1e5/
        │   └── CellType~ATAC~2.1e5.bw
        ├── CellType~ATAC~4.2e5/
        │   └── CellType~ATAC~4.2e5.bw
        └── CellType~ATAC~bulk/
            └── CellType~ATAC~1.7e8.bw  (actual depth in filename)
        
        Args:
            include_bulk: Also generate bulk BigWig from full BAM
            
        Returns:
            List of generated BigWig files
        """
        depths = self.get_target_depths()
        logger.info(f"\n{'='*70}")
        logger.info(f"Starting ATAC-seq preprocessing")
        logger.info(f"{'='*70}")
        logger.info(f"Cell type: {self.cell_type}")
        logger.info(f"Genome: {self.genome}")
        logger.info(f"Target depths: {len(depths)}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Structure: {self.cell_type}~ATAC~<depth>/")
        logger.info(f"{'='*70}\n")
        
        bigwig_files = []
        
        # Process subsampled depths
        for i, depth in enumerate(depths, 1):
            try:
                logger.info(f"\n[{i}/{len(depths)}] Processing depth: {depth:,}")
                logger.info(f"{'='*50}")
                bam = self.subsample_bam(depth)
                bw = self.bam_to_bigwig(bam, depth)
                bigwig_files.append(bw)
            except Exception as e:
                logger.error(f"Failed at depth {depth}: {e}")
                continue
        
        # Generate bulk BigWig
        if include_bulk:
            try:
                logger.info(f"\n[Bulk] Generating bulk BigWig...")
                logger.info(f"{'='*50}")
                
                # Get actual total reads
                bulk_depth = self._get_total_reads()
                logger.info(f"Bulk depth: {bulk_depth:,} reads")
                
                # Create bulk directory with 'bulk' label
                bulk_dir = self._get_depth_dir('bulk')
                
                # Generate filename with actual depth
                bulk_filename = self._generate_filename(bulk_depth, 'bam')
                bulk_bam = bulk_dir / bulk_filename
                
                # Copy original BAM to bulk directory
                logger.info("Copying input BAM to bulk directory...")
                subprocess.run(['cp', str(self.input_bam), str(bulk_bam)], check=True)
                logger.info(f"✓ Copied to: {bulk_bam.relative_to(self.output_dir)}")
                
                # Convert to BigWig (use actual depth for final filename)
                bulk_bw = self.bam_to_bigwig(bulk_bam, 'bulk')
                
                # Rename BigWig to include actual depth in filename
                # The bam_to_bigwig uses the BAM filename stem, so we need to rename
                final_bw_name = self._generate_filename(bulk_depth, 'bw')
                final_bw_path = bulk_dir / final_bw_name
                
                if bulk_bw != final_bw_path:
                    bulk_bw.rename(final_bw_path)
                    logger.info(f"✓ Renamed to: {final_bw_path.name}")
                    bulk_bw = final_bw_path
                
                bigwig_files.append(bulk_bw)
                
            except Exception as e:
                logger.error(f"Failed to generate bulk: {e}")
        
        # Summary
        logger.info(f"\n{'='*70}")
        logger.info(f"✓ Processing Complete")
        logger.info(f"{'='*70}")
        logger.info(f"Generated {len(bigwig_files)} BigWig files:")
        for bw in bigwig_files:
            logger.info(f"  • {bw.relative_to(self.output_dir)}")
        logger.info(f"\nOutput directory: {self.output_dir}")
        logger.info(f"{'='*70}\n")
        
        return bigwig_files


def preprocess_atac(
    input_bam: str,
    cell_type: str,
    output_dir: str,
    chrom_sizes: str,
    genome: str = 'hg38',
    depths: Optional[List[Union[int, float]]] = None,
    min_depth: float = 2e5,
    max_depth: float = 2e7,
    step: float = 2e4,
    include_bulk: bool = True,
    **kwargs
) -> List[Path]:
    """
    Convenience function for ATAC-seq preprocessing.
    
    Example usage:
        # Using depth range
        bigwigs = preprocess_atac(
            input_bam='GM12878.bam',
            cell_type='GM12878',
            output_dir='data/ATAC/hg38',
            chrom_sizes='hg38.chrom.sizes',
            genome='hg38',
            min_depth=2e5,
            max_depth=2e7,
            step=2e4
        )
        
        # Using custom depth list
        bigwigs = preprocess_atac(
            input_bam='K562.bam',
            cell_type='K562',
            output_dir='data/ATAC/hg38',
            chrom_sizes='hg38.chrom.sizes',
            genome='hg38',
            depths=[1e5, 2.1e5, 5e5, 1e6]
        )
    
    Args:
        input_bam: Input BAM file path
        cell_type: Cell line name (e.g., 'GM12878', 'K562')
        output_dir: Output directory
        chrom_sizes: Chromosome sizes file
        genome: Genome build ('hg38', 'mm10', etc.)
        depths: Custom depth list (if provided, overrides range parameters)
        min_depth: Minimum depth for range mode (default: 2e5)
        max_depth: Maximum depth for range mode (default: 2e7)
        step: Step size for range mode (default: 2e4)
        include_bulk: Generate bulk BigWig (default: True)
        **kwargs: Additional arguments for ATACPreprocessor
        
    Returns:
        List of generated BigWig file paths
    """
    preprocessor = ATACPreprocessor(
        input_bam=input_bam,
        cell_type=cell_type,
        output_dir=output_dir,
        chrom_sizes=chrom_sizes,
        genome=genome,
        depth_mode='list' if depths else 'range',
        depth_list=depths,
        min_depth=min_depth,
        max_depth=max_depth,
        step=step,
        **kwargs
    )
    
    return preprocessor.process(include_bulk=include_bulk)


