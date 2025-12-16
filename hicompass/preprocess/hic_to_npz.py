#!/usr/bin/env python3
"""
Hi-C to NPZ conversion for Hi-Compass training data preparation.

Step 2: Convert cool files to NPZ format with diagonal compression
This step extracts diagonals from Hi-C contact matrices and saves them
in a compressed format suitable for efficient model training.

Note:
    This step should be run AFTER preprocess-hic-norm. The normalization
    step sets all weights to 1.0, so the --balance flag has no effect here.
"""

import numpy as np
import cooler
from pathlib import Path
from typing import Union, Optional, List
import logging

logger = logging.getLogger(__name__)


class HiCToNPZConverter:
    """
    Converts Hi-C cool files to NPZ format with diagonal compression.
    
    Extracts diagonals from contact matrices and saves them as NPZ files,
    one per chromosome. This format is optimized for Hi-Compass training.
    
    Parameters are fixed for Hi-Compass:
    - Resolution: 10kb (required)
    - Window size: 256 diagonals (minimum recommended)
    - Balance: True (but has no effect after norm step)
    - Dtype: float16 (for compression)
    """
    
    def __init__(
        self,
        input_cool: Union[str, Path],
        output_dir: Union[str, Path],
        resolution: int = 10000,
        window_size: int = 256,
        balance: bool = True,
        chr_list: Optional[List[str]] = None
    ):
        """
        Initialize Hi-C to NPZ converter.
        
        Args:
            input_cool: Input cool/mcool file path (preferably after norm step)
            output_dir: Output directory for NPZ files
            resolution: Resolution in bp (default: 10000, strongly recommended)
            window_size: Number of diagonals to extract (default: 256, minimum)
            balance: Use balanced matrix (default: True, no effect after norm)
            chr_list: List of chromosomes to process. If None, processes all
            
        Note:
            If using a normalized cool file from preprocess-hic-norm, the balance
            parameter has no effect since all weights are set to 1.0.
        """
        self.input_cool = Path(input_cool)
        self.output_dir = Path(output_dir)
        self.resolution = resolution
        self.window_size = window_size
        self.balance = balance
        self.chr_list = chr_list
        
        # Fixed parameters optimized for Hi-Compass
        self.dtype = np.float16
        self.add_chr_prefix = True
        
        # Validate and warn
        self._check_parameters()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Validate inputs
        self._validate_inputs()
    
    def _check_parameters(self):
        """Check parameters and issue warnings if needed."""
        
        # Check resolution
        if self.resolution != 10000:
            logger.warning(
                f"\n{'='*70}\n"
                f"   WARNING: Resolution is {self.resolution}bp, not 10kb!\n"
                f"   Hi-Compass is trained and optimized for 10kb resolution.\n"
                f"   Using other resolutions may lead to poor performance.\n"
                f"{'='*70}\n"
            )
        
        # Check window size
        if self.window_size < 256:
            logger.warning(
                f"\n{'='*70}\n"
                f"   WARNING: Window size is {self.window_size}, less than 256!\n"
                f"   Hi-Compass training requires at least 256 diagonals.\n"
                f"   Smaller windows may result in insufficient information\n"
                f"   for the model to learn meaningful patterns.\n"
                f"{'='*70}\n"
            )
        
        # Info about balance parameter
        if self.balance:
            logger.info(
                "\n Note: If input is from preprocess-hic-norm, the --balance flag\n"
                "   has no effect since all weights are already set to 1.0."
            )
    
    def _validate_inputs(self):
        """Validate inputs."""
        if not self.input_cool.exists():
            raise FileNotFoundError(f"Input cool file not found: {self.input_cool}")
        
        # Check cooler library
        try:
            import cooler
        except ImportError:
            raise ImportError("cooler library not found. Install: pip install cooler")
        
        # Determine if this is a regular cool or mcool
        try:
            # Try loading as regular cool first
            clr = cooler.Cooler(str(self.input_cool))
            self.is_mcool = False
            self.cooler_path = str(self.input_cool)
            
        except (ValueError, KeyError):
            # Might be an mcool, need to specify resolution
            self.is_mcool = True
            self.cooler_path = f'{self.input_cool}::resolutions/{self.resolution}'
            
            try:
                clr = cooler.Cooler(self.cooler_path)
            except Exception as e:
                raise ValueError(
                    f"Failed to load cool file. If this is an mcool file, "
                    f"resolution {self.resolution} must exist in the file. Error: {e}"
                )
        
        # Validate resolution
        actual_res = clr.binsize
        if actual_res != self.resolution:
            raise ValueError(
                f"Cool file resolution ({actual_res}) does not match "
                f"specified resolution ({self.resolution})"
            )
        
        logger.info(
            f"✓ Cool file validated: {'mcool' if self.is_mcool else 'cool'}, "
            f"resolution={actual_res}"
        )
        
        # Get chromosome list if not provided
        if self.chr_list is None:
            self.chr_list = list(clr.chromnames)
            logger.info(f"Processing all chromosomes: {len(self.chr_list)} total")
        
        # Validate chromosomes
        available_chroms = set(clr.chromnames)
        missing = set(self.chr_list) - available_chroms
        if missing:
            logger.warning(f"Chromosomes not in cool file: {missing}")
            self.chr_list = [c for c in self.chr_list if c in available_chroms]
    
    @staticmethod
    def _compress_diagonals(matrix, window_size: int, dtype) -> dict:
        """
        Extract diagonals from matrix and compress.
        
        Args:
            matrix: Sparse or dense contact matrix
            window_size: Number of diagonals to extract
            dtype: NumPy dtype for compression
            
        Returns:
            Dictionary mapping diagonal index (str) to diagonal array
        """
        diag_dict = {}
        
        for d in range(window_size):
            # Upper diagonal
            diag_dict[str(d)] = np.nan_to_num(
                matrix.diagonal(d).astype(dtype)
            )
            
            # Lower diagonal (mirror)
            diag_dict[str(-d)] = np.nan_to_num(
                matrix.diagonal(-d).astype(dtype)
            )
        
        return diag_dict
    
    def _format_chrom_name(self, chrom: str) -> str:
        """Format chromosome name for output file."""
        # Always add 'chr' prefix if not present for UCSC compatibility
        formatted = chrom if chrom.startswith('chr') else f'chr{chrom}'
        return f"{formatted}.npz"
    
    def process_chromosome(self, cooler_obj: cooler.Cooler, chrom: str) -> Path:
        """
        Process a single chromosome and save to NPZ.
        
        Args:
            cooler_obj: Cooler object
            chrom: Chromosome name
            
        Returns:
            Path to output NPZ file
        """
        logger.info(f"Processing {chrom}...")
        
        # Fetch matrix
        matrix = cooler_obj.matrix(balance=self.balance, sparse=True).fetch(chrom)/10
        
        # Extract diagonals
        diag_dict = self._compress_diagonals(matrix, self.window_size, self.dtype)
        
        # Generate output filename
        output_filename = self._format_chrom_name(chrom)
        output_path = self.output_dir / output_filename
        
        # Save to NPZ (compressed)
        np.savez_compressed(output_path, **diag_dict)
        
        # Log file size
        file_size_mb = output_path.stat().st_size / (1024 * 1024)
        logger.info(f"✓ Saved {output_filename} ({file_size_mb:.2f} MB)")
        
        return output_path
    
    def process(self) -> List[Path]:
        """
        Process all chromosomes and save to NPZ files.
        
        Returns:
            List of paths to generated NPZ files
        """
        logger.info(f"\n{'='*70}")
        logger.info(f"Converting Hi-C to NPZ format")
        logger.info(f"{'='*70}")
        logger.info(f"Input:  {self.input_cool}")
        logger.info(f"Output: {self.output_dir}")
        logger.info(f"Resolution: {self.resolution} bp")
        logger.info(f"Window size: {self.window_size} diagonals")
        logger.info(f"Balance: {self.balance}")
        logger.info(f"Chromosomes: {len(self.chr_list)}")
        logger.info(f"Storage dtype: {self.dtype}")
        logger.info(f"{'='*70}\n")
        
        # Load cooler
        cooler_obj = cooler.Cooler(self.cooler_path)
        
        # Process each chromosome
        output_files = []
        for i, chrom in enumerate(self.chr_list, 1):
            try:
                logger.info(f"[{i}/{len(self.chr_list)}] {chrom}")
                output_path = self.process_chromosome(cooler_obj, chrom)
                output_files.append(output_path)
                
            except Exception as e:
                logger.error(f"Failed to process {chrom}: {e}")
                raise
        
        # Summary
        total_size = sum(f.stat().st_size for f in output_files)
        total_size_mb = total_size / (1024 * 1024)
        
        logger.info(f"\n{'='*70}")
        logger.info(f"✓ Conversion complete!")
        logger.info(f"  Files generated: {len(output_files)}")
        logger.info(f"  Total size: {total_size_mb:.2f} MB")
        logger.info(f"  Output directory: {self.output_dir}")
        logger.info(f"{'='*70}\n")
        
        return output_files


def hic_to_npz(
    input_cool: str,
    output_dir: str,
    **kwargs
) -> List[Path]:
    """
    Convenience function for Hi-C to NPZ conversion. An essential fucntion for further training.
    
    Args:
        input_cool: Input cool/mcool file path
        output_dir: Output directory for NPZ files
        **kwargs: Additional arguments for HiCToNPZConverter
        
    Returns:
        List of paths to generated NPZ files
    """
    converter = HiCToNPZConverter(
        input_cool=input_cool,
        output_dir=output_dir,
        **kwargs
    )
    
    return converter.process()