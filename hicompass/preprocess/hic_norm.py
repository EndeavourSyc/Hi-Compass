#!/usr/bin/env python3
"""
Hi-C normalization preprocessing for Hi-Compass training data preparation.

Step 1: Contrast stretching normalization of cool files
This step prepares Hi-C contact maps for model training by applying
contrast stretching to standardize signal intensity.
"""

import cooler
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Union
from functools import reduce
from skimage import exposure
import logging

logger = logging.getLogger(__name__)


# Default chromosome sizes (hg38)
HG38_CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636,
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345,
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

# Default chromosome sizes (mm10)
MM10_CHROM_SIZES = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
    'chrY': 91744698
}


def load_chrom_sizes(chrom_sizes_file: Union[str, Path]) -> Dict[str, int]:
    """
    Load chromosome sizes from file.
    
    Supports multiple formats:
    1. Two-column TSV (chrom\tsize)
    2. UCSC chrom.sizes format
    3. Cool/mcool file (extracts from bins)
    
    Args:
        chrom_sizes_file: Path to chromosome sizes file or cool file
        
    Returns:
        Dictionary mapping chromosome names to sizes
        
    Examples:
        Format 1 (TSV):
            chr1    248956422
            chr2    242193529
            
        Format 2 (cool file):
            Automatically extracts from cool bins table
    """
    path = Path(chrom_sizes_file)
    
    if not path.exists():
        raise FileNotFoundError(f"Chrom sizes file not found: {path}")
    
    # Try loading as cool file first
    if path.suffix in ['.cool', '.mcool']:
        try:
            if path.suffix == '.mcool':
                # For mcool, need to specify resolution
                raise ValueError(
                    "For mcool files, please extract chromosome sizes manually "
                    "or specify resolution in the path (e.g., file.mcool::resolutions/10000)"
                )
            
            clr = cooler.Cooler(str(path))
            bins_df = clr.bins()[:]
            
            # Group by chromosome and get max end position
            chrom_sizes = {}
            for chrom in bins_df['chrom'].unique():
                chrom_bins = bins_df[bins_df['chrom'] == chrom]
                chrom_sizes[chrom] = int(chrom_bins['end'].max())
            
            logger.info(f"Loaded chromosome sizes from cool file: {len(chrom_sizes)} chromosomes")
            return chrom_sizes
            
        except Exception as e:
            logger.warning(f"Failed to load as cool file: {e}")
    
    # Load as TSV file
    try:
        df = pd.read_csv(
            path,
            sep='\t',
            header=None,
            names=['chrom', 'size'],
            comment='#'
        )
        
        chrom_sizes = dict(zip(df['chrom'], df['size'].astype(int)))
        logger.info(f"Loaded chromosome sizes from file: {len(chrom_sizes)} chromosomes")
        return chrom_sizes
        
    except Exception as e:
        raise ValueError(f"Failed to parse chromosome sizes file: {e}")


def extract_chrom_sizes_from_cool(cool_file: Union[str, Path]) -> Dict[str, int]:
    """
    Extract chromosome sizes directly from a cool file.
    
    Args:
        cool_file: Path to cool file
        
    Returns:
        Dictionary mapping chromosome names to sizes
    """
    try:
        clr = cooler.Cooler(str(cool_file))
        bins_df = clr.bins()[:]
        
        chrom_sizes = {}
        for chrom in bins_df['chrom'].unique():
            chrom_bins = bins_df[bins_df['chrom'] == chrom]
            chrom_sizes[chrom] = int(chrom_bins['end'].max())
        
        return chrom_sizes
        
    except Exception as e:
        raise ValueError(f"Failed to extract chromosome sizes from cool file: {e}")


class HiCNormalizer:
    """
    Normalizes Hi-C cool files using contrast stretching.
    
    This preprocessing step standardizes Hi-C contact map intensities
    to prepare them for Hi-Compass model training.
    """
    
    def __init__(
        self,
        input_cool: Union[str, Path],
        output_cool: Union[str, Path],
        resolution: int = 10000,
        genome: str = 'hg38',
        chrom_sizes: Optional[Union[Dict[str, int], str, Path]] = None,
        percentile_min: float = 2.0,
        percentile_max: float = 98.0,
        sample_region_start: int = 4000,
        sample_region_size: int = 256,
        balance: bool = True,
        chr_list: Optional[List[str]] = None,
        target_scale: float = 100.0,
        force_resolution: bool = False
    ):
        """
        Initialize Hi-C normalizer.
        
        Args:
            input_cool: Input cool file path
            output_cool: Output cool file path (contrast-stretched)
            resolution: Resolution in bp (default: 10000)
            genome: Genome assembly ('hg38', 'mm10', or 'custom')
            chrom_sizes: Chromosome sizes. Can be:
                - Dict[str, int]: chromosome name -> size
                - str/Path: path to chrom.sizes file or cool file
                - None: uses genome default or extracts from input_cool
            percentile_min: Lower percentile for contrast stretching (default: 2.0)
            percentile_max: Upper percentile for contrast stretching (default: 98.0)
            sample_region_start: Start bin for sampling percentiles (default: 4000)
            sample_region_size: Size of sample region in bins (default: 256)
            balance: Use balanced matrix (default: True)
            chr_list: List of chromosomes to process. If None, uses all chromosomes
            target_scale: Target scale after normalization (default: 100.0)
            force_resolution: Skip resolution check warning (default: False)
            
        Note:
            The sample_region (start=4000, size=256) was empirically determined
            to work well for 98% of Hi-C datasets at 10kb resolution for both
            human and mouse genomes. These parameters are exposed for flexibility
            but the defaults are strongly recommended.
        """
        self.input_cool = Path(input_cool)
        self.output_cool = Path(output_cool)
        self.resolution = resolution
        self.balance = balance
        self.percentile_min = percentile_min
        self.percentile_max = percentile_max
        self.sample_region_start = sample_region_start
        self.sample_region_size = sample_region_size
        self.target_scale = target_scale
        self.force_resolution = force_resolution
        
        # Validate resolution (warning only, not blocking)
        if resolution != 10000 and not force_resolution:
            logger.warning(
                f"⚠️  Resolution {resolution} is not 10kb. "
                f"Hi-Compass is optimized for 10kb resolution only. "
                f"Results may be suboptimal. Use --force-resolution to suppress this warning."
            )
        
        # Load chromosome sizes
        self.chrom_sizes = self._load_chrom_sizes(genome, chrom_sizes)
        
        # Set chromosome list
        if chr_list is not None:
            self.chr_list = chr_list
        else:
            # Use all chromosomes from chrom_sizes
            self.chr_list = list(self.chrom_sizes.keys())
            logger.info(f"Processing all chromosomes: {', '.join(self.chr_list[:5])}...")
        
        # Filter chrom_sizes to only include chr_list
        self.chrom_sizes = {
            k: v for k, v in self.chrom_sizes.items() 
            if k in self.chr_list
        }
        
        # Validate input
        self._validate_inputs()
    
    def _load_chrom_sizes(
        self,
        genome: str,
        chrom_sizes: Optional[Union[Dict[str, int], str, Path]]
    ) -> Dict[str, int]:
        """Load chromosome sizes from various sources."""
        
        # Priority 1: User-provided chrom_sizes
        if chrom_sizes is not None:
            # If it's already a dict
            if isinstance(chrom_sizes, dict):
                logger.info(f"Using provided chromosome sizes dict: {len(chrom_sizes)} chromosomes")
                return chrom_sizes
            
            # If it's a file path
            else:
                logger.info(f"Loading chromosome sizes from: {chrom_sizes}")
                return load_chrom_sizes(chrom_sizes)
        
        # Priority 2: Built-in genome assemblies
        if genome == 'hg38':
            logger.info("Using built-in hg38 chromosome sizes")
            return HG38_CHROM_SIZES.copy()
        elif genome == 'mm10':
            logger.info("Using built-in mm10 chromosome sizes")
            return MM10_CHROM_SIZES.copy()
        
        # Priority 3: Extract from input cool file
        elif genome == 'custom':
            logger.info("Extracting chromosome sizes from input cool file...")
            try:
                chrom_sizes = extract_chrom_sizes_from_cool(self.input_cool)
                logger.info(f"Extracted {len(chrom_sizes)} chromosomes from cool file")
                return chrom_sizes
            except Exception as e:
                raise ValueError(
                    f"Failed to extract chromosome sizes from cool file. "
                    f"Please provide chrom_sizes explicitly. Error: {e}"
                )
        
        else:
            raise ValueError(
                f"Unknown genome: {genome}. "
                f"Use 'hg38', 'mm10', or 'custom' (with --chrom-sizes provided)"
            )
    
    def _validate_inputs(self):
        """Validate inputs."""
        if not self.input_cool.exists():
            raise FileNotFoundError(f"Input cool file not found: {self.input_cool}")
        
        # Check cooler library
        try:
            import cooler
            from skimage import exposure
        except ImportError as e:
            raise ImportError(f"Required library not found: {e}")
        
        # Validate cool file and resolution
        try:
            clr = cooler.Cooler(str(self.input_cool))
            actual_res = clr.binsize
            
            if actual_res != self.resolution:
                raise ValueError(
                    f"Cool file resolution ({actual_res}) does not match "
                    f"specified resolution ({self.resolution})"
                )
            
            logger.info(f"✓ Cool file validated: resolution={actual_res}")
            
            # Validate chromosomes exist in cool file
            cool_chroms = set(clr.chromnames)
            missing_chroms = set(self.chr_list) - cool_chroms
            if missing_chroms:
                logger.warning(
                    f"⚠️  Chromosomes not found in cool file: {missing_chroms}. "
                    f"These will be skipped."
                )
                self.chr_list = [c for c in self.chr_list if c in cool_chroms]
            
        except Exception as e:
            raise ValueError(f"Invalid cool file: {e}")
    
    @staticmethod
    def _get_chr_stack(
        chr_list: List[str],
        chr_name: str,
        chrom_sizes: Dict[str, int],
        resolution: int
    ) -> int:
        """Get cumulative bin offset for a chromosome."""
        chr_before_list = chr_list[:chr_list.index(chr_name)]
        if len(chr_before_list) == 0:
            return 0
        
        stack = 0
        for before_chr in chr_before_list:
            stack += int(chrom_sizes[before_chr] / resolution) + 1
        return stack
    
    @staticmethod
    def _make_cooler_bins_chr(
        chr_name: str,
        length: int,
        resolution: int
    ) -> pd.DataFrame:
        """Create bins DataFrame for a chromosome."""
        total_bin = []
        for i in range(0, length, resolution):
            total_bin.append([chr_name, i, min(length, i + resolution)])
        
        return pd.DataFrame(total_bin, columns=['chrom', 'start', 'end'])
    
    @staticmethod
    def _pseudo_weight(cool_path: Path, weight: float = 1.0):
        """Add pseudo weight column to cool file."""
        cooler_obj = cooler.Cooler(str(cool_path))
        
        stats = {
            "min_nnz": 0,
            "min_count": 0,
            "mad_max": 0,
            "cis_only": False,
            "ignore_diags": 2,
            "converged": False,
            "divisive_weights": False,
        }
        
        with cooler_obj.open("r+") as grp:
            if 'weight' in grp['bins']:
                grp['bins']['weight'].attrs.update(stats)
    
    def process_chromosome(
        self,
        cooler_obj: cooler.Cooler,
        chr_name: str,
        sample_percentiles: Optional[tuple] = None
    ) -> tuple:
        """Process a single chromosome with contrast stretching."""
        # Create bins
        bins_df = self._make_cooler_bins_chr(
            chr_name=chr_name,
            length=self.chrom_sizes[chr_name],
            resolution=self.resolution
        )
        
        # Fetch matrix
        matrix = cooler_obj.matrix(balance=self.balance, sparse=False).fetch(chr_name)
        
        # Clean diagonal and adjacent diagonals
        np.fill_diagonal(matrix, 0)
        np.fill_diagonal(matrix[1:, :], 0)
        np.fill_diagonal(matrix[:, 1:], 0)
        
        # Clip and clean
        matrix = np.clip(matrix, a_min=0, a_max=None)
        matrix = np.nan_to_num(matrix)
        
        # Compute percentiles if not provided
        if sample_percentiles is None:
            sample_end = self.sample_region_start + self.sample_region_size
            
            if sample_end <= matrix.shape[0]:
                sample_region = matrix[
                    self.sample_region_start:sample_end,
                    self.sample_region_start:sample_end
                ]
                p_min, p_max = np.percentile(
                    sample_region,
                    (self.percentile_min, self.percentile_max)
                )
            else:
                logger.warning(
                    f"{chr_name}: Sample region exceeds matrix size, "
                    f"using full matrix for percentiles"
                )
                p_min, p_max = np.percentile(
                    matrix,
                    (self.percentile_min, self.percentile_max)
                )
        else:
            p_min, p_max = sample_percentiles
        
        logger.info(
            f"{chr_name}: Percentiles ({self.percentile_min}%, {self.percentile_max}%) = "
            f"({p_min:.2f}, {p_max:.2f})"
        )
        
        # Contrast stretching on upper triangle
        matrix = np.triu(matrix)
        matrix_rescaled = exposure.rescale_intensity(
            matrix,
            in_range=(p_min, p_max)
        ) * self.target_scale
        
        # Convert to sparse format with global bin IDs
        rows, cols = np.nonzero(matrix_rescaled)
        stack = self._get_chr_stack(
            chr_list=self.chr_list,
            chr_name=chr_name,
            chrom_sizes=self.chrom_sizes,
            resolution=self.resolution
        )
        
        rows += stack
        cols += stack
        counts = matrix_rescaled[np.nonzero(matrix_rescaled)]
        
        pixels_df = pd.DataFrame({
            'bin1_id': rows,
            'bin2_id': cols,
            'count': counts
        })
        
        return bins_df, pixels_df
    
    def process(self):
        """Process all chromosomes and create output cool file."""
        logger.info(f"Processing {self.input_cool}")
        logger.info(f"Chromosomes: {', '.join(self.chr_list)}")
        logger.info(f"Resolution: {self.resolution}")
        logger.info(
            f"Contrast stretching: {self.percentile_min}%-{self.percentile_max}% "
            f"→ [0, {self.target_scale}]"
        )
        
        cooler_obj = cooler.Cooler(str(self.input_cool))
        
        chr_data = {}
        for chr_name in self.chr_list:
            try:
                logger.info(f"Processing {chr_name}...")
                bins_df, pixels_df = self.process_chromosome(cooler_obj, chr_name)
                chr_data[chr_name] = (bins_df, pixels_df)
            except Exception as e:
                logger.error(f"Failed to process {chr_name}: {e}")
                raise
        
        logger.info("Combining chromosomes...")
        
        all_bins = [data[0] for data in chr_data.values()]
        all_pixels = [data[1] for data in chr_data.values()]
        
        combined_bins = pd.concat(all_bins, axis=0, ignore_index=True)
        combined_bins.columns = ['chrom', 'start', 'end']
        combined_bins['weight'] = 1.0
        
        combined_pixels = pd.concat(all_pixels, axis=0, ignore_index=True)
        combined_pixels.columns = ['bin1_id', 'bin2_id', 'count']
        
        invalid = combined_pixels['bin1_id'] > combined_pixels['bin2_id']
        if invalid.any():
            first_error = invalid.idxmax()
            raise ValueError(
                f"Invalid pixel at index {first_error}: "
                f"bin1_id > bin2_id (should be upper triangle)"
            )
        
        logger.info(f"Creating output cool file: {self.output_cool}")
        self.output_cool.parent.mkdir(parents=True, exist_ok=True)
        
        cooler.create_cooler(
            cool_uri=str(self.output_cool),
            bins=combined_bins,
            pixels=combined_pixels
        )
        
        self._pseudo_weight(self.output_cool)
        
        logger.info(f"✓ Successfully created: {self.output_cool}")
        return self.output_cool


def normalize_hic(
    input_cool: str,
    output_cool: str,
    resolution: int = 10000,
    genome: str = 'hg38',
    **kwargs
) -> Path:
    """Convenience function for Hi-C normalization preprocessing."""
    normalizer = HiCNormalizer(
        input_cool=input_cool,
        output_cool=output_cool,
        resolution=resolution,
        genome=genome,
        **kwargs
    )
    
    return normalizer.process()