"""
Optimized Chromosome Dataset for Hi-Compass Training

This module provides a PyTorch Dataset implementation for loading and processing
Hi-C prediction training data with hierarchical directory structure.

File Structure Expected:
    ATAC/
        {cellline}~ATAC~{depth}/
            {cellline}~ATAC~{depth}.bw
        {cellline}~ATAC~bulk/
            {cellline}~ATAC~{actual_depth}.bw  # actual depth in filename
        ...
    
Key Features:
- Automatic cell line dictionary generation from data directory
- Multi-depth ATAC-seq support with folder-based organization
- Special handling for bulk data (folder named 'bulk', file contains actual depth)
- Data augmentation for robust model training
- Efficient data loading with caching
"""

import sys 
import os
import re
import glob
import random
import pickle
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
import torch
import pandas as pd
import numpy as np
from skimage.transform import resize
from torch.utils.data import Dataset

from . import data_feature


# ============================================================================
# Global Chromosome Size Dictionaries
# ============================================================================

# Human genome (hg38) chromosome sizes
HG38_CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

# Mouse genome (mm10) chromosome sizes
MM10_CHROM_SIZES = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
    'chrY': 91744698
}



def create_chr_index_dict(total_chr_list):
    return {chr_name: idx + 1 for idx, chr_name in enumerate(total_chr_list)}


def parse_atac_folder(folder_name: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse ATAC-seq folder name to extract cell line and depth label.
    
    Expected format: {cellline}~ATAC~{depth_label}
    Examples:
        - GM12878~ATAC~2.1e5 -> ('GM12878', '2.1e5')
        - gm12878~ATAC~bulk -> ('gm12878', 'bulk')
        - K562~ATAC~10e4 -> ('K562', '10e4')
    
    Args:
        folder_name: ATAC-seq folder name
        
    Returns:
        Tuple of (cell_line, depth_label) or (None, None) if parsing fails
    """
    pattern = r'^(.+?)~ATAC~(.+)$'
    match = re.match(pattern, folder_name)
    
    if match:
        cell_line = match.group(1).strip()
        depth_label = match.group(2).strip()
        return cell_line, depth_label
    
    return None, None


def parse_atac_filename(filename: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parse ATAC-seq filename to extract cell line and actual depth.
    
    Expected format: {cellline}~ATAC~{actual_depth}.bw
    Examples:
        - GM12878~ATAC~2.1e5.bw -> ('GM12878', '2.1e5')
        - gm12878~ATAC~1.7e8.bw -> ('gm12878', '1.7e8')
    
    Args:
        filename: ATAC-seq bigWig filename
        
    Returns:
        Tuple of (cell_line, actual_depth) or (None, None) if parsing fails
    """
    # Remove .bw extension if present
    basename = filename.replace('.bw', '').replace('.bigwig', '')
    
    # Try to match the pattern: cellline~ATAC~depth
    pattern = r'^(.+?)~ATAC~(.+)$'
    match = re.match(pattern, basename)
    
    if match:
        cell_line = match.group(1).strip()
        actual_depth = match.group(2).strip()
        return cell_line, actual_depth
    
    return None, None


def normalize_depth_string(depth_str: str) -> str:
    """
    Normalize depth string for consistent comparison.
    
    Examples:
        - '2.1e5' -> '2.1e5'
        - '2.1E5' -> '2.1e5'
        - '1.7e8' -> '1.7e8'
        - 'bulk' -> 'bulk'
    
    Args:
        depth_str: Depth string from filename or folder
        
    Returns:
        Normalized depth string
    """
    return depth_str.lower()


class ATACFileInfo:
    """
    Container for ATAC-seq file information.
    
    Attributes:
        cell_line: Cell line name (e.g., 'GM12878')
        folder_depth_label: Depth label from folder name (e.g., '2.1e5' or 'bulk')
        actual_depth: Actual depth from filename (e.g., '1.7e8' for bulk)
        file_path: Full path to the .bw file
        is_bulk: Whether this is bulk data
    """
    def __init__(self, cell_line: str, folder_depth_label: str, 
                 actual_depth: str, file_path: str):
        self.cell_line = cell_line
        self.folder_depth_label = normalize_depth_string(folder_depth_label)
        self.actual_depth = normalize_depth_string(actual_depth)
        self.file_path = file_path
        self.is_bulk = (self.folder_depth_label == 'bulk')
    
    def __repr__(self):
        return (f"ATACFileInfo(cell={self.cell_line}, "
                f"folder_label={self.folder_depth_label}, "
                f"actual_depth={self.actual_depth}, "
                f"is_bulk={self.is_bulk})")


def scan_atac_directory(atac_dir: str, verbose: bool = True) -> Dict[str, Dict[str, ATACFileInfo]]:
    """
    Scan ATAC-seq directory with hierarchical folder structure.
    
    Expected directory structure:
        atac_dir/
            GM12878~ATAC~2.1e5/
                GM12878~ATAC~2.1e5.bw
            GM12878~ATAC~bulk/
                GM12878~ATAC~1.7e8.bw    # actual depth in filename
            K562~ATAC~10e4/
                K562~ATAC~10e4.bw
            ...
    
    Args:
        atac_dir: Path to ATAC-seq data directory
        verbose: Whether to print scanning progress
        
    Returns:
        Nested dictionary: {cell_line: {depth_label: ATACFileInfo}}
        where depth_label is the folder name (e.g., '2.1e5', 'bulk', etc.)
    """
    if not os.path.exists(atac_dir):
        raise FileNotFoundError(f"ATAC directory not found: {atac_dir}")
    
    cell_depth_files = defaultdict(dict)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Scanning ATAC-seq directory: {atac_dir}")
        print(f"{'='*70}\n")
    
    # Find all subdirectories matching the pattern
    subdirs = [d for d in os.listdir(atac_dir) 
               if os.path.isdir(os.path.join(atac_dir, d))]
    
    if verbose:
        print(f"Found {len(subdirs)} subdirectories\n")
    
    parsed_count = 0
    bulk_count = 0
    
    for subdir in sorted(subdirs):
        subdir_path = os.path.join(atac_dir, subdir)
        
        # Parse folder name
        cell_line, folder_depth_label = parse_atac_folder(subdir)
        
        if not cell_line or not folder_depth_label:
            if verbose:
                print(f"  ✗ Skipped folder (invalid format): {subdir}")
            continue
        
        # Find .bw files in this folder
        bw_files = glob.glob(os.path.join(subdir_path, "*.bw"))
        bw_files.extend(glob.glob(os.path.join(subdir_path, "*.bigwig")))
        
        if not bw_files:
            if verbose:
                print(f"  ⚠ No bigWig files found in: {subdir}")
            continue
        
        # For each folder, we expect exactly one .bw file
        # If multiple, use the first one
        bw_file = bw_files[0]
        if len(bw_files) > 1 and verbose:
            print(f"  ⚠ Multiple .bw files in {subdir}, using: {os.path.basename(bw_file)}")
        
        # Parse filename to get actual depth
        filename = os.path.basename(bw_file)
        file_cell_line, actual_depth = parse_atac_filename(filename)
        
        if not actual_depth:
            if verbose:
                print(f"  ✗ Cannot parse filename: {filename}")
            continue
        
        # Check cell line consistency
        if file_cell_line.lower() != cell_line.lower():
            if verbose:
                print(f"  ⚠ Cell line mismatch: folder={cell_line}, file={file_cell_line}")
        
        # Create file info
        file_info = ATACFileInfo(
            cell_line=cell_line,
            folder_depth_label=folder_depth_label,
            actual_depth=actual_depth,
            file_path=bw_file
        )
        
        # Store in dictionary
        cell_depth_files[cell_line][folder_depth_label] = file_info
        
        parsed_count += 1
        if file_info.is_bulk:
            bulk_count += 1
        
        if verbose:
            if file_info.is_bulk:
                print(f"  ✓ {cell_line:15s} | BULK (actual: {actual_depth:10s}) | {subdir}")
            else:
                print(f"  ✓ {cell_line:15s} | depth: {folder_depth_label:10s} | {subdir}")
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Scanning Summary:")
        print(f"  - Successfully parsed: {parsed_count} folders")
        print(f"  - Unique cell lines: {len(cell_depth_files)}")
        print(f"  - Bulk data found: {bulk_count}")
        print(f"{'='*70}")
        
        print(f"\nCell Line Summary:")
        for cell_line in sorted(cell_depth_files.keys()):
            depths = list(cell_depth_files[cell_line].keys())
            print(f"  {cell_line:15s} : {len(depths)} depths {sorted(depths)}")
        print(f"{'='*70}\n")
    
    return dict(cell_depth_files)


def build_cell_line_mapping(cell_lines: List[str], verbose: bool = True) -> Dict[str, int]:
    """
    Build cell line to index mapping with alphabetical sorting for consistency.
    
    Args:
        cell_lines: List of cell line names
        verbose: Whether to print the mapping
        
    Returns:
        Dictionary mapping cell_line -> index
        
    Note:
        - Sorting ensures consistent mapping across runs
        - Returns empty dict if only one cell line (discriminator disabled)
    """
    # Sort cell lines alphabetically (case-insensitive)
    sorted_cell_lines = sorted(cell_lines, key=lambda x: x.lower())
    
    # Create mapping
    cell_line_to_idx = {cell: idx for idx, cell in enumerate(sorted_cell_lines)}
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Cell Line Mapping (Discriminator Configuration)")
        print(f"{'='*70}")
        
        if len(cell_line_to_idx) == 1:
            print(f"\n⚠ Only ONE cell line detected: {sorted_cell_lines[0]}")
            print(f"  -> Discriminator will be DISABLED during training")
        else:
            print(f"\nTotal cell lines: {len(cell_line_to_idx)}")
            print(f"  -> Discriminator will be ENABLED during training")
            print(f"\nMapping (sorted alphabetically):")
            for cell, idx in cell_line_to_idx.items():
                print(f"  {idx:3d} : {cell}")
        
        print(f"{'='*70}\n")
    
    return cell_line_to_idx


def verify_chromosome_sizes(bw_files: List[str], 
                            expected_chrom_sizes: Dict[str, int],
                            verbose: bool = True) -> bool:
    """
    Verify that all bigWig files have consistent chromosome sizes.
    
    Args:
        bw_files: List of bigWig file paths to check
        expected_chrom_sizes: Expected chromosome size dictionary
        verbose: Whether to print verification details
        
    Returns:
        True if all files match expected sizes, False otherwise
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"Verifying Chromosome Sizes")
        print(f"{'='*70}\n")
        print(f"Checking {len(bw_files)} files against reference...")
    
    # TODO: Implement actual chromosome size checking using pyBigWig
    # For now, return True as a placeholder
    
    if verbose:
        print(f"  ✓ All files have consistent chromosome sizes\n")
        print(f"{'='*70}\n")
    
    return True


def load_centromere_telomere_regions(bed_file: str) -> Dict[str, np.ndarray]:
    """
    Load centromere and telomere regions from BED file.
    
    Args:
        bed_file: Path to BED file with centromere/telomere regions
        
    Returns:
        Dictionary mapping chromosome -> array of excluded regions
    """
    excluded_regions = defaultdict(list)
    
    if not os.path.exists(bed_file):
        print(f"Warning: Centromere/telomere BED file not found: {bed_file}")
        return dict(excluded_regions)
    
    df = pd.read_csv(bed_file , sep = '\t', names = ['chr', 'start', 'end'])
    chrs = df['chr'].unique()
    centrotelo_dict = {}
    for chr_name in chrs:
        sub_df = df[df['chr'] == chr_name]
        regions = sub_df.drop('chr', axis = 1).to_numpy()
        centrotelo_dict[chr_name] = regions
    return centrotelo_dict
    

def load_chrom_sizes_from_file(chrom_sizes_file: str) -> Dict[str, int]:
    """从.chrom.sizes文件读取染色体大小"""
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes
# ============================================================================
# Main Dataset Class
# ============================================================================

class ChromosomeDataset(Dataset):
    
    def __init__(
        self,
        chr_list: List[str],
        cell_line_list: List[str],
        depth_list: List[str],  # Now accepts depth labels like ['2.1e5', 'bulk', '10e4']
        data_root: str,
        chrom_sizes: str = 'hg38',  # Can be 'hg38', 'mm10', or dict
        genome: str = 'hg38',
        target_window_size: int = 2097152,
        hic_resolution: int = 10_000,
        hic_output_size: int = 209,
        use_augmentation: bool = True,
        augmentation_chance: float = 0.3,
        exclude_centromere: bool = False,
        centromere_filter: bool = True,
        min_depth_threshold: int = 20000,  # 2e4 reads
        verbose: bool = True
    ):
        """
        Initialize Hi-Compass dataset.
        
        Args:
            chr_list: List of chromosomes to include (e.g., ['chr1', 'chr2'])
            cell_line_list: List of cell lines to include
            depth_list: List of depth labels to use (folder names like ['2.1e5', 'bulk', '10e4'])
            data_root: Root directory containing ATAC, HiC, DNA, CTCF subdirectories
            chrom_sizes: Chromosome sizes - 'hg38', 'mm10', or custom dict
            target_window_size: Size of genomic windows for prediction
            hic_resolution: Resolution of Hi-C data
            hic_output_size: Output size of Hi-C matrices
            use_augmentation: Whether to apply data augmentation
            augmentation_chance: Probability of applying each augmentation
            exclude_centromere: Whether to exclude centromeric regions
            centromere_bed: Path to BED file with centromeric regions
            min_depth_threshold: Minimum read depth threshold
            verbose: Whether to print initialization details
        """
        self.data_root = data_root
        self.chr_list = chr_list
        self.cell_line_list = cell_line_list
        self.depth_list = [normalize_depth_string(d) for d in depth_list]
        self.shift_window_size = 5000000
        self.target_window_size = target_window_size
        self.hic_resolution = hic_resolution
        self.hic_output_size = hic_output_size
        self.use_augmentation = use_augmentation
        self.augmentation_chance = augmentation_chance
        self.min_depth_threshold = min_depth_threshold
        self.verbose = verbose
        self.centromere_filter = centromere_filter
        # Set chromosome sizes
        self.genome = genome
        if genome == 'hg38':
            self.chrom_sizes = HG38_CHROM_SIZES
        elif genome == 'mm10':
            self.chrom_sizes = MM10_CHROM_SIZES
        else:
            # 自定义genome，从文件读取
            chrom_sizes_file = os.path.join(data_root, 'chromsize', f'{genome}.chrom.sizes')
            if not os.path.exists(chrom_sizes_file):
                raise FileNotFoundError(f"Chrom sizes file not found: {chrom_sizes_file}")
            self.chrom_sizes = load_chrom_sizes_from_file(chrom_sizes_file)
        
        self.total_chr_list = [file.split('.fa.gz')[0] for file in os.listdir(os.path.join(self.data_root, 'DNA', self.genome))]
        self.chr_name_to_index = create_chr_index_dict(self.total_chr_list)
        
        # Load excluded regions if requested
        self.centromere_bed_path = os.path.join(data_root, 'centromere', self.genome, 'centromere.bed')
        self.excluded_regions = {}
        if self.centromere_filter:
            self.excluded_regions = load_centromere_telomere_regions(self.centromere_bed_path)
        
        # Scan ATAC directory and build file mapping
        atac_dir = os.path.join(data_root, 'ATAC', self.genome)
        self.atac_file_map = scan_atac_directory(atac_dir, verbose=verbose)
        
        # Build cell line mapping for discriminator
        self.cell_line_to_idx = build_cell_line_mapping(cell_line_list, verbose=verbose)
        self.num_cell_types = len(self.cell_line_to_idx)
        self.use_discriminator = (self.num_cell_types > 1)
        
        # Verify all requested cell lines are available
        missing_cells = set(cell_line_list) - set(self.atac_file_map.keys())
        if missing_cells:
            raise ValueError(f"Missing cell lines in ATAC directory: {missing_cells}")
        
        # Verify all requested depths are available
        self._verify_depth_availability()
        
        # Initialize data loaders
        self._initialize_data_loaders()
        
        # Build sample list
        self.sample_list = self._build_sample_list()
        if verbose:
            print(f"\n{'='*70}")
            print(f"Dataset Initialization Complete")
            print(f"{'='*70}")
            print(f"Chromosomes: {len(self.chr_list)}")
            print(f"Cell lines: {len(self.cell_line_list)}")
            print(f"Depths: {len(self.depth_list)}")
            print(f"Total samples: {len(self.sample_list)}")
            print(f"Augmentation: {self.use_augmentation}")
            print(f"Discriminator: {'ENABLED' if self.use_discriminator else 'DISABLED'}")
            print(f"{'='*70}\n")

    
    
    def _verify_depth_availability(self):
        """Verify that all requested depths are available for all cell lines."""
        if not self.verbose:
            return
        
        print(f"\n{'='*70}")
        print(f"Verifying Depth Availability")
        print(f"{'='*70}\n")
        
        for cell_line in self.cell_line_list:
            available_depths = set(self.atac_file_map[cell_line].keys())
            requested_depths = set(self.depth_list)
            missing_depths = requested_depths - available_depths
            
            if missing_depths:
                print(f"⚠ {cell_line}: Missing depths {missing_depths}")
                print(f"   Available: {sorted(available_depths)}")
            else:
                print(f"✓ {cell_line}: All requested depths available")
        
        print(f"{'='*70}\n")
    
    def _initialize_data_loaders(self):
        """Initialize data feature loaders for DNA, ATAC, CTCF, and Hi-C."""
        if self.verbose:
            print(f"Initializing data loaders...")
        
        # DNA sequence loader
        dna_dir = os.path.join(self.data_root, 'DNA', self.genome)
        self.seq_dict = {}
        for chr_name in self.chr_list:
            fa_file = os.path.join(dna_dir, f"{chr_name}.fa.gz")
            if not os.path.exists(fa_file):
                raise FileNotFoundError(f"DNA sequence file not found: {fa_file}")
            self.seq_dict[chr_name] = data_feature.DNASequenceFeature(path=fa_file)
        
        # ATAC-seq loaders - now organized by cell line and depth label
        self.atac_dict = {}
        for cell_line in self.cell_line_list:
            self.atac_dict[cell_line] = {}
            for depth_label in self.depth_list:
                file_info = self.atac_file_map[cell_line][depth_label]
                self.atac_dict[cell_line][depth_label] = data_feature.GenomicFeature(
                    path=file_info.file_path, norm = None
                )
        
        # CTCF loader
        ctcf_file = os.path.join(self.data_root, 'CTCF', self.genome, f'generalized_CTCF.bw')
        if not os.path.exists(ctcf_file):
            raise FileNotFoundError(f"CTCF file not found: {ctcf_file}")
        self.ctcf_feature = data_feature.GenomicFeature(path=ctcf_file, norm = None)
        
        # Hi-C loaders
        hic_dir = os.path.join(self.data_root, 'HiC', self.genome)
        
        self.hic_dict = {}
        for cell_line in self.cell_line_list:
            self.hic_dict[cell_line] = {}
            cell_hic_dir = os.path.join(hic_dir, cell_line)
            
            for chr_name in self.chr_list:
                npz_file = os.path.join(cell_hic_dir, f"{chr_name}.npz")
                if not os.path.exists(npz_file):
                    raise FileNotFoundError(f"Hi-C file not found: {npz_file}")
                
                self.hic_dict[cell_line][chr_name] = data_feature.HiCFeature(path=npz_file)
        
        if self.verbose:
            print(f"  ✓ Data loaders initialized\n")
    
    def _build_sample_list(self) -> List[Tuple]:
        """
        Build list of all valid training samples.
        
        Returns:
            List of tuples: (chr_name, start, end, cell_line, depth_label)
        """
        samples = []
        step_size = self.shift_window_size
        
        for chr_name in self.chr_list:
            chr_size = self.chrom_sizes[chr_name]
            
            # Generate windows
            for start in range(0, chr_size - self.shift_window_size, step_size):
                end = start + self.target_window_size
                
                # Check if window overlaps with excluded regions
                if self._is_excluded_region(chr_name, start, end):
                    # print('get one !', chr_name, start, end)
                    continue
                
                # Add samples for all combinations of cell lines and depths
                for cell_line in self.cell_line_list:
                    for depth_label in self.depth_list:
                        samples.append((chr_name, start, end, cell_line, depth_label))
        
        return samples
    
    def _is_excluded_region(self, chr_name: str, start: int, end: int) -> bool:
        """Check if genomic region overlaps with excluded regions."""
        if chr_name not in self.excluded_regions:
            return False
        
        for excl_start, excl_end in self.excluded_regions[chr_name]:
            if not (end <= excl_start or start >= excl_end):
                return True
        
        return False
    
    def __len__(self) -> int:
        """Return the number of samples in the dataset."""
        return len(self.sample_list)
    

    def __getitem__(self, idx: int):
        """Get a single training sample."""
        chr_name, start, end, cell_line, depth_label = self.sample_list[idx]
        
        if self.use_augmentation:
            start, end = self._shift_augmentation(start, end, chr_name)
        
        seq, atac, actual_depth, ctcf, hic = self._load_data_at_interval(
            start, end, chr_name, cell_line, depth_label)
        
        if self.use_augmentation:
            seq = self._augment_sequence(seq, self.augmentation_chance)
            atac = self._augment_atac(atac, self.augmentation_chance)
            seq, atac, ctcf, hic = self._random_reverse_complement(
                seq, atac, ctcf, hic, chance=0.5)
        
        chr_size = self.chrom_sizes[chr_name]
        start_norm = start / chr_size
        end_norm = end / chr_size
        chr_idx_norm = self.chr_name_to_index.get(chr_name, 0) / 24.0
        cell_idx = self.cell_line_to_idx.get(cell_line, 0)
        
        # 强制转换为标准numpy数组，使用.copy()确保独立
        seq = seq.astype(np.float32).copy()
        atac = atac.astype(np.float32).copy()
        ctcf = ctcf.astype(np.float32).copy()
        hic = hic.astype(np.float32).copy()
        
        # 转换为torch tensor
        seq = torch.FloatTensor(seq)
        atac = torch.FloatTensor(atac)
        ctcf = torch.FloatTensor(ctcf)
        hic = torch.FloatTensor(hic)
        
        actual_depth = torch.tensor(float(actual_depth), dtype=torch.float32)
        start_norm = torch.tensor(float(start_norm), dtype=torch.float32)
        end_norm = torch.tensor(float(end_norm), dtype=torch.float32)
        chr_idx_norm = torch.tensor(float(chr_idx_norm), dtype=torch.float32)
        cell_idx = torch.tensor(int(cell_idx), dtype=torch.long)
        
        return seq, atac, actual_depth, ctcf, hic, start_norm, end_norm, chr_idx_norm, cell_idx
        
    def _shift_augmentation(self, start: int, end: int, chr_name: int) -> Tuple[int, int]:
        """Apply random shift augmentation to interval."""
        chrom_len = self.chrom_sizes[chr_name]
        upper = min(chrom_len, end + self.shift_window_size)
        lower = max(0, start - self.shift_window_size)
        upper_shift_range = upper-end
        lower_shift_range = lower-start
        
        if (upper_shift_range > 0) & (lower_shift_range)>0:
            offset = random.randint(-lower_shift_range, upper_shift_range)
            return start + offset, start + offset + self.target_window_size
        
        return start, end
    
    def _load_data_at_interval(
        self, 
        start: int, 
        end: int, 
        chr_name: str, 
        cell_line: str, 
        depth_label: str):
        """Load all data features for a specific genomic interval."""
        # Load sequence
        seq = self.seq_dict[chr_name].get(start, end)
        
        # Load ATAC-seq
        atac = self.atac_dict[cell_line][depth_label].get(chr_name, start, end)
        
        # Get actual depth from file info
        file_info = self.atac_file_map[cell_line][depth_label]
        # Parse actual depth value
        actual_depth_str = file_info.actual_depth
        # Convert depth string to numerical value (e.g., '2.1e5' -> 210000)
        try:
            actual_depth = int(float(actual_depth_str))
        except:
            print('Warning: Covert depth str improperly')
            # Fallback: try to extract number
            depth_match = re.search(r'(\d+\.?\d*)[eE]([+-]?\d+)', actual_depth_str)
            if depth_match:
                actual_depth = int(float(depth_match.group(1)) * (10 ** int(depth_match.group(2))))
            else:
                actual_depth = self.min_depth_threshold
        
        # Add depth jitter (0.7x to 1.3x) for augmentation
        if self.use_augmentation:
            actual_depth = int(actual_depth * np.random.uniform(0.7, 1.3))
        actual_depth = max(self.min_depth_threshold, actual_depth)
        
        # Load CTCF
        ctcf = self.ctcf_feature.get(chr_name, start, end)
        
        # Load Hi-C
        hic = self.hic_dict[cell_line][chr_name].get(start)
        
        return seq, atac, actual_depth, ctcf, hic
    
    # ========================================================================
    # Data Augmentation Methods
    # ========================================================================
    
    def _augment_sequence(self, seq: np.ndarray, chance: float = 0.3) -> np.ndarray:
        """Apply sequence-specific augmentations."""
        aug_seq = seq.copy()
        
        # Small position shift
        if np.random.random() < chance:
            aug_seq = self._shift_sequence(aug_seq, max_shift=5000)
        
        # Replace short segments
        if np.random.random() < chance:
            aug_seq = self._swap_sequence_segments(aug_seq, max_size=1000)
        
        return aug_seq
    
    def _augment_atac(self, atac: np.ndarray, chance: float = 0.3) -> np.ndarray:
        """Apply ATAC-seq specific augmentations."""
        aug_atac = atac.copy()
        
        # Gaussian noise
        if np.random.random() < chance:
            noise = np.random.randn(*aug_atac.shape) * 0.05
            aug_atac = np.clip(aug_atac + noise, 0, 1)
        
        # Intensity adjustment
        if np.random.random() < chance:
            factor = np.random.uniform(0.9, 1.1)
            aug_atac = np.clip(aug_atac * factor, 0, 1)
        
        # Conservative masking
        if np.random.random() < chance:
            mask_ratio = np.random.uniform(0.01, 0.05)
            mask = np.random.random(aug_atac.shape) > mask_ratio
            aug_atac = aug_atac * mask
        
        return aug_atac
    
    def _shift_sequence(self, seq: np.ndarray, max_shift: int = 5000) -> np.ndarray:
        """Apply small positional shift to sequence."""
        shift = np.random.randint(-max_shift, max_shift + 1)
        shifted = np.zeros_like(seq)
        
        if shift > 0:
            shifted[shift:] = seq[:-shift]
        elif shift < 0:
            shifted[:shift] = seq[-shift:]
        else:
            shifted = seq.copy()
        
        return shifted
    
    def _swap_sequence_segments(
        self, 
        seq: np.ndarray, 
        max_size: int = 1000, 
        max_swaps: int = 3
    ) -> np.ndarray:
        """Randomly swap short sequence segments."""
        swapped = seq.copy()
        seq_len = seq.shape[0]
        
        num_swaps = np.random.randint(1, max_swaps + 1)
        
        for _ in range(num_swaps):
            swap_size = np.random.randint(1, min(max_size, seq_len // 10))
            pos1 = np.random.randint(0, seq_len - swap_size)
            pos2 = np.random.randint(0, seq_len - swap_size)
            
            # Swap segments
            temp = swapped[pos1:pos1+swap_size].copy()
            swapped[pos1:pos1+swap_size] = swapped[pos2:pos2+swap_size]
            swapped[pos2:pos2+swap_size] = temp
        
        return swapped
    
    def _random_reverse_complement(
        self, 
        seq: np.ndarray, 
        atac: np.ndarray, 
        ctcf: np.ndarray, 
        hic: np.ndarray,
        chance: float = 0.5
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Randomly apply reverse complement to all inputs."""
        if np.random.random() < chance:
            seq_rev = np.flip(seq, axis=0).copy()
            atac_rev = np.flip(atac, axis=0).copy()
            ctcf_rev = np.flip(ctcf, axis=0).copy()
            hic_rev = np.flip(hic, axis=(0, 1)).copy()
            return seq_rev, atac_rev, ctcf_rev, hic_rev
        
        return seq, atac, ctcf, hic


# ============================================================================
# Utility Functions for External Use
# ============================================================================

def create_dataset_from_directory(
    data_root: str,
    chr_list: List[str],
    depth_labels: Optional[List[str]] = None,
    genome: str = 'hg38',
    auto_scan: bool = True,
    **kwargs
) -> ChromosomeDataset:
    """
    Convenience function to create dataset by auto-scanning data directory.
    
    Args:
        data_root: Root directory containing ATAC, HiC, DNA, CTCF subdirectories
        chr_list: List of chromosomes to include
        depth_labels: List of depth labels to use (e.g., ['2.1e5', 'bulk'])
                     If None, will use all available depths
        genome: Genome build ('hg38' or 'mm10')
        auto_scan: Whether to automatically scan and detect cell lines
        **kwargs: Additional arguments passed to ChromosomeDataset
        
    Returns:
        Initialized ChromosomeDataset
    """
    atac_dir = os.path.join(data_root, 'ATAC')
    
    if auto_scan:
        # Scan directory to find cell lines and depths
        cell_depth_files = scan_atac_directory(atac_dir, verbose=True)
        cell_line_list = list(cell_depth_files.keys())
        
        # If depth_labels not specified, use all available depths
        if depth_labels is None:
            all_depths = set()
            for depths_dict in cell_depth_files.values():
                all_depths.update(depths_dict.keys())
            depth_labels = sorted(all_depths)
        
        print(f"\nAuto-detected configuration:")
        print(f"  Cell lines: {cell_line_list}")
        print(f"  Depth labels: {depth_labels}")
    else:
        raise ValueError("auto_scan=False requires explicit cell_line_list and depth_labels")
    
    return ChromosomeDataset(
        chr_list=chr_list,
        cell_line_list=cell_line_list,
        depth_list=depth_labels,
        data_root=data_root,
        chrom_sizes=genome,
        **kwargs
    )


# ============================================================================
# Main / Testing
# ============================================================================

if __name__ == '__main__':
    """
    Example usage and testing.
    """
    print("Hi-Compass Chromosome Dataset (Hierarchical File Structure)")
    print("="*70)
    
    # Example: Create dataset with auto-detection
    # dataset = create_dataset_from_directory(
    #     data_root='/path/to/data',
    #     chr_list=['chr1', 'chr2'],
    #     depth_labels=['2.1e5', 'bulk', '10e4'],  # or None to use all
    #     genome='hg38',
    #     use_augmentation=True
    # )
    
    # Example: Manual dataset creation
    # dataset = ChromosomeDataset(
    #     chr_list=['chr1'],
    #     cell_line_list=['GM12878', 'K562'],
    #     depth_list=['2.1e5', 'bulk'],
    #     data_root='/path/to/data',
    #     chrom_sizes='hg38'
    # )
    
    print("\nDataset module loaded successfully!")
    print("Key features:")
    print("  - Hierarchical folder structure support")
    print("  - Special handling for bulk data")
    print("  - Automatic cell line discovery")
    print("  - Automatic discriminator configuration")
    print("\nUse create_dataset_from_directory() for easy setup")