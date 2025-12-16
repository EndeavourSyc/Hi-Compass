# PredictDataset.py
import sys 
import os
import random
import pickle
import pandas as pd
import numpy as np

from skimage.transform import resize
from torch.utils.data import Dataset

from . import data_feature


# ============================================================================
# Chromosome Size Dictionaries
# ============================================================================

HG38_CHROM_SIZES = {
    'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
    'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 
    'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 
    'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415
}

MM10_CHROM_SIZES = {
    'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
    'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
    'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
    'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
    'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
    'chrY': 91744698
}


def load_chrom_sizes_from_file(chrom_sizes_file: str):
    """Load chromosome sizes from .chrom.sizes file"""
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def proc_centrotelo(bed_dir):
    """Load centromere/telomere regions from BED file"""
    if bed_dir is None or not os.path.exists(bed_dir):
        return {}
    
    df = pd.read_csv(bed_dir, sep='\t', names=['chr', 'start', 'end'])
    chrs = df['chr'].unique()
    centrotelo_dict = {}
    for chr_name in chrs:
        sub_df = df[df['chr'] == chr_name]
        regions = sub_df.drop('chr', axis=1).to_numpy()
        centrotelo_dict[chr_name] = regions
    return centrotelo_dict


class PredictDataset(Dataset):
    """
    Dataset for Hi-C prediction from genomic features.
    
    Args:
        chr_name_list: List of chromosome names to predict
        atac_path_list: List of ATAC-seq bigWig file paths
        depth: Sequencing depth for normalization
        dna_dir_path: Directory containing DNA sequence files
        general_ctcf_bw_path: Path to CTCF bigWig file
        omit_regions_file_path: Path to BED file with excluded regions (optional)
        genome: Genome build ('hg38', 'mm10', or custom)
        chrom_sizes_file: Path to .chrom.sizes file (required if genome not hg38/mm10)
        stride: Stride for sliding windows (in bins)
    """
    
    def __init__(
        self, 
        chr_name_list,
        atac_path_list,
        depth,
        dna_dir_path,
        general_ctcf_bw_path,
        omit_regions_file_path=None,
        genome='hg38',
        chrom_sizes_file=None,
        stride=158
    ):
        self.res = 10000
        self.bins = 209.7152
        self.sample_bins = 500
        self.target_window_size = 2097152
        self.atac_path_list = atac_path_list
        self.dna_dir_path = dna_dir_path
        self.general_ctcf_bw_path = general_ctcf_bw_path
        self.stride = stride 
        self.depth = depth
        self.chr_name_list = chr_name_list
        self.genome = genome
        
        # Load chromosome sizes
        if genome == 'hg38':
            self.chrom_sizes = HG38_CHROM_SIZES
        elif genome == 'mm10':
            self.chrom_sizes = MM10_CHROM_SIZES
        else:
            if chrom_sizes_file is None:
                raise ValueError(
                    f"Custom genome '{genome}' requires chrom_sizes_file parameter"
                )
            if not os.path.exists(chrom_sizes_file):
                raise FileNotFoundError(f"Chromosome sizes file not found: {chrom_sizes_file}")
            self.chrom_sizes = load_chrom_sizes_from_file(chrom_sizes_file)
        
        # Load excluded regions (optional)
        self.omit_regions_dict = proc_centrotelo(omit_regions_file_path)
        
        # Load genomic features
        self.seq_dict = {}
        for chr_name in self.chr_name_list:
            fa_file = os.path.join(self.dna_dir_path, f'{chr_name}.fa.gz')
            self.seq_dict[chr_name] = data_feature.DNASequenceFeature(path=fa_file)
        
        self.atac_list = [
            data_feature.GenomicFeature(path=path, norm=None) 
            for path in self.atac_path_list
        ]
        self.ctcf_feature = data_feature.GenomicFeature(
            path=self.general_ctcf_bw_path, 
            norm=None
        )
        
        # Build intervals
        self.all_intervals_list = []
        for chr_name in self.chr_name_list:
            intervals = self.get_intervals_chr(
                chr_name=chr_name,
                chr_size=self.chrom_sizes[chr_name],
                omit_regions=self.omit_regions_dict.get(chr_name, np.array([]))
            )
            self.all_intervals_list.append(intervals)
        
        self.intervals = np.concatenate(self.all_intervals_list, axis=0)

    def __len__(self):
        return len(self.intervals)

    def __getitem__(self, idx):
        start, end, chr_name = self.intervals[idx]
        start = int(start)
        end = int(end)
        
        # Fix window size to target_window_size
        start, end = self.shift_fix(self.target_window_size, start)
        start = int(start)
        end = int(end)
        
        seq = self.seq_dict[chr_name].get(start, end)
        
        features_atac_list = [
            atac.get(chr_name, start, end) 
            for atac in self.atac_list
        ]
        atac = np.sum(np.array(features_atac_list), axis=0)
        
        ctcf = self.ctcf_feature.get(chr_name, start, end)
        
        real_depth = int(self.depth)
        real_depth = min(len(ctcf), real_depth)
        
        return seq, atac, real_depth, ctcf, start, end, chr_name

    def shift_fix(self, target_size, start):
        """Fix window size to target_size"""
        offset = 0
        return start + offset, start + offset + target_size

    def get_intervals_chr(self, chr_name, chr_size, omit_regions):
        """Generate prediction intervals for a chromosome"""
        chr_bins = chr_size / self.res
        data_size = (chr_bins - self.sample_bins) / self.stride
        starts = np.arange(0, data_size).reshape(-1, 1) * self.stride
        intervals_bin = np.append(starts, starts + self.sample_bins, axis=1)
        intervals = intervals_bin * self.res
        intervals = intervals.astype(int)
        
        # Filter only if omit_regions is provided
        if len(omit_regions) > 0:
            intervals = self.filter(intervals, omit_regions)
        
        chr_name_array = np.full(len(intervals), chr_name).reshape(-1, 1)
        intervals = np.append(intervals, chr_name_array, axis=1)
        return intervals

    def filter(self, intervals, omit_regions):
        """Filter out intervals overlapping with excluded regions"""
        valid_intervals = []
        for start, end in intervals: 
            start_cond = start <= omit_regions[:, 1]
            end_cond = omit_regions[:, 0] <= end
            if sum(start_cond * end_cond) == 0:
                valid_intervals.append([start, end])
        return valid_intervals