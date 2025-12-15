# utils.py
import numpy as np
import os
import pandas as pd
from functools import reduce
import time
import cooler


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


def load_chrom_sizes_from_file(chrom_sizes_file):
    """Load chromosome sizes from .chrom.sizes file"""
    chrom_sizes = {}
    with open(chrom_sizes_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom_sizes[parts[0]] = int(parts[1])
    return chrom_sizes


def get_chrom_sizes(genome, chrom_sizes_file=None):
    """
    Get chromosome sizes for a genome build.
    
    Args:
        genome: Genome build ('hg38', 'mm10', or custom)
        chrom_sizes_file: Path to .chrom.sizes file (required if genome not hg38/mm10)
        
    Returns:
        Dictionary mapping chromosome names to sizes
    """
    if genome == 'hg38':
        return HG38_CHROM_SIZES
    elif genome == 'mm10':
        return MM10_CHROM_SIZES
    else:
        if chrom_sizes_file is None:
            raise ValueError(
                f"Custom genome '{genome}' requires chrom_sizes_file parameter"
            )
        if not os.path.exists(chrom_sizes_file):
            raise FileNotFoundError(f"Chromosome sizes file not found: {chrom_sizes_file}")
        return load_chrom_sizes_from_file(chrom_sizes_file)


def positive(matrix):
    min_value = matrix.min()
    if min_value >= 0:
        return matrix
    else:
        return matrix - min_value


def norm_hic(input_matrix, max_num=6):
    now_max = max(np.diagonal(input_matrix))
    return input_matrix / ((now_max + 0.0001) / (max_num + 0.0001))


def make_cooler_bins_chr(chr_name, length, resolution):
    total_bin = []
    for i in range(0, length, resolution):
        total_bin.append([chr_name, i, min(length, i + resolution)])
    df = pd.DataFrame(total_bin, columns=['chrom', 'start', 'end'])
    return df


def get_chr_stack(chr_list, chr_name, chrom_sizes, resolution=10000):
    """
    Calculate the bin offset for a chromosome in the genome-wide matrix.
    
    Args:
        chr_list: List of all chromosomes
        chr_name: Current chromosome name
        chrom_sizes: Dictionary of chromosome sizes
        resolution: Hi-C resolution
        
    Returns:
        Number of bins before this chromosome
    """
    chr_before_list = chr_list[:chr_list.index(chr_name)]
    if len(chr_before_list) == 0:
        return 0
    else:
        stack = 0
        for before_chr in chr_before_list:
            stack += int(chrom_sizes[before_chr] / resolution) + 1
        return stack


def merge_hic_segment(hic_list, save_path, chrom_sizes, window_size=2097152, resolution=10000):
    """
    Merge Hi-C segments into a single cooler file.
    
    Args:
        hic_list: Dictionary mapping chromosome names to lists of segments
        save_path: Output path for cooler file
        chrom_sizes: Dictionary of chromosome sizes
        window_size: Prediction window size in bp
        resolution: Hi-C resolution in bp
    """
    chr_list = list(hic_list.keys())
    bins = int(window_size / resolution)
    chr_hic_dict = {}
    
    for chr_num in chr_list:
        chr_hic_dict[chr_num] = [
            make_cooler_bins_chr(
                chr_name=chr_num,
                length=chrom_sizes[chr_num],
                resolution=resolution
            )
        ]
        sub_list = hic_list[chr_num]
        large_pic = np.zeros((
            int(chrom_sizes[chr_num] / resolution), 
            int(chrom_sizes[chr_num] / resolution)
        ))
        
        for segment in sub_list:
            sub_start_bin = int(segment[0] / resolution)
            large_pic[sub_start_bin:sub_start_bin + bins, 
                     sub_start_bin:sub_start_bin + bins] = segment[2]
        
        # Sparsify matrix
        rows, cols = np.nonzero(large_pic)
        stack = get_chr_stack(
            chr_list=chr_list,
            chr_name=chr_num,
            chrom_sizes=chrom_sizes,
            resolution=resolution
        )
        rows += stack
        cols += stack
        counts = large_pic[np.nonzero(large_pic)]
        large_pic_sp = np.column_stack((rows, cols, counts))
        large_pic_sp = pd.DataFrame(large_pic_sp)
        chr_hic_dict[chr_num].append(large_pic_sp)
    
    # Merge all chromosomes
    total_bin_list = []
    total_sp_list = []
    
    for chr_num in chr_hic_dict:
        total_bin_list.append(chr_hic_dict[chr_num][0])
        total_sp_list.append(chr_hic_dict[chr_num][1])
    
    bin_df = reduce(lambda x, y: pd.concat([x, y], axis=0), total_bin_list)
    bin_df.reset_index(drop=True, inplace=True)
    bin_df.columns = ['chrom', 'start', 'end']
    
    sp_df = reduce(lambda x, y: pd.concat([x, y], axis=0), total_sp_list)
    sp_df.reset_index(drop=True, inplace=True)
    sp_df.columns = ['bin1_id', 'bin2_id', 'count']
    
    cooler.create_cooler(cool_uri=save_path, bins=bin_df, pixels=sp_df)


def pseudo_weight(mcool_path, save_path, weight=1):
    """
    Add pseudo weight to cooler file for visualization.
    
    Args:
        mcool_path: Input cooler file path
        save_path: Output cooler file path
        weight: Pseudo weight value
    """
    c = cooler.Cooler(mcool_path)
    bins_df = c.bins()[:]
    bins_df['weight'] = weight
    cooler.create_cooler(cool_uri=save_path, bins=bins_df, pixels=c.pixels()[:])
    cooler_obj = cooler.Cooler(save_path)
    stats = {
        "min_nnz": 0,
        "min_count": 0,
        "mad_max": 0,
        "cis_only": False,
        "ignore_diags": 2,
        "converged": False,
        "divisive_weights": False,
    }
    with cooler_obj.open("r+") as cljr:
        cljr["bins"]['weight'].attrs.update(stats)