import numpy as np
import os
import pandas as pd
from functools import reduce
import time
import cooler
chrome_size_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
 'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 
 'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 
 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 
 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983, 
 'chr22': 50818468}


def positive(matrix):
    min_value = matrix.min()
    if min_value>=0:
        return matrix
    else:
        return matrix-min_value

def norm_hic(input_matrix, max_num=6):
    now_max = max(np.diagonal(input_matrix))
    return input_matrix/((now_max+0.0001)/(max_num+0.0001))


def make_cooler_bins_chr(chr_name, length, resolution):
    total_bin = []
    for i in range(0, length, resolution):
        total_bin.append([chr_name, i, min(length, i+resolution)])
    df = pd.DataFrame(total_bin, columns=['chrom', 'start', 'end'])
    return df

def get_chr_stack(chr_list, chr_name, chrome_size_dict=chrome_size_dict, resolution=10000):
    chrome_size_dict = {key: chrome_size_dict[key] for key in chr_list if key in chrome_size_dict}
    chr_before_list = chr_list[:chr_list.index(chr_name)]
    if len(chr_before_list) == 0:
        return 0
    else:
        stack = 0
        for before_chr in chr_before_list:
            stack+= int(chrome_size_dict[before_chr]/resolution) + 1
        return stack

def merge_hic_segment(hic_list, save_path, window_size=2097152, resolution=10000):
    chr_list = [i for i in hic_list]
    bins = int(window_size/resolution)
    chr_hic_dict = {}
    for chr_num in chr_list:
        chr_hic_dict[chr_num] = [make_cooler_bins_chr(chr_name=chr_num, length=chrome_size_dict[chr_num], resolution=resolution)]
        sub_list = hic_list[chr_num]
        # print(sub_list)
        large_pic = np.zeros((int(chrome_size_dict[chr_num]/resolution), 
                              int(chrome_size_dict[chr_num]/resolution)))
        for segement in sub_list:
            # print(segement[0])
            sub_start_bin = int(segement[0]/resolution)
            large_pic[sub_start_bin: sub_start_bin+bins, sub_start_bin: sub_start_bin+bins] = segement[2]
        # 将单chr的矩阵稀疏化
        rows, cols = np.nonzero(large_pic)
        stack = get_chr_stack(chr_list=chr_list, chr_name=chr_num, chrome_size_dict=chrome_size_dict, resolution=resolution)
        rows += stack
        cols += stack
        counts = large_pic[np.nonzero(large_pic)]
        large_pic_sp = np.column_stack((rows, cols, counts))
        large_pic_sp = pd.DataFrame(large_pic_sp)
        chr_hic_dict[chr_num].append(large_pic_sp)
    # 每个key对应的是一个list, [bin表, 单chr的hic三列矩阵]
    total_bin_list = []
    total_sp_list = []
    start_bin = 0
    total_bin_len = 0
    for chr_num in chr_hic_dict:
        total_bin_list.append(chr_hic_dict[chr_num][0])
        total_sp_list.append(chr_hic_dict[chr_num][1])
        # total_bin_len += chr_hic_dict[chr_num][1].shape[0]
    bin_df = reduce(lambda x, y: pd.concat([x,y], axis = 0), total_bin_list)
    bin_df.reset_index()
    bin_df.columns = ['chrom', 'start', 'end']
    sp_df = reduce(lambda x, y: pd.concat([x,y], axis = 0), total_sp_list)
    sp_df.reset_index()
    sp_df.columns = ['bin1_id', 'bin2_id', 'count']
    cooler.create_cooler(cool_uri=save_path, 
                         bins=bin_df, pixels=sp_df)

def cal_depth(filename):
    with open(filename, 'r') as file:
        depth = 0
        for _ in file:
            depth += 1
    return depth

def cal_size(filename):
    stats = os.stat(filename)
    return stats.st_size  


def pseudo_weight(mcool_path, save_path, weight=1):
    """
    mcool_path: path to mcool for processing
    res: resolution
    save_path: path to saved cool file
    weight: the pseudo weight
    """
    c = cooler.Cooler(f'{mcool_path}')
    bins_df = c.bins()[:]
    bins_df['weight'] = weight
    cooler.create_cooler(cool_uri=f'{save_path}', bins=bins_df, pixels=c.pixels()[:])
    cooler_obj = cooler.Cooler(f'{save_path}')
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