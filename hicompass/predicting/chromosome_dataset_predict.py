import sys 
import os
import random
import pickle
import pandas as pd
import numpy as np

from skimage.transform import resize
from torch.utils.data import Dataset

from . import data_feature


class ChromosomeDataset(Dataset):
    '''
    Dataloader that provide sequence, features, and HiC data. Assume input
    folder strcuture.
    Args:
        data_root (str): Directory including sequence features, and HiC matrix 
            as subdirectories.
        chr_name (str): Name of the represented chromosome (e.g. chr1)
            as ``root/DNA/chr1/DNA`` for DNA as an example.
        omit_regions (list of tuples): start and end of excluded regions
    '''
    def __init__(self, chr_name_list, atac_path_list, depth,
        general_ctcf_bw_path = '/cluster/home/Yuanchen/project/scHiC/dataset/CTCF/ctcf_count_norm.bw',
        dna_dir_path='/cluster/home/Yuanchen/project/scHiC/dataset/DNA/hg38', 
        omit_regions_file_path='/cluster/home/Yuanchen/project/scHiC/dataset/centrotelo.bed',
        use_aug=False, stride=158):
        self.use_aug = use_aug
        self.res = 10000 # 10kb resolution
        self.bins = 209.7152 # 209.7152 bins 2097152 bp
        self.sample_bins = 500
        self.atac_path_list = atac_path_list
        self.dna_dir_path = dna_dir_path
        self.general_ctcf_bw_path = general_ctcf_bw_path
        self.stride = stride 
        self.depth = depth
        self.chr_name_list = chr_name_list
        self.omit_regions_dict = proc_centrotelo(omit_regions_file_path)
        print(f'Loading chromosome {self.chr_name_list}...')
        self.seq_dict = {}
        for chr_name in self.chr_name_list:
            self.seq_dict[chr_name] = data_feature.DNASequenceFeature(path=f'{self.dna_dir_path}/{chr_name}.fa.gz')
        self.atac_list = [data_feature.GenomicFeature(path=path, norm=None) for path in self.atac_path_list]
        self.ctcf_feature = data_feature.GenomicFeature(path=self.general_ctcf_bw_path, norm= None)
        self.all_intervals_list = []
        for chr_name in self.chr_name_list:
            self.all_intervals_list.append(self.get_intervals_chr_ct_dep(seq=self.seq_dict[chr_name],
                                                                            chr_name=chr_name,
                                                                            omit_regions=self.omit_regions_dict[chr_name]))

        self.intervals = np.concatenate(self.all_intervals_list, axis=0)


    def __getitem__(self, idx):
        start, end, chr_name = self.intervals[idx]
        # start, end, chr_name = self.intervals[idx]
        start = int(start)
        end = int(end)
        depth = int(self.depth)
        target_size = int(self.bins * self.res)
        total_seq = self.seq_dict[chr_name]
        total_atac_list = self.atac_list
        total_ctcf = self.ctcf_feature
        start, end = self.shift_fix(target_size, start)
        start = int(start)
        end = int(end)
        seq, atac, real_depth, ctcf = self.get_data_at_chr_interval(start=start, end=end, chr_name=chr_name,
                                                                total_seq=total_seq,
                                                                atac_list=total_atac_list,
                                                                ctcf=total_ctcf,
                                                                depth=depth)
        return seq, atac, real_depth, ctcf, start, end, chr_name

    def __len__(self):
        return len(self.intervals)


    def get_data_at_chr_interval(self, start, end, chr_name, total_seq, atac_list, ctcf, depth):
        '''
        used in get item
        '''
        # Sequence processing
        start = int(start)
        end = int(end)
        # print('get_data_at_chr_interval',start, end)
        seq = total_seq.get(start, end)
        # Features processing
        features_atac_list = [item.get(chr_name, start, end) for item in atac_list]
        features_atac = np.sum(np.array(features_atac_list), axis=0)
        # features_atac = atac_list[0].get(chr_name, start, end)
        features_ctcf = ctcf.get(chr_name, start, end)
        real_depth = int(depth)
        real_depth = min(len(features_ctcf), real_depth)
        return seq, features_atac, real_depth, features_ctcf

    
    def get_intervals_chr_ct_dep(self, seq, chr_name, omit_regions):
        # chr_num = int(chr_name.split('chr')[-1])
        chr_bins = len(seq) / self.res
        data_size = (chr_bins - self.sample_bins) / self.stride
        starts = np.arange(0, data_size).reshape(-1, 1) * self.stride
        intervals_bin = np.append(starts, starts + self.sample_bins, axis=1)
        intervals = intervals_bin * self.res
        intervals = intervals.astype(int)
        intervals = self.filter(intervals, omit_regions)
        chr_name_array = np.full(len(intervals), chr_name).reshape(-1, 1)
        intervals = np.append(intervals, chr_name_array, axis=1)
        return intervals
    

    def get_active_intervals(self):
        '''
        Get intervals for sample data: [[start, end]]
        '''
        chr_bins = len(self.seq) / self.res
        data_size = (chr_bins - self.sample_bins) / self.stride
        starts = np.arange(0, data_size).reshape(-1, 1) * self.stride
        intervals_bin = np.append(starts, starts + self.sample_bins, axis=1)
        intervals = intervals_bin * self.res
        return intervals.astype(int)

    def filter(self, intervals, omit_regions):
        valid_intervals = []
        for start, end in intervals: 
            start_cond = start <= omit_regions[:, 1]
            end_cond = omit_regions[:, 0] <= end
            if sum(start_cond * end_cond) == 0:
                valid_intervals.append([start, end])
        return valid_intervals

    def encode_seq(self, seq):
        ''' 
        encode dna to onehot (n x 5)
        '''
        seq_emb = np.zeros((len(seq), 5))
        seq_emb[np.arange(len(seq)), seq] = 1
        return seq_emb


    def shift_fix(self, target_size, start):
        offset = 0
        return start + offset , start + offset + target_size
    
def proc_centrotelo(bed_dir):
    ''' Take a bed file indicating location, output a dictionary of items 
    by chromosome which contains a list of 2 value lists (range of loc)
    '''
    df = pd.read_csv(bed_dir , sep='\t', names=['chr', 'start', 'end'])
    chrs = df['chr'].unique()
    centrotelo_dict = {}
    for chr_name in chrs:
        sub_df = df[df['chr'] == chr_name]
        regions = sub_df.drop('chr', axis = 1).to_numpy()
        centrotelo_dict[chr_name] = regions
    return centrotelo_dict




