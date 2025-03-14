import sys 
import os
import random
import pickle
import pandas as pd
import numpy as np

from skimage.transform import resize
from torch.utils.data import Dataset

import data_feature as data_feature

cell_line_class = {'a549': 0, 'imr90': 1, 'h1_hesc': 2, 'adipocyte': 3, 
              'transverse_colon': 4, 'k562': 5, 'hffc6': 6, 'gm12878': 7, 
              'sun449': 8,'retina_1': 9,'retina_2':9, 'retina_3': 9, 'retina_4' : 9, 
              'gastrocnemius': 10, 'osteogenesis':11, 'neuron':12, 'huh7':13, 
              'huh1': 14, 'hmec': 15, 'hep3b': 16, 'dnd41': 17, 
              'ct27_stem' : 18, 'ct27_evt': 18, 
              'colon_5':19, 'colon_4':19, 'colon_3':19, 'colon_2':19, 'colon_1':19, 
              'adipogenesis':20, 'foresking_fibroblasts': 21, 'Tcell': 22}
chr_name_class = {'chr1':1, 'chr2': 2, 'chr3': 3, 'chr4': 4,
                  'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9,
                  'chr10': 10, 'chr11': 11, 'chr12': 12, 'chr13': 13, 'chr14': 14,
                  'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19,
                  'chr20': 20, 'chr21': 21, 'chr22': 22, 'chrX': 23, 'chrY': 24}
chrome_size_dict = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
 'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 
 'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 
 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 
 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983, 
 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}

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
    def __init__(self, chr_name_list, atac_path_list, depth, chrom_max_num=24, omit_regions_file_path='/cluster/home/Yuanchen/project/scHiC/dataset/centrotelo.bed',use_aug=False, stride=158):
        self.use_aug = use_aug
        self.res = 10000 # 10kb resolution
        self.bins = 209.7152 # 209.7152 bins 2097152 bp
        self.sample_bins = 500
        self.chrom_max_num = chrom_max_num
        self.atac_path_list = atac_path_list
        self.stride = stride # bins
        self.depth = depth
        self.chr_name_list = chr_name_list
        self.add_item = '' #_1k_mean0_norm
        self.omit_regions_dict = proc_centrotelo(omit_regions_file_path)
        print(f'Loading chromosome {self.chr_name_list}...')
        self.seq_dict = {}
        for chr_name in self.chr_name_list:
            self.seq_dict[chr_name] = data_feature.SequenceFeature(path=f'/cluster/home/Yuanchen/project/scHiC/dataset/DNA/hg38/{chr_name}.fa.gz')
        
        self.atac_list = [data_feature.GenomicFeature(path=path, norm=None) for path in self.atac_path_list]
        self.ctcf_feature = data_feature.GenomicFeature(path=f'/cluster/home/Yuanchen/project/scHiC/dataset/CTCF/ctcf_count_norm.bw', norm= None)
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
        # total_mat_ref = self.mat_ref_dict[chr_name]
        total_atac_list = self.atac_list
        total_ctcf = self.ctcf_feature
        # Shift Augmentations
        if self.use_aug:
            start, end = self.shift_aug(target_size, start, end)
            start = int(start)
            end = int(end)
        else:
            start, end = self.shift_fix(target_size, start)
            start = int(start)
            end = int(end)
        # print(start)
        seq, atac, real_depth, ctcf = self.get_data_at_chr_interval(start=start, end=end, chr_name=chr_name,
                                                                total_seq=total_seq,
                                                                atac_list=total_atac_list,
                                                                ctcf=total_ctcf,
                                                                depth=depth)
        # print('get_item', start, end)
        return seq, atac, real_depth, ctcf, start, start/chrome_size_dict[chr_name], end, end/chrome_size_dict[chr_name], chr_name, chr_name_class[chr_name]/self.chrom_max_num

    def __len__(self):
        return len(self.intervals)

    def gaussian_noise(self, inputs, std = 1, chance = 0.5):
        prob = np.random.uniform()
        if prob <= chance:
            noise = np.random.randn(*inputs.shape) * std
            outputs = inputs + noise
            return outputs
        else:
            return inputs

    def reverse(self, seq, atac, ctcf, mat, chance = 0.5):
        '''
        Reverse sequence and matrix
        '''
        r_bool = np.random.rand(1)
        if r_bool < chance:
            seq_r = np.flip(seq, 0).copy() # n x 5 shape
            atac_r = np.flip(atac, 0).copy()
            ctcf_r = np.flip(ctcf, 0).copy()
            mat_r = np.flip(mat, [0, 1]).copy() # n x n

            # Complementary sequence
            # seq_r = self.complement(seq_r)


        else:
            seq_r = seq
            atac_r = atac
            ctcf_r = ctcf
            mat_r = mat
        return seq_r, atac_r, ctcf_r, mat_r

    def complement(self, seq, chance = 0.5):
        '''
        Complimentary sequence
        '''
        r_bool = np.random.rand(1)
        if r_bool < chance:
            seq_comp = np.concatenate([seq[:, 1:2],
                                       seq[:, 0:1],
                                       seq[:, 3:4],
                                       seq[:, 2:3],
                                       seq[:, 4:5]], axis = 1)
        else:
            seq_comp = seq
        return seq_comp

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
            # Way smaller than omit or way larger than omit
            start_cond = start <= omit_regions[:, 1]
            end_cond = omit_regions[:, 0] <= end
            #import pdb; pdb.set_trace()
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

    def shift_aug(self, target_size, start, end):
        '''
        All unit are in basepairs
        '''
        offset = random.choice(range(end - start - target_size))
        return start + offset , start + offset + target_size

    def shift_fix(self, target_size, start):
        offset = 0
        return start + offset , start + offset + target_size
    
    def get_unique_elements(dictionary):
        unique_elements = set()
        for key in dictionary:
            unique_elements.update(set(dictionary[key]))
        return list(unique_elements)

    def check_length(self):
        assert len(self.seq.seq) == self.genomic_features[0].length(self.chr_name), f'Sequence {len(self.seq)} and First feature {self.genomic_features[0].length(self.chr_name)} have different length.' 
        assert abs(len(self.seq) / self.res - len(self.mat)) < 2, f'Sequence {len(self.seq) / self.res} and Hi-C {len(self.mat)} have different length.'

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




