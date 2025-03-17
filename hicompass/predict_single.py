import torch
from torch.utils.data import DataLoader
from hicompass.chromosome_dataset_predict import ChromosomeDataset
import os
import pandas as pd
from hicompass.HiCompass_models import ConvTransModel
import numpy as np
from skimage.transform import resize
from functools import reduce
import time
# time.sleep(1800 * 8)
from hicompass.utils import chrome_size_dict, make_cooler_bins_chr, get_chr_stack, merge_hic_segment, cal_depth, cal_size
torch.multiprocessing.set_sharing_strategy('file_system')
device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu") # 推理就只用单卡就行

chr_list = ['chr'+ str(i) for i in range(1, 23)]
stride = 50 # keep 50 diag
cell_type= 'Tcell'
cell_type = 'CD8_Naive~2'
cell_type_name = cell_type.split('~')[0]
cell_type_list = [
    # 'transverse_colon', 'sun449', 'osteogenesis', 'k562', 'imr90', 'huh7','huh1', 'hmec', 'hffc6', 'hep3b', 
                  'h1_hesc','gm12878', 'gastrocnemius', 'foresking_fibroblasts', 'dnd41', 'ct27_stem', 'ct27_evt', 'a549'
                  ]
for cell_type_name in cell_type_list:
    if not os.path.isdir(f'/cluster/home/Yuanchen/project/scHiC/dataset/ATAC/{cell_type_name}'):
        continue
    model_weight_path = '/cluster/home/Yuanchen/project/scHiC/weight/TJ/cross/model_2.pth'
    if not os.path.exists(f'/cluster/home/Yuanchen/project/scHiC/paper_figure/fig2/cell_line_predict_result_80/{cell_type_name}'):
        os.system(f'mkdir /cluster/home/Yuanchen/project/scHiC/paper_figure/fig2/cell_line_predict_result_80/{cell_type_name}/')
    save_path = f'/cluster/home/Yuanchen/project/scHiC/paper_figure/fig2/cell_line_predict_result_80/{cell_type_name}/{cell_type_name}.cool'
    bw_path = f"/cluster/home/Yuanchen/project/scHiC/dataset/ATAC/{cell_type_name}/{cell_type_name}_ATAC_80e4/{cell_type_name}_ATAC_80e4.bw"
    # bw_path = f'/cluster/home/Yuanchen/project/scHiC/dataset/scATAC/output/atac_v1_pbmc_10k/outs/hic_cluster_n=10_new/{cell_type}/{cell_type}.bw'
    depth = 80e4


    model = ConvTransModel()
    model.to(device)
    checkpoint = torch.load(model_weight_path, map_location=device)
    model_weights = checkpoint
    model.load_state_dict(model_weights, strict=True)
    model.eval()

    cache_barcode_path_list = [bw_path]
    startTime = time.time()
    dataset = ChromosomeDataset(chr_name_list=chr_list, 
                                    atac_path_list=cache_barcode_path_list,
                                    stride=stride, 
                                    depth=depth,use_aug=False)
    print("load dataset in %f s" % (time.time() - startTime))
    dataloader = DataLoader(dataset, batch_size=2, shuffle=False, num_workers=16)
    output_dict = {}
    for chr_name in chr_list:
        output_dict[chr_name] = []
    startTime = time.time()
    for step, data in enumerate(dataloader):
        seq, atac, real_depth, ctcf, start, start_ratio, end, end_ratio, chr_name, chr_name_ratio = data
        sub_time = time.time()
        mat_pred, pred_cls, insulation = model(seq.to(device), atac.to(device), real_depth.to(device), ctcf.to(device), start_ratio.to(device), end_ratio.to(device), chr_name_ratio.to(device)) 
        for i in range(seq.shape[0]):
            result = mat_pred[i].cpu().detach().numpy()
            result = np.clip(result, a_max=10, a_min=0) * 10
            result = np.triu(result)
            np.fill_diagonal(result, 0)
            output_dict[str(chr_name[i])].append([start[i].cpu(), end[i].cpu(), result])
    print("Run model in %f s" % (time.time() - startTime))
    startTime = time.time()
    merge_hic_segment(output_dict,save_path=save_path, window_size=2097152, resolution=10000)
    print("Run merge in %f s" % (time.time() - startTime))