#!/usr/bin/env python3
import argparse
import torch
from torch.utils.data import DataLoader
import os
import numpy as np
import time
import cooler
from hicompass.chromosome_dataset_predict import ChromosomeDataset
from hicompass.HiCompass_models import ConvTransModel
from hicompass.utils import chrome_size_dict, merge_hic_segment, pseudo_weight


def configure_parser(parser):
    """
    Configure the argument parser for the predict command
    """
    # Required arguments
    parser.add_argument('--cell_type', type=str, required=True,
                        help='Cell type name for the prediction')
    parser.add_argument('--atac_path', type=str, required=True,
                        help='Path to the ATAC-seq .bw file')
    parser.add_argument('--ctcf_path', type=str, required=True,
                        help='Path to the general CTCF ChIP-seq .bw file')
    parser.add_argument('--dna_dir_path', type=str, required=True,
                        help='Path to the DNA sequence fa.gz dir file')
    parser.add_argument('--model_path', type=str, required=True,
                        help='Path to the model weights')

    # Optional arguments with defaults
    parser.add_argument('--output_path', type=str,
                        help='Output directory path (default: ./output/<cell_type>)')
    parser.add_argument('--device', type=str, default='cpu',
                        help='decide gpu or cpu usage (e.g., "cuda:0" or "cpu"), default cpu')
    parser.add_argument('--chromosomes', type=str, default='1-22',
                        help='Chromosomes to process (e.g., "1,3,5" or "1-22")')
    parser.add_argument('--depth', type=float, default=80e4,
                        help='Sequencing depth parameter')
    parser.add_argument('--stride', type=int, default=50,
                        help='Stride value for prediction (default: 50)')
    parser.add_argument('--batch_size', type=int, default=2,
                        help='Batch size for prediction (default: 2)')
    parser.add_argument('--num_workers', type=int, default=16,
                        help='Number of worker threads for data loading (default: 16)')
    parser.add_argument('--resolution', type=int, default=10000,
                        help='Resolution of the output Hi-C matrix in bp (default: 10000)')
    parser.add_argument('--window_size', type=int, default=2097152,
                        help='Window size for prediction in bp (default: 2097152)')

    return parser


def parse_chromosome_input(chrom_str):
    """Parse chromosome input to get list of chromosomes."""
    chromosomes = []

    if '-' in chrom_str:
        # Range format: "1-22"
        start, end = map(int, chrom_str.split('-'))
        chromosomes = [f'chr{i}' for i in range(start, end + 1)]
    else:
        # List format: "1,2,3,X,Y"
        chrom_list = chrom_str.split(',')
        chromosomes = [f'chr{c}' for c in chrom_list]

    return chromosomes


def run(args):
    """Run prediction with the specified arguments."""
    # Configure torch multiprocessing
    torch.multiprocessing.set_sharing_strategy('file_system')

    # Set device
    device = torch.device(args.device if torch.cuda.is_available() and "cuda" in args.device else "cpu")
    print(f"Using device: {device}")

    # Parse chromosomes
    chr_list = parse_chromosome_input(args.chromosomes)
    print(f"Processing chromosomes: {', '.join(chr_list)}")

    # Set output path
    if args.output_path is None:
        output_dir = f"./output/{args.cell_type}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        save_path = f"{output_dir}/{args.cell_type}.cool"
    else:
        output_dir = args.output_path
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        save_path = f"{output_dir}/{args.cell_type}.cool"

    print(f"Output will be saved to: {save_path}")

    # Load model
    print(f"Loading model from {args.model_path}")
    model = ConvTransModel()
    model.to(device)
    checkpoint = torch.load(args.model_path, map_location=device)
    model.load_state_dict(checkpoint, strict=True)
    model.eval()

    # Create dataset
    print(f"Creating dataset with ATAC-seq data from {args.atac_path}")
    start_time = time.time()
    dataset = ChromosomeDataset(
        chr_name_list=chr_list,
        atac_path_list=[args.atac_path],
        dna_dir_path=args.dna_dir_path,
        general_ctcf_bw_path=args.ctcf_path,
        stride=args.stride,
        depth=args.depth,
        use_aug=False
    )
    print(f"Dataset loaded in {time.time() - start_time:.2f} seconds")

    # Create dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers
    )

    # Initialize output dictionary
    output_dict = {chr_name: [] for chr_name in chr_list}

    # Run prediction
    start_time = time.time()
    with torch.no_grad():
        for step, data in enumerate(dataloader):
            if step % 10 == 0:
                print(f"Processing batch {step}/{len(dataloader)}")

            seq, atac, real_depth, ctcf, start, start_ratio, end, end_ratio, chr_name, chr_name_ratio = data

            # Forward pass
            mat_pred, pred_cls, insulation = model(
                seq.to(device),
                atac.to(device),
                real_depth.to(device),
                ctcf.to(device),
                start_ratio.to(device),
                end_ratio.to(device),
                chr_name_ratio.to(device)
            )

            # Process results
            for i in range(seq.shape[0]):
                result = mat_pred[i].cpu().detach().numpy()
                result = np.clip(result, a_max=10, a_min=0) * 10
                result = np.triu(result)
                np.fill_diagonal(result, 0)
                output_dict[str(chr_name[i])].append([start[i].cpu(), end[i].cpu(), result])

    print(f"Prediction completed in {time.time() - start_time:.2f} seconds")

    # Merge segments and save
    start_time = time.time()
    print(f"Merging segments and saving to {save_path}")
    merge_hic_segment(
        output_dict,
        save_path=save_path,
        window_size=args.window_size,
        resolution=args.resolution
    )
    print(f"Merging and saving completed in {time.time() - start_time:.2f} seconds")

    # Add pseudo weight
    pseudo_weight(save_path, save_path, weight=1)
    print(f'Pseudo weight=1 added')
    print(f"All operations completed. Results saved to {save_path}")