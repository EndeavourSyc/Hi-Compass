# predicting.py
#!/usr/bin/env python3
"""Hi-Compass Prediction Command"""

import os
import time
import argparse
import torch
import numpy as np
from torch.utils.data import DataLoader

from ..predicting.PredictModel import PredictModel
from ..predicting.PredictDataset import PredictDataset
from ..predicting.utils import merge_hic_segment, pseudo_weight, get_chrom_sizes


def configure_parser(parser):
    """Configure argument parser for prediction command"""
    
    # Required arguments
    parser.add_argument(
        '--model-path',
        type=str,
        required=True,
        help='Path to trained model checkpoint'
    )
    parser.add_argument(
        '--atac-path',
        type=str,
        required=True,
        nargs='+',
        help='Path(s) to ATAC-seq bigWig file(s)'
    )
    parser.add_argument(
        '--ctcf-path',
        type=str,
        required=True,
        help='Path to CTCF ChIP-seq bigWig file'
    )
    parser.add_argument(
        '--dna-dir',
        type=str,
        required=True,
        help='Directory containing DNA sequences (chr*.fa.gz)'
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help='Output path for predicted Hi-C (.cool file)'
    )
    
    # Genome configuration
    parser.add_argument(
        '--genome',
        type=str,
        default='hg38',
        help='Genome build (hg38, mm10, or custom). Default: hg38'
    )
    parser.add_argument(
        '--chrom-sizes',
        type=str,
        default=None,
        help='Path to .chrom.sizes file (required if genome is not hg38/mm10)'
    )
    parser.add_argument(
        '--centromere-bed',
        type=str,
        default=None,
        help='Path to centromere/telomere BED file (optional)'
    )
    
    # Prediction parameters
    parser.add_argument(
        '--chromosomes',
        type=str,
        default='1-22',
        help='Chromosomes to predict (e.g., "1,2,3" or "1-22"). Default: 1-22'
    )
    parser.add_argument(
        '--depth',
        type=float,
        default=80000,
        help='Sequencing depth for normalization. Default: 80000'
    )
    parser.add_argument(
        '--stride',
        type=int,
        default=50,
        help='Stride for sliding windows (bins). Default: 50'
    )
    
    # Computational parameters
    parser.add_argument(
        '--device',
        type=str,
        default='cpu',
        help='Device (cpu, cuda, cuda:0). Default: cpu'
    )
    parser.add_argument(
        '--batch-size',
        type=int,
        default=2,
        help='Batch size. Default: 2'
    )
    parser.add_argument(
        '--num-workers',
        type=int,
        default=16,
        help='Number of data loading workers. Default: 16'
    )
    parser.add_argument(
        '--no-pseudo-weight',
        action='store_true',
        help='Do not add pseudo weight to output'
    )


def parse_chromosome_input(chrom_str):
    """Parse chromosome input string to list of chromosome names"""
    chromosomes = []
    
    if '-' in chrom_str:
        start, end = map(int, chrom_str.split('-'))
        chromosomes = [f'chr{i}' for i in range(start, end + 1)]
    else:
        chrom_list = chrom_str.split(',')
        chromosomes = [f'chr{c.strip()}' for c in chrom_list]
    
    return chromosomes


def run(args):
    """Run Hi-C prediction"""
    start_time = time.time()
    
    # Set device
    device = torch.device(
        args.device if torch.cuda.is_available() and 'cuda' in args.device else 'cpu'
    )
    
    print(f"\nHi-Compass Prediction")
    print(f"Device: {device}")
    
    # Parse chromosomes
    chr_list = parse_chromosome_input(args.chromosomes)
    print(f"Chromosomes: {', '.join(chr_list)}")
    print(f"Genome: {args.genome}")
    
    # Load chromosome sizes
    chrom_sizes = get_chrom_sizes(args.genome, args.chrom_sizes)
    
    # Verify input files
    for atac_path in args.atac_path:
        if not os.path.exists(atac_path):
            raise FileNotFoundError(f"ATAC file not found: {atac_path}")
    
    if not os.path.exists(args.ctcf_path):
        raise FileNotFoundError(f"CTCF file not found: {args.ctcf_path}")
    
    if not os.path.exists(args.dna_dir):
        raise FileNotFoundError(f"DNA directory not found: {args.dna_dir}")
    
    if args.centromere_bed and not os.path.exists(args.centromere_bed):
        raise FileNotFoundError(f"Centromere BED not found: {args.centromere_bed}")
    
    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"Output: {args.output}\n")
    
    # Load model
    print(f"Loading model...")
    model = PredictModel()
    model.to(device)
    
    checkpoint = torch.load(args.model_path, map_location=device)
    
    if isinstance(checkpoint, dict) and 'state_dict' in checkpoint:
        state_dict = checkpoint['state_dict']
    else:
        state_dict = checkpoint
    
    # Remove 'model.' prefix
    state_dict = {k.replace('model.', ''): v for k, v in state_dict.items()}
    
    # Remove discriminator weights: pos_embed, resnet, fc
    discriminator_prefixes = ['pos_embed.', 'resnet.', 'fc.']
    filtered_state_dict = {
        k: v for k, v in state_dict.items()
        if not any(k.startswith(prefix) for prefix in discriminator_prefixes)
    }
    
    model.load_state_dict(filtered_state_dict, strict=True)
    model.eval()
    print(f"Model loaded ({time.time() - start_time:.1f}s)\n")
    
    # Create dataset
    print(f"Loading dataset...")
    dataset_start = time.time()
    
    dataset = PredictDataset(
        chr_name_list=chr_list,
        atac_path_list=args.atac_path,
        dna_dir_path=args.dna_dir,
        general_ctcf_bw_path=args.ctcf_path,
        omit_regions_file_path=args.centromere_bed,
        genome=args.genome,
        chrom_sizes_file=args.chrom_sizes,
        stride=args.stride,
        depth=args.depth
    )
    
    print(f"Dataset loaded: {len(dataset)} intervals ({time.time() - dataset_start:.1f}s)\n")
    
    # Create dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=False,
        num_workers=args.num_workers
    )
    
    # Initialize output
    output_dict = {chr_name: [] for chr_name in chr_list}
    
    # Run prediction
    print(f"Predicting...")
    pred_start = time.time()
    
    with torch.no_grad():
        for step, data in enumerate(dataloader):
            seq, atac, real_depth, ctcf, start, end, chr_name = data
            
            mat_pred = model(
                seq.to(device),
                atac.to(device),
                real_depth.to(device),
                ctcf.to(device)
            )
            
            for i in range(seq.shape[0]):
                result = mat_pred[i].cpu().detach().numpy()
                result = np.clip(result, a_max=10, a_min=0) * 10
                result = np.triu(result)
                np.fill_diagonal(result, 0)
                output_dict[str(chr_name[i])].append([
                    start[i].cpu(), 
                    end[i].cpu(), 
                    result
                ])
            
            if (step + 1) % 50 == 0:
                elapsed = time.time() - pred_start
                progress = (step + 1) / len(dataloader) * 100
                print(f"  {step+1}/{len(dataloader)} ({progress:.1f}%) - {elapsed:.0f}s")
    
    print(f"Prediction complete ({time.time() - pred_start:.1f}s)\n")
    
    # Merge and save
    print(f"Saving to {args.output}...")
    merge_start = time.time()
    
    merge_hic_segment(
        output_dict,
        save_path=args.output,
        chrom_sizes=chrom_sizes,
        window_size=2097152,
        resolution=10000
    )
    
    print(f"Saved ({time.time() - merge_start:.1f}s)")
    
    # Add pseudo weight
    if not args.no_pseudo_weight:
        pseudo_weight(args.output, args.output, weight=1)
        print(f"Pseudo weight added")
    
    total_time = time.time() - start_time
    print(f"\nTotal time: {total_time:.1f}s ({total_time/60:.1f}min)")
    print(f"Output: {args.output}\n")