#!/usr/bin/env python3
"""
Hi-C to NPZ conversion command for Hi-Compass.

This is step 2 of Hi-C preprocessing, converting normalized cool files
to NPZ format for training. Should be run after preprocess-hic-norm.
"""

import logging
from ..preprocess import HiCToNPZConverter

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def configure_parser(parser):
    """Configure argument parser for preprocess-hic-to-npz command."""
    
    # Required arguments
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input cool or mcool file (preferably from preprocess-hic-norm)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory for NPZ files'
    )
    
    # Optional arguments (with strong defaults)
    parser.add_argument(
        '--resolution', '-r',
        type=int,
        default=10000,
        help='Resolution in bp (default: 10000, only 10kb recommended)'
    )
    parser.add_argument(
        '--window', '-w',
        type=int,
        default=256,
        help='Number of diagonals to extract (default: 256, minimum recommended)'
    )
    parser.add_argument(
        '--no-balance',
        action='store_true',
        help='Do not use balanced matrix (no effect if input is from norm step)'
    )
    parser.add_argument(
        '--chr-list',
        nargs='+',
        help='List of chromosomes to process (default: all chromosomes)'
    )


def run(args):
    """Run Hi-C to NPZ conversion."""
    
    try:
        converter = HiCToNPZConverter(
            input_cool=args.input,
            output_dir=args.output,
            resolution=args.resolution,
            window_size=args.window,
            balance=not args.no_balance,
            chr_list=args.chr_list
        )
        
        output_files = converter.process()
        
        print(f"✓ Successfully converted to NPZ format")
        print(f"  Output: {args.output}")
        print(f"  Files: {len(output_files)}")
        
    except Exception as e:
        logger.error(f"Conversion failed: {e}")
        raise