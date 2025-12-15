#!/usr/bin/env python3
"""
Hi-C normalization preprocessing command for Hi-Compass.
"""

import logging
from ..preprocess import HiCNormalizer

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def configure_parser(parser):
    """Configure argument parser for preprocess-hic-norm command."""
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input cool file'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output cool file (contrast-stretched)'
    )
    parser.add_argument(
        '--genome', '-g',
        default='hg38',
        choices=['hg38', 'mm10', 'custom'],
        help='Genome assembly (default: hg38). Use "custom" for other species.'
    )
    parser.add_argument(
        '--chrom-sizes', '-s',
        help='Chromosome sizes file (TSV: chrom<tab>size) or cool file to extract from. '
             'Required when --genome=custom, optional otherwise.'
    )
    parser.add_argument(
        '--resolution', '-r',
        type=int,
        default=10000,
        help='Resolution in bp (default: 10000)'
    )
    parser.add_argument(
        '--force-resolution',
        action='store_true',
        help='Skip resolution warning for non-10kb data'
    )
    parser.add_argument(
        '--percentile-min',
        type=float,
        default=2.0,
        help='Lower percentile for contrast stretching (default: 2.0)'
    )
    parser.add_argument(
        '--percentile-max',
        type=float,
        default=98.0,
        help='Upper percentile for contrast stretching (default: 98.0, recommended)'
    )
    parser.add_argument(
        '--sample-start',
        type=int,
        default=4000,
        help='Start bin for percentile sampling (default: 4000, recommended)'
    )
    parser.add_argument(
        '--sample-size',
        type=int,
        default=256,
        help='Sample region size in bins (default: 256, recommended)'
    )
    parser.add_argument(
        '--no-balance',
        action='store_true',
        help='Do not use balanced matrix'
    )
    parser.add_argument(
        '--chr-list',
        nargs='+',
        help='List of chromosomes to process (default: all chromosomes)'
    )


def run(args):
    """Run Hi-C normalization preprocessing."""
    
    # Validate chrom-sizes requirement
    if args.genome == 'custom' and not args.chrom_sizes:
        logger.error(
            "Error: --chrom-sizes is required when --genome=custom\n"
            "Provide either:\n"
            "  1. A chrom.sizes file (TSV format: chrom<tab>size)\n"
            "  2. A cool file to extract chromosome sizes from"
        )
        raise ValueError("Missing required argument: --chrom-sizes")
    
    try:
        normalizer = HiCNormalizer(
            input_cool=args.input,
            output_cool=args.output,
            resolution=args.resolution,
            genome=args.genome,
            chrom_sizes=args.chrom_sizes,
            percentile_min=args.percentile_min,
            percentile_max=args.percentile_max,
            sample_region_start=args.sample_start,
            sample_region_size=args.sample_size,
            balance=not args.no_balance,
            chr_list=args.chr_list,
            force_resolution=args.force_resolution
        )
        
        output_path = normalizer.process()
        
        print(f"\n✓ Successfully normalized Hi-C data")
        print(f"Output: {output_path}")
        
    except Exception as e:
        logger.error(f"Normalization failed: {e}")
        raise