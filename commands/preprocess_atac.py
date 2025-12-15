#!/usr/bin/env python3
"""
ATAC-seq preprocessing command for Hi-Compass.
"""

import logging
from pathlib import Path
import json

from ..preprocess import ATACPreprocessor

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def configure_parser(parser):
    """Configure argument parser for preprocess-atac command."""
    
    parser.add_argument(
        '--input', '-i',
        required=True,
        help='Input BAM/SAM file (bulk ATAC-seq)'
    )
    parser.add_argument(
        '--cell-type', '-c',
        required=True,
        help='Cell type name, only used for output file name (e.g., GM12878, K562 or whatever you like)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output directory'
    )
    parser.add_argument(
        '--chrom-sizes', '-s',
        required=True,
        help='Chromosome sizes file (e.g., hg38.chrom.sizes)'
    )
    
    # Depth specification
    depth_group = parser.add_mutually_exclusive_group()
    depth_group.add_argument(
        '--depths', '-d',
        help='Custom depths (comma-separated, @file.txt, or @file.json)'
    )
    depth_group.add_argument(
        '--min-depth',
        type=float,
        default=2e5,
        help='Minimum depth for range mode (default: 2e5)'
    )
    
    parser.add_argument(
        '--max-depth',
        type=float,
        default=2e7,
        help='Maximum depth for range mode (default: 2e7)'
    )
    parser.add_argument(
        '--step',
        type=float,
        default=2e4,
        help='Step size for range mode (default: 2e4)'
    )
    parser.add_argument(
        '--no-bulk',
        action='store_true',
        help='Skip bulk BigWig generation'
    )
    parser.add_argument(
        '--keep-intermediate',
        action='store_true',
        help='Keep intermediate files'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed (default: 42)'
    )


def parse_depths(depth_str: str):
    """Parse depths from string, file, or JSON."""
    if depth_str.startswith('@'):
        file_path = Path(depth_str[1:])
        if not file_path.exists():
            raise FileNotFoundError(f"Depth file not found: {file_path}")
        
        if file_path.suffix == '.json':
            with open(file_path) as f:
                return json.load(f)
        else:
            with open(file_path) as f:
                return [float(line.strip()) for line in f if line.strip()]
    
    depth_str = depth_str.strip('[]')
    return [float(d.strip()) for d in depth_str.split(',')]


def run(args):
    """Run ATAC-seq preprocessing."""
    
    # Parse depth list if provided
    depth_list = None
    depth_mode = 'range'
    
    if args.depths:
        depth_list = parse_depths(args.depths)
        depth_mode = 'list'
        logger.info(f"Using {len(depth_list)} custom depths")
    
    # Initialize preprocessor
    try:
        preprocessor = ATACPreprocessor(
            input_bam=args.input,
            cell_type=args.cell_type,
            output_dir=args.output,
            chrom_sizes=args.chrom_sizes,
            depth_mode=depth_mode,
            depth_list=depth_list,
            min_depth=args.min_depth,
            max_depth=args.max_depth,
            step=args.step,
            seed=args.seed,
            keep_intermediate=args.keep_intermediate
        )
        
        # Process
        bigwig_files = preprocessor.process(include_bulk=not args.no_bulk)
        
        print(f"\n✓ Successfully generated {len(bigwig_files)} BigWig files")
        print(f"Output directory: {args.output}")
        
    except Exception as e:
        logger.error(f"Preprocessing failed: {e}")
        raise