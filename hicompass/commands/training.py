#!/usr/bin/env python3
"""
Training command for Hi-Compass.
"""

import logging
import os
import torch
from argparse import Namespace

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def configure_parser(parser):
    """
    Configure argument parser for training command.
    
    Args:
        parser: ArgumentParser or subparser to configure
    """
    # Required arguments
    parser.add_argument(
        '--data-root',
        required=True,
        help='Root directory containing ATAC/HiC/DNA/CTCF subdirectories'
    )
    parser.add_argument(
        '--cell-type',
        nargs='+',
        required=True,
        help='Training cell types (e.g., gm12878 k562)'
    )
    parser.add_argument(
        '--train-chr',
        nargs='+',
        required=True,
        help='Training chromosomes (supports: chr1-chr5, 1-5, train_default, all)'
    )
    parser.add_argument(
        '--valid-chr',
        nargs='+',
        required=True,
        help='Validation chromosomes (supports: chr19-chr20, 19-20, valid_default)'
    )
    
    # Optional arguments with defaults
    parser.add_argument(
        '--genome',
        default='hg38',
        help='Genome build (hg38, mm10, or custom). For custom genome, provide {genome}.chrom.sizes in your_data_root/chromsize'
    )
    parser.add_argument(
        '--cell-type-valid',
        nargs='+',
        default=None,
        help='Validation cell types (default: same as training)'
    )
    parser.add_argument(
        '--train-depth',
        nargs='+',
        default=['bulk'],
        help='Training ATAC-seq depths'
    )
    parser.add_argument(
        '--valid-depth',
        nargs='+',
        default=['bulk'],
        help='Validation ATAC-seq depths'
    )
    
    # Training hyperparameters
    parser.add_argument(
        '--batch-size',
        type=int,
        default=2,
        help='Batch size per GPU (default: 2)'
    )
    parser.add_argument(
        '--num-workers',
        type=int,
        default=4,
        help='Number of data loading workers (default: 4)'
    )
    parser.add_argument(
        '--max-epochs',
        type=int,
        default=100,
        help='Maximum training epochs (default: 100)'
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=114514,
        help='Random seed (default: 114514)'
    )
    
    # GPU configuration
    parser.add_argument(
        '--gpu-id',
        nargs='+',
        type=int,
        default=None,
        help='GPU IDs to use (default: auto-detect)'
    )
    
    # Paths
    parser.add_argument(
        '--save-path',
        default='checkpoints',
        help='Directory to save checkpoints (default: checkpoints)'
    )
    parser.add_argument(
        '--ckpt-path',
        default=None,
        help='Path to checkpoint for resuming training (default: None)'
    )
    
    # Data processing options
    parser.add_argument(
        '--use-augmentation',
        action='store_true',
        help='Enable data augmentation'
    )
    parser.add_argument(
        '--centromere-filter',
        action='store_true',
        help='Filter centromeric regions'
    )
    
    return parser


def run(args):
    """
    Run Hi-Compass Training.
    
    Args:
        args: Parsed arguments from CLI (argparse.Namespace)
    """
    logger.info("="*70)
    logger.info("Hi-Compass Training")
    logger.info("="*70)
    
    try:
        # Set PyTorch optimizations
        torch.set_float32_matmul_precision('high') 
        torch.backends.cudnn.enabled = True
        
        # Validate arguments
        validate_args(args)
        
        # Convert CLI args to training args format
        train_args = convert_args(args)

        
        # Import and run training
        from ..train import init_training
        init_training(train_args)
        
        logger.info("✓ Training completed successfully")
        
    except Exception as e:
        logger.error(f"✗ Training failed: {e}", exc_info=True)
        raise


def validate_args(args):
    """Validate training arguments."""
    logger.info("Validating arguments...")
    
    # Check data root exists
    if not os.path.exists(args.data_root):
        raise FileNotFoundError(f"Data root not found: {args.data_root}")
    
    # Set default cell_type_valid if not provided
    if args.cell_type_valid is None:
        args.cell_type_valid = args.cell_type
        logger.info(f"Using training cell types for validation: {args.cell_type}")
    
    # Check required directories
    required_dirs = {
        'ATAC': os.path.join(args.data_root, 'ATAC', args.genome),
        'HiC': os.path.join(args.data_root, 'HiC', args.genome),
        'DNA': os.path.join(args.data_root, 'DNA', args.genome),
        'CTCF': os.path.join(args.data_root, 'CTCF', args.genome),
    }
    
    missing_dirs = []
    for name, path in required_dirs.items():
        if not os.path.exists(path):
            missing_dirs.append(f"{name}: {path}")
            logger.warning(f"⚠ {name} directory not found: {path}")
    
    if missing_dirs:
        logger.warning(f"Missing {len(missing_dirs)} required directories")
    
    # GPU configuration
    if args.gpu_id is None:
        if torch.cuda.is_available():
            args.gpu_id = [0]
            logger.info("✓ Auto-detected GPU 0")
        else:
            raise RuntimeError("No GPU available and no --gpu-id specified")
    else:
        logger.info(f"✓ Using GPUs: {args.gpu_id}")
    
    # Check chromosome sizes file for custom genome
    if args.genome not in ['hg38', 'mm10']:
        chrom_sizes_file = os.path.join(args.data_root, 'chromsize', f'{args.genome}.chrom.sizes')
        if not os.path.exists(chrom_sizes_file):
            raise FileNotFoundError(
                f"Custom genome '{args.genome}' requires chromosome sizes file at: {chrom_sizes_file}"
            )
        logger.info(f"✓ Found chromosome sizes for custom genome: {args.genome}")
    
    # Parse chromosome arguments (support simplified format)
    logger.info("Parsing chromosome arguments...")
    args.train_chr = parse_chrom_arg(args.train_chr, args.genome, args.data_root)
    args.valid_chr = parse_chrom_arg(args.valid_chr, args.genome, args.data_root)
    
    logger.info(f"  Training: {len(args.train_chr)} chromosomes")
    logger.info(f"    {args.train_chr}")
    logger.info(f"  Validation: {len(args.valid_chr)} chromosomes")
    logger.info(f"    {args.valid_chr}")
    
    # Check for overlap
    overlap = set(args.train_chr) & set(args.valid_chr)
    if overlap:
        logger.warning(f"⚠ Train and validation chromosomes overlap: {overlap}")
    
    logger.info("✓ Arguments validated")


def convert_args(args):
    """
    Convert CLI arguments to training script format.
    
    CLI format → HicompassTrain.py format
    - train_chr → train_chrom
    - max_epochs → trainer_max_epochs
    etc.
    """
    logger.info("Converting arguments to training format...")
    
    train_args = Namespace()
    
    # Basic configuration
    train_args.run_seed = args.seed
    train_args.run_save_path = args.save_path
    train_args.data_root = args.data_root
    train_args.genome = args.genome
    
    # Cell types and chromosomes
    train_args.cell_type = args.cell_type
    train_args.cell_type_valid = args.cell_type_valid
    train_args.train_depth = args.train_depth
    train_args.valid_depth = args.valid_depth
    train_args.train_chrom = args.train_chr  # Note: train_chr → train_chrom
    train_args.valid_chrom = args.valid_chr  # Note: valid_chr → valid_chrom
    
    # GPU configuration
    train_args.gpu_id = args.gpu_id
    train_args.num_gpu = len(args.gpu_id)
    train_args.trainer_num_gpu = train_args.num_gpu
    
    # Training parameters
    train_args.trainer_max_epochs = args.max_epochs
    train_args.trainer_save_top_n = 20  # Fixed default
    train_args.save_step_period = 100   # Fixed default
    
    # DataLoader configuration
    train_args.dataloader_batch_size = args.batch_size
    train_args.dataloader_num_workers = args.num_workers
    train_args.dataloader_ddp_disabled = True  # Enable DDP
    
    # Checkpoint
    train_args.ckpt_path = args.ckpt_path
    
    # Data processing
    train_args.use_augmentation = args.use_augmentation
    train_args.centromere_filter = args.centromere_filter
    
    logger.info("✓ Arguments converted")
    return train_args


def print_config(args):
    """Print training configuration."""
    print("\n" + "="*70)
    print("Training Configuration")
    print("="*70)
    print(f"Data root:          {args.data_root}")
    print(f"Genome:             {args.genome}")
    print(f"Train cells:        {args.cell_type}")
    print(f"Valid cells:        {args.cell_type_valid}")
    print(f"Train chromosomes:  {args.train_chrom[:3]}{'...' if len(args.train_chrom) > 3 else ''} ({len(args.train_chrom)} total)")
    print(f"Valid chromosomes:  {args.valid_chrom} ({len(args.valid_chrom)} total)")
    print(f"Train depths:       {args.train_depth}")
    print(f"Valid depths:       {args.valid_depth}")
    print(f"GPUs:               {args.gpu_id}")
    print(f"Batch size:         {args.dataloader_batch_size}")
    print(f"Num workers:        {args.dataloader_num_workers}")
    print(f"Max epochs:         {args.trainer_max_epochs}")
    print(f"Augmentation:       {args.use_augmentation}")
    print(f"Centromere filter:  {args.centromere_filter}")
    print(f"Save path:          {args.run_save_path}")
    if args.ckpt_path:
        print(f"Resume from:        {args.ckpt_path}")
    print("="*70 + "\n")


def parse_chrom_arg(chrom_arg, genome='hg38', data_root=None):
    """
    Parse chromosome arguments with simplified input support.
    
    Supports:
    1. Range: chr1-chr5 or 1-5 → ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']
    2. Preset (hg38/mm10 only): train_default, valid_default, all, autosome
    3. Exclude: all ^chr19 ^chr20 → all except chr19, chr20
    4. Mix: chr1-chr10 chrX → combined
    
    For custom genomes:
    - Presets not available, use explicit chromosome lists
    - 'all' reads from {genome}.chrom.sizes file
    - Range and exclude still work
    
    Args:
        chrom_arg: List of chromosome arguments
        genome: Genome build (hg38/mm10/custom)
        data_root: Data root directory (required for custom genome with 'all')
        
    Returns:
        Sorted list of chromosome names
    """
    import re
    
    # Preset groups (only for hg38/mm10)
    PRESETS = {
        'hg38': {
            'all': [f'chr{i}' for i in range(1, 23)] + ['chrX'],
            'autosome': [f'chr{i}' for i in range(1, 23)],
            'train_default': [f'chr{i}' for i in range(1, 19)],
            'valid_default': ['chr19', 'chr20'],
            'test_default': ['chr21', 'chr22'],
        },
        'mm10': {
            'all': [f'chr{i}' for i in range(1, 20)] + ['chrX'],
            'autosome': [f'chr{i}' for i in range(1, 20)],
            'train_default': [f'chr{i}' for i in range(1, 16)],
            'valid_default': ['chr17', 'chr18'],
            'test_default': ['chr19'],
        }
    }
    
    def normalize(chrom):
        """Add chr prefix if needed."""
        chrom = chrom.strip()
        return chrom if chrom.startswith('chr') else f'chr{chrom}'
    
    def expand_range(start, end):
        """Expand range like chr1-chr5 or 1-5."""
        try:
            start_num = int(re.search(r'(\d+)', start).group(1))
            end_num = int(re.search(r'(\d+)', end).group(1))
            return [f'chr{i}' for i in range(start_num, end_num + 1)]
        except:
            return []
    
    def sort_key(chrom):
        """Sorting key for chromosomes."""
        c = chrom.replace('chr', '')
        if c.isdigit():
            return (0, int(c))
        elif c == 'X':
            return (1, 0)
        elif c == 'Y':
            return (1, 1)
        return (2, c)
    
    def get_all_chroms_from_file(genome, data_root):
        """Read all chromosomes from chrom.sizes file for custom genome."""
        if not data_root:
            raise ValueError(f"data_root required for custom genome '{genome}' with 'all' preset")
        
        chrom_sizes_file = os.path.join(data_root, 'chromsize', f'{genome}.chrom.sizes')
        if not os.path.exists(chrom_sizes_file):
            raise FileNotFoundError(f"Chrom sizes file not found: {chrom_sizes_file}")
        
        chroms = []
        with open(chrom_sizes_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    chroms.append(parts[0])
        return chroms
    
    # Check if genome is custom (not hg38/mm10)
    is_custom_genome = genome not in PRESETS
    
    result = set()
    exclude = set()
    
    for arg in chrom_arg:
        # Handle exclusion (^chr19)
        if arg.startswith('^'):
            chrom = arg[1:]
            # Check if it's a preset
            if not is_custom_genome and chrom in PRESETS[genome]:
                exclude.update(PRESETS[genome][chrom])
            else:
                exclude.add(normalize(chrom))
            continue
        
        # Handle 'all' keyword
        if arg == 'all':
            if is_custom_genome:
                # For custom genome, read from chrom.sizes file
                all_chroms = get_all_chroms_from_file(genome, data_root)
                result.update(all_chroms)
            else:
                # For hg38/mm10, use preset
                result.update(PRESETS[genome]['all'])
            continue
        
        # Preset group (only for hg38/mm10)
        if not is_custom_genome and arg in PRESETS[genome]:
            result.update(PRESETS[genome][arg])
        
        # Range: chr1-chr5 or 1-5
        elif '-' in arg and arg.count('-') == 1:
            start, end = arg.split('-')
            chroms = expand_range(start, end)
            if chroms:
                result.update(chroms)
            else:
                result.add(normalize(arg))
        
        # Single chromosome
        else:
            # If it's a preset name but custom genome, warn user
            if is_custom_genome and arg in ['train_default', 'valid_default', 'test_default', 'autosome']:
                logger.warning(f"⚠ Preset '{arg}' not available for custom genome '{genome}'. Treating as chromosome name.")
            result.add(normalize(arg))
    
    # Apply exclusions
    result -= exclude
    
    # Sort and return
    return sorted(result, key=sort_key)