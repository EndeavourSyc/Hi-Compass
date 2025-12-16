#!/usr/bin/env python3
import argparse
import sys
from hicompass.commands import preprocess_atac, preprocess_hic_norm, preprocess_hic_to_npz, training, predicting

def main():
    """Hi-Compass main command"""
    parser = argparse.ArgumentParser(
        description='Hi-Compass: Cell-type-specific chromatin interaction prediction',
        usage='hicompass <command> [<args>]'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    
    # Preprocess
    # Preprocess ATAC command
    preprocess_atac_parser = subparsers.add_parser(
        'preprocess-atac',
        help='Preprocess ATAC-seq data (BAM to multi-depth BigWig)'
    )
    preprocess_atac.configure_parser(preprocess_atac_parser)
    
    # Preprocess Hi-C normalization command
    preprocess_hic_norm_parser = subparsers.add_parser(
        'preprocess-hic-norm',
        help='Normalize Hi-C data (contrast stretching, step 1 of 2)'
    )
    preprocess_hic_norm.configure_parser(preprocess_hic_norm_parser)
    
    # Preprocess Hi-C to NPZ command
    preprocess_hic_to_npz_parser = subparsers.add_parser(
        'preprocess-hic-to-npz',
        help='Convert Hi-C to NPZ format (diagonal extraction, step 2 of 2)'
    )
    preprocess_hic_to_npz.configure_parser(preprocess_hic_to_npz_parser)
    
    # Preprocess DNA sequence command
    
    
    #Train
    training_parser = subparsers.add_parser(
        'training',
        help='Train Hi-Compass with given processed dataset'
    )
    training.configure_parser(training_parser)
    
    #Predict
    predicting_parser = subparsers.add_parser(
        'predicting',
        help='Predict Hi-C contact matrices from ATAC-seq and genomic features with a cooler style output'
    )
    predicting.configure_parser(predicting_parser)
    # Parse arguments
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        sys.exit(1)
    
    # Execute command
    if args.command == 'preprocess-atac':
        preprocess_atac.run(args)
    elif args.command == 'preprocess-hic-norm':
        preprocess_hic_norm.run(args)
    elif args.command == 'preprocess-hic-to-npz':
        preprocess_hic_to_npz.run(args)
    elif args.command == 'training':
        training.run(args)
    elif args.command == 'predicting':
        predicting.run(args)

if __name__ == "__main__":
    main()