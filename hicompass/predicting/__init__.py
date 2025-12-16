"""
Hi-Compass Predicting Module

This module provides functionality for Hi-C prediction from genomic features.
"""

from .PredictDataset import PredictDataset
from .PredictModel import PredictModel


__all__ = [
    'PredictDataset',
    'PredictModel',
]
