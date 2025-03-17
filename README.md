# Hi-Compass

Hi-Compass is a tool for predicting chromatin interactions (Hi-C matrices) from ATAC-seq data.

## Installation

### Prerequisites

**Important Note**: Hi-Compass requires PyTorch (version >= 1.13.1), but does not install it automatically since PyTorch installation methods vary depending on your system, CUDA version, and hardware. Please install PyTorch manually following the instructions at [pytorch.org](https://pytorch.org/get-started/locally/) before installing Hi-Compass.

### Installing via pip (recommended)

```bash
pip install hicompass
```

## Dependencies

Hi-Compass depends on the following Python packages:

- torch (>= 1.13.1) - **must be installed manually**
- numpy
- pandas
- pyBigWig
- scikit-image
- cooler

These dependencies (except for PyTorch) will be automatically installed when using pip.

## Usage

### Command Line Interface

After installation, you can use the `hicompass` command:

```bash
hicompass predict --cell_type CELL_TYPE --atac_path PATH_TO_ATAC_BW_FILE --ctcf_path PATH_TO_CTCF_BW_FILE --dna_dir_path PATH_TO_DNA_DIR --model_path PATH_TO_MODEL
```

#### predict subcommand parameters

Required parameters:
- `--cell_type`: Name of the cell type for prediction
- `--atac_path`: Path to the ATAC-seq .bw file
- `--ctcf_path`: Path to the general CTCF ChIP-seq .bw file
- `--dna_dir_path`: Path to the DNA sequence fa.gz directory
- `--model_path`: Path to the model weights file

Optional parameters:
- `--output_path`: Output directory path (default: ./output/<cell_type>)
- `--device`: Device to use for computation (e.g., "cuda:0" or "cpu", default: "cpu")
- `--chromosomes`: Chromosomes to process (e.g., "1,3,5" or "1-22", default: "1-22")
- `--depth`: Sequencing depth parameter (default: 800000)
- `--stride`: Stride value for prediction (default: 50)
- `--batch_size`: Batch size for prediction (default: 2)
- `--num_workers`: Number of worker threads for data loading (default: 16)
- `--resolution`: Resolution of the output Hi-C matrix in bp (default: 10000)
- `--window_size`: Window size for prediction in bp (default: 2097152)

### Python API

You can also use Hi-Compass in Python scripts:

```python
from hicompass.models import ConvTransModel
from hicompass.chromosome_dataset import ChromosomeDataset
from hicompass.utils import merge_hic_segment

# Create dataset
dataset = ChromosomeDataset(
    chr_name_list=['chr1', 'chr2'],
    atac_path_list=['/path/to/atac.bw'],
    dna_dir_path='/path/to/dna/sequences',
    general_ctcf_bw_path='/path/to/ctcf.bw',
    stride=50,
    depth=800000
)

# Load model
model = ConvTransModel()
model.load_state_dict(torch.load('/path/to/model.pth'))

# Further processing...
```

## Examples

Process a specific cell type with required parameters:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth
```

Process only chromosomes 1, 3, and 5:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth \
                  --chromosomes 1,3,5
```

Specify a custom output path and use GPU:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth \
                  --output_path /my/custom/output/directory \
                  --device cuda:0
```

## Output

The tool generates a .cool file containing the predicted Hi-C matrix for the specified cell type.
This file can be visualized using tools like Juicebox, HiCExplorer or just matplotlib.pyplot.imshow.


## License

Hi-Compass is released under the MIT License.