# Hi-Compass

Hi-Compass is a tool for predicting chromatin interactions (Hi-C matrices) from ATAC-seq data.

## Installation

### Prerequisites

**Important Note**: Hi-Compass requires PyTorch (version >= 1.13.1), but does not install it automatically since PyTorch installation methods vary depending on your system, CUDA version, and hardware. Please install PyTorch manually following the instructions at [pytorch.org](https://pytorch.org/get-started/locally/) before using Hi-Compass.

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

### Source data 
https://zenodo.org/records/15037064
Including basic model weight and other constant input data. A demo GM12878 ATAC-seq bw file with depth=840000 is also included.

## Usage

### Command Line Interface

After installation, you can use the `hicompass` command:

```bash
hicompass predict --cell_type CELL_TYPE --atac_path PATH_TO_ATAC_BW_FILE --ctcf_path PATH_TO_CTCF_BW_FILE --dna_dir_path PATH_TO_DNA_DIR --omit_regions_path OMIT_REGIONS_PATH --model_path PATH_TO_MODEL
```

#### predict subcommand parameters

Required parameters:
- `--cell_type`: Name of the cell type for prediction
- `--atac_path`: Path to the ATAC-seq .bw file
- `--ctcf_path`: Path to the general CTCF ChIP-seq .bw file
- `--dna_dir_path`: Path to the DNA sequence fa.gz directory
- `--model_path`: Path to the model weights file
- `--omit_regions_path`: Path to the chromosome omit regions file. This file should be a three-column CSV file with tab ('\t') as the separator. The three columns represent chromosome, start, and end positions, respectively.

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
from hicompass.utils import merge_hic_segment, pseudo_weight

# Data source
atac_path = '/path/to/atac.bw'
chr_name_list = ['chr1', 'chr2']
dna_dir_path = '/path/to/dna_sequences_dir/'
ctcf_path = '/path/to/ctcf.bw'
save_path = '/path/to/save_cooler.cool'

# Load model
model = ConvTransModel()
model.load_state_dict(torch.load('/path/to/model.pth'))

    # Create dataset
    print(f"Creating dataset with ATAC-seq data from {args.atac_path}")
    start_time = time.time()
    dataset = ChromosomeDataset(
        chr_name_list=chr_name_list,
        atac_path_list=[atac_path],
        dna_dir_path=dna_dir_path,
        general_ctcf_bw_path=ctcf_path,
        stride=50,
        depth=800000
    )

    # Create dataloader
    dataloader = DataLoader(
        dataset,
        batch_size=16,
        shuffle=False,
        num_workers=4
    )

    # Initialize output dictionary
    output_dict = {chr_name: [] for chr_name in chr_list}

    # Run prediction
    with torch.no_grad():
        for step, data in enumerate(dataloader):
            seq, atac, real_depth, ctcf, start, start_ratio, end, end_ratio, chr_name, chr_name_ratio = data
            mat_pred, pred_cls, _ = model(
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

    # Merge segments and save
    merge_hic_segment(
        output_dict,
        save_path=save_path,
        window_size=2097152,
        resolution=10000
    )

    # Add pseudo weight
    pseudo_weight(save_path, save_path, weight=1)
```

## Examples

Process a specific cell type with required parameters:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --omit_regions_path /path/to/omit_region_file \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth
```

Process only chromosomes 1, 3, and 5:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --omit_regions_path /path/to/omit_region_file \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth \
                  --chromosomes 1,3,5
                    
```

Specify a custom output path and use GPU:

```bash
hicompass predict --cell_type gm12878 \
                  --atac_path /path/to/gm12878_ATAC.bw \
                  --omit_regions_path /path/to/omit_region_file \
                  --ctcf_path /path/to/ctcf.bw \
                  --dna_dir_path /path/to/dna/sequences \
                  --model_path /path/to/model.pth \
                  --output_path /my/custom/output/directory \
                  --device cuda:0
```

## Output

The tool generates a .cool file containing the predicted Hi-C matrix for the specified cell type.
This file can be visualized using tools like Juicebox, HiCExplorer or just matplotlib.pyplot.imshow.

## Advanced tutorial

An advanced Hi-Compass usage example demonstrating a de novo scATAC-seq-based prediction and further analysis workflow: from performing meta cell Hi-C prediction to clustering results, loop calling, and loop annotation. The complete tutorial is available at http://wulab.bjmu.edu.cn/hicompass. 

## License

Hi-Compass is released under the MIT License.
