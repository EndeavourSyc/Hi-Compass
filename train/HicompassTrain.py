import sys
import torch
import argparse
import numpy as np
import pytorch_lightning as pl
import pytorch_lightning.callbacks as callbacks
from .HicompassModel import ConvTransModel
from .blocks import EncoderSplit 
from .HicompassDataset import ChromosomeDataset
import os
import torch.nn as nn
import matplotlib.pyplot as plt
from piq import InformationWeightedSSIMLoss
import re
# Human genome (hg38) chromosome sizes

def main():
    torch.set_float32_matmul_precision('high') 
    torch.backends.cudnn.enabled = True
    args = init_parser()
    init_training(args)

def init_parser():
    parser = argparse.ArgumentParser(description='Hi-Compass Training Module.')

    parser.add_argument('--seed', dest='run_seed', default=114514, type=int)
    parser.add_argument('--save-path', dest='run_save_path', default='checkpoints')
    parser.add_argument('--data-root', dest='data_root', required=True)
    parser.add_argument('--genome', dest='genome', default='hg38')

    parser.add_argument('--cell-type', dest='cell_type', nargs='+', required=True, type=str)
    parser.add_argument('--cell-type-valid', dest='cell_type_valid', nargs='+', required=True, type=str)
    parser.add_argument('--train-depth', dest='train_depth', nargs='+', required=True, type=str)
    parser.add_argument('--valid-depth', dest='valid_depth', nargs='+', required=True, type=str)
    parser.add_argument('--train-chrom', dest='train_chrom', nargs='+', required=True, type=str)
    parser.add_argument('--valid-chrom', dest='valid_chrom', nargs='+', required=True, type=str)
    
    parser.add_argument('--gpu-id', dest='gpu_id', nargs='+', default=None, type=int)
    parser.add_argument('--num-gpu', dest='num_gpu', default=None, type=int)
    
    parser.add_argument('--max-epochs', dest='trainer_max_epochs', default=100, type=int)
    parser.add_argument('--save-top-n', dest='trainer_save_top_n', default=20, type=int)
    parser.add_argument('--save-step-period', dest='save_step_period', default=100, type=int)

    parser.add_argument('--batch-size', dest='dataloader_batch_size', default=2, type=int)
    parser.add_argument('--num-workers', dest='dataloader_num_workers', default=16, type=int)
    parser.add_argument('--ddp-disabled', dest='dataloader_ddp_disabled', action='store_false')
    parser.add_argument('--ckpt-path', dest='ckpt_path', default=None)
    
    parser.add_argument('--use-augmentation', dest='use_augmentation', action='store_true')
    parser.add_argument('--centromere-filter', dest='centromere_filter', action='store_true')
    
    args = parser.parse_args()
    args = _validate_gpu_config(args)
    
    return args

def _validate_gpu_config(args):
    if args.gpu_id is None and args.num_gpu is None:
        args.gpu_id = [0]
        args.trainer_num_gpu = 1
    elif args.gpu_id is not None and args.num_gpu is not None:
        if len(args.gpu_id) != args.num_gpu:
            raise ValueError(f"GPU config mismatch: {len(args.gpu_id)} vs {args.num_gpu}")
        args.trainer_num_gpu = len(args.gpu_id)
    elif args.gpu_id is not None:
        args.trainer_num_gpu = len(args.gpu_id)
    elif args.num_gpu is not None:
        args.gpu_id = list(range(args.num_gpu))
        args.trainer_num_gpu = args.num_gpu
    
    
    return args

def init_training(args):
    os.makedirs(f'{args.run_save_path}/models', exist_ok=True)
    os.makedirs(f'{args.run_save_path}/csv', exist_ok=True)
    os.makedirs(f'{args.run_save_path}/figure_result', exist_ok=True)
    os.makedirs(f'{args.run_save_path}/model_self', exist_ok=True)
    
    checkpoint_callback = callbacks.ModelCheckpoint(
        dirpath=f'{args.run_save_path}/models',
        save_top_k=args.trainer_save_top_n,
        monitor='val_loss'
    )

    lr_monitor = callbacks.LearningRateMonitor(logging_interval='epoch')
    csv_logger = pl.loggers.CSVLogger(save_dir=f'{args.run_save_path}/csv')

    pl_module = TrainModule(args)
    

    if args.trainer_num_gpu > 1:
        from pytorch_lightning.strategies import DDPStrategy
        strategy = DDPStrategy(find_unused_parameters=True,)
        devices = args.gpu_id
    else:
        strategy = 'auto'
        devices = args.gpu_id[0]
    
    pl_trainer = pl.Trainer(
        strategy=strategy,
        accelerator="gpu", 
        devices=devices,
        gradient_clip_val=1,
        gradient_clip_algorithm='norm', 
        logger=csv_logger,
        callbacks=[checkpoint_callback, lr_monitor],
        max_epochs=args.trainer_max_epochs,
    )
    
    trainloader = pl_module.get_dataloader(args, mode='train')
    valloader = pl_module.get_dataloader(args, mode='val')
    
    pl_trainer.fit(model=pl_module, train_dataloaders=trainloader, val_dataloaders=valloader)



class TrainModule(pl.LightningModule):

    def __init__(self, args):
        super().__init__()
        self.args = args
        self.HG38_CHROM_SIZES = {
        'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555,
        'chr5': 181538259, 'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 
        'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309, 
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 
        'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
        self.MM10_CHROM_SIZES = {
            'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116,
            'chr5': 151834684, 'chr6': 149736546, 'chr7': 145441459, 'chr8': 129401213,
            'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
            'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768,
            'chr17': 94987271, 'chr18': 90702639, 'chr19': 61431566, 'chrX': 171031299,
            'chrY': 91744698}
        self.sorted_cell_lines = sorted(list(set(args.cell_type + args.cell_type_valid)) , key=lambda x: x.lower())
        self.cell_line_to_idx = {cell: idx for idx, cell in enumerate(self.sorted_cell_lines)}
        self.sorted_chroms = sorted(list(set(args.train_chrom + args.valid_chrom)), key=lambda x: int(re.search(r'\d+', x).group()))
        self.chr_name_to_index = {chr_name: idx + 1 for idx, chr_name in enumerate(self.sorted_chroms)}
        self.num_cell_types = len(args.cell_type)
        if self.num_cell_types>1:
            self.use_discriminator = True
        else:
            self.use_discriminator = False
        if args.genome in ['hg38', 'HG38']:
            self.chrom_sizes = self.HG38_CHROM_SIZES
        elif args.genome in ['mm10', 'MM10']:
            self.chrom_sizes = self.MM10_CHROM_SIZES
        else:
            chrom_sizes_file = os.path.join(args.data_root, 'chromsize', f'{args.genome}.chrom.sizes')
            if not os.path.exists(chrom_sizes_file):
                raise FileNotFoundError(f"You are using a customed genome dataset rather than hg38 or mm10, but the chrom sizes file not found: {chrom_sizes_file}")
            self.chrom_sizes = self.load_chrom_sizes_from_file(self, chrom_sizes_file)
        
        self.model = self.get_model(args)
        
        self.save_hyperparameters()
        
        self.use_full_loss = False
        self.mse_threshold = 1.5
        self.epoch_mse_losses = []
        
        self.train_loss_tracker = {ct: {'mse': [], 'sim': [], 'cls': []} 
                                   for ct in self.cell_line_to_idx.keys()}
        self.val_loss_tracker = {ct: {'mse': [], 'sim': []} 
                                 for ct in self.cell_line_to_idx.keys()}
        
        self.csv_file = f"{args.run_save_path}/loss_by_cell_type_train.csv"
        self.csv_file_valid = f"{args.run_save_path}/loss_by_cell_type_valid.csv"
        
        self._init_csv_files()
        
        print(f"\n{'='*70}")
        print(f"Training Configuration:")
        print(f"{'='*70}")
        print(f"Cell line mapping: {self.cell_line_to_idx}")
        print(f"Chromosome mapping size: {len(self.chr_name_to_index)}")
        print(f"Discriminator enabled: {self.use_discriminator}")
        print(f"Number of cell types in training set: {self.num_cell_types}")
        print(f"{'='*70}\n")
        
        
    def load_chrom_sizes_from_file(chrom_sizes_file):
        chrom_sizes = {}
        with open(chrom_sizes_file, 'r') as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2:
                    chrom_sizes[parts[0]] = int(parts[1])
        return chrom_sizes
    
    def _init_csv_files(self):
        with open(self.csv_file, 'w') as f:
            header = ['Epoch'] + [f"{ct}_mse,{ct}_sim,{ct}_cls" for ct in self.train_loss_tracker.keys()]
            f.write(','.join(header) + '\n')
        
        with open(self.csv_file_valid, 'w') as f:
            header = ['Epoch'] + [f"{ct}_mse,{ct}_sim" for ct in self.val_loss_tracker.keys()]
            f.write(','.join(header) + '\n')

    def forward(self, seq, atac, real_depth, ctcf, start_ratio, end_ratio, chr_name_ratio):
        return self.model(seq, atac, real_depth, ctcf, start_ratio, end_ratio, chr_name_ratio)

    def training_step(self, batch, batch_idx):
        seq, atac, real_depth, ctcf, mat, start_ratio, end_ratio, chr_name_ratio, cell_type = batch
        mat = mat.float()
        
        outputs, cls, ins = self(seq, atac, real_depth, ctcf, start_ratio, end_ratio, chr_name_ratio)
        outputs = torch.nan_to_num(outputs, nan=0.0)
        mat = torch.nan_to_num(mat, nan=0.0)
        
        criterion_mse = nn.MSELoss()
        
        if not self.use_full_loss:
            if self.use_discriminator:
                for param in self.model.resnet.parameters():
                    param.requires_grad = False
                for param in self.model.fc.parameters():
                    param.requires_grad = False
            
            loss_mse = criterion_mse(outputs, mat)
            loss = loss_mse
            
            self.epoch_mse_losses.append(loss_mse.item())
            self._track_losses(outputs, mat, cell_type, loss_mse, torch.tensor(0.0), torch.tensor(0.0), mode='train')
            
            metrics = {
                'train_loss': loss, 
                'mse_loss': loss_mse, 
                'sim_loss': torch.tensor(0.0, device=self.device),
                'cls_loss': torch.tensor(0.0, device=self.device),
                'training_stage': torch.tensor(0.0, device=self.device)
            }
        else:
            if self.use_discriminator:
                for param in self.model.resnet.parameters():
                    param.requires_grad = True
                for param in self.model.fc.parameters():
                    param.requires_grad = True
            
            criterion_sim = InformationWeightedSSIMLoss(data_range=10)
            criterion_cls = nn.CrossEntropyLoss()
            
            loss_mse = criterion_mse(outputs, mat)
            
            try:
                loss_sim = criterion_sim(outputs.unsqueeze(1), mat.unsqueeze(1))
                if torch.isnan(loss_sim):
                    loss_sim = loss_mse
            except:
                loss_sim = loss_mse
            
            if self.use_discriminator:
                loss_cls = criterion_cls(cls, cell_type)
                if torch.isnan(loss_cls):
                    loss_cls = torch.tensor(0.0, device=self.device)
            else:
                loss_cls = torch.tensor(0.0, device=self.device)
            
            loss = 0.4 * loss_mse + 0.4 * loss_sim + 0.2 * loss_cls
            
            self.epoch_mse_losses.append(loss_mse.item())
            self._track_losses(outputs, mat, cell_type, loss_mse, loss_sim, loss_cls, mode='train')
            
            metrics = {
                'train_loss': loss, 
                'mse_loss': loss_mse, 
                'sim_loss': loss_sim,
                'cls_loss': loss_cls, 
                'training_stage': torch.tensor(1.0, device=self.device)
            }
        
        self.log_dict(metrics, batch_size=atac.shape[0], prog_bar=True, sync_dist=True)
        
        if batch_idx % self.args.save_step_period == 0:
            if self.trainer.is_global_zero:
                self._save_visualization(outputs, mat, cell_type, real_depth, 
                                        start_ratio, end_ratio, chr_name_ratio, 
                                        batch_idx, mode='train')
        
        return loss

    def validation_step(self, batch, batch_idx):
        seq, atac, real_depth, ctcf, mat, start_ratio, end_ratio, chr_name_ratio, cell_type = batch
        mat = mat.float()
        
        outputs, cls, ins = self(seq, atac, real_depth, ctcf, start_ratio, end_ratio, chr_name_ratio)
        outputs = torch.nan_to_num(outputs, nan=0.0)
        mat = torch.nan_to_num(mat, nan=0.0)
        
        criterion_mse = nn.MSELoss()
        criterion_sim = InformationWeightedSSIMLoss(data_range=10)
        
        loss_mse = criterion_mse(outputs, mat)
        
        try:
            loss_sim = criterion_sim(outputs.unsqueeze(1), mat.unsqueeze(1))
            if torch.isnan(loss_sim):
                loss_sim = loss_mse
        except:
            loss_sim = loss_mse
        
        loss = loss_mse
        
        self._track_losses(outputs, mat, cell_type, loss_mse, loss_sim, torch.tensor(0.0, device=self.device), mode='val')
        
        metrics = {'val_loss': loss, 'val_mse': loss_mse, 'val_sim': loss_sim}
        self.log_dict(metrics, batch_size=atac.shape[0], prog_bar=True, sync_dist=True)
        
        if batch_idx % max(1, self.args.save_step_period // 10) == 0:
            self._save_visualization(outputs, mat, cell_type, real_depth,
                                    start_ratio, end_ratio, chr_name_ratio,
                                    batch_idx, mode='val')
        
        return loss

    def _track_losses(self, outputs, mat, cell_type, loss_mse, loss_sim, loss_cls, mode='train'):
        tracker = self.train_loss_tracker if mode == 'train' else self.val_loss_tracker
        idx_to_cell = {idx: cell for cell, idx in self.cell_line_to_idx.items()}
        
        for ct_idx in torch.unique(cell_type):
            ct_idx_val = ct_idx.item()
            if ct_idx_val not in idx_to_cell:
                continue
            
            cell_name = idx_to_cell[ct_idx_val]
            if cell_name not in tracker:
                continue
            
            mask = cell_type == ct_idx
            if not mask.any():
                continue
            
            outputs_ct = outputs[mask]
            mat_ct = mat[mask]
            
            criterion_mse = nn.MSELoss()
            mse_ct = criterion_mse(outputs_ct, mat_ct)
            tracker[cell_name]['mse'].append(mse_ct.item())
            
            if mode == 'train':
                tracker[cell_name]['sim'].append(loss_sim.item() if torch.is_tensor(loss_sim) else loss_sim)
                tracker[cell_name]['cls'].append(loss_cls.item() if torch.is_tensor(loss_cls) else loss_cls)
            else:
                tracker[cell_name]['sim'].append(loss_sim.item() if torch.is_tensor(loss_sim) else loss_sim)


    def _save_visualization(self, outputs, mat, cell_type, real_depth, 
                        start_ratio, end_ratio, chr_name_ratio, batch_idx, mode='train'):
        if self.trainer.is_global_zero:
            try:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
                
                idx_to_chr = {idx: chr_name for chr_name, idx in self.chr_name_to_index.items()}
                chr_idx = int(chr_name_ratio[0].item() * len(self.chr_name_to_index))
                chr_name = idx_to_chr.get(chr_idx, 'chrUnknown')
                
                chr_size = self.chrom_sizes.get(chr_name, 200000000)
                start = int(start_ratio[0].item() * chr_size)
                end = int(end_ratio[0].item() * chr_size)
                
                idx_to_cell = {idx: cell for cell, idx in self.cell_line_to_idx.items()}
                cell_name = idx_to_cell.get(cell_type[0].item(), 'Unknown')
                
                output_array = outputs[0].detach().cpu()
                output_array = torch.clamp(output_array, 0, 10).numpy()
                
                im1 = ax1.imshow(output_array, cmap='Reds')
                stage_text = "MSE Only" if not self.use_full_loss else "Full Loss"
                ax1.set_title(f'Prediction\n{chr_name}:{start}-{end}\nStage: {stage_text}')
                ax1.axis('off')
                fig.colorbar(im1, ax=ax1)
                
                mat_array = mat[0].detach().cpu().numpy()
                im2 = ax2.imshow(mat_array, cmap='Reds')
                depth_str = f"{int(real_depth[0].item())//10000}e4"
                ax2.set_title(f'Ground Truth\n{cell_name}_{depth_str}')
                ax2.axis('off')
                fig.colorbar(im2, ax=ax2)
                
                fig.suptitle(f'Epoch {self.current_epoch} - {mode.capitalize()} Batch {batch_idx}')
                plt.tight_layout()
                
                save_path = f"{self.args.run_save_path}/figure_result/epoch{self.current_epoch}_{mode}_batch{batch_idx}.png"
                plt.savefig(save_path, dpi=100, bbox_inches='tight')
                plt.close('all')
            except Exception as e:
                print(f"Warning: Failed to save visualization: {e}")
                plt.close('all')

    def on_train_epoch_end(self):
        if not self.use_full_loss and len(self.epoch_mse_losses) > 0:
            avg_mse = sum(self.epoch_mse_losses) / len(self.epoch_mse_losses)
            
            if self.trainer.is_global_zero:
                print(f"\nEpoch {self.current_epoch} - Average MSE: {avg_mse:.4f}")
            
            if avg_mse < self.mse_threshold:
                if self.trainer.is_global_zero:
                    print(f"Switching to full loss training")
                self.use_full_loss = True
                
                for param in self.model.parameters():
                    param.requires_grad = True
        
        self.epoch_mse_losses = []
        
        if self.trainer.is_global_zero:
            self._log_cell_type_losses(mode='train')
            torch.save(self.model.state_dict(), 
                      f'{self.args.run_save_path}/model_self/model_epoch{self.current_epoch}.pth')
        
        self.train_loss_tracker = {ct: {'mse': [], 'sim': [], 'cls': []} 
                                   for ct in self.cell_line_to_idx.keys()}

    def on_validation_epoch_end(self):
        if self.trainer.is_global_zero:
            self._log_cell_type_losses(mode='val')
        
        self.val_loss_tracker = {ct: {'mse': [], 'sim': []} 
                                 for ct in self.cell_line_to_idx.keys()}

    def _log_cell_type_losses(self, mode='train'):
        tracker = self.train_loss_tracker if mode == 'train' else self.val_loss_tracker
        csv_file = self.csv_file if mode == 'train' else self.csv_file_valid
        
        with open(csv_file, 'a') as f:
            row = [f'{mode.capitalize()}_Epoch_{self.current_epoch}']
            
            for ct in tracker.keys():
                mse_mean = np.mean(tracker[ct]['mse']) if tracker[ct]['mse'] else 0
                sim_mean = np.mean(tracker[ct]['sim']) if tracker[ct]['sim'] else 0
                
                if mode == 'train':
                    cls_mean = np.mean(tracker[ct]['cls']) if tracker[ct]['cls'] else 0
                    row.append(f"{mse_mean:.4f},{sim_mean:.4f},{cls_mean:.4f}")
                else:
                    row.append(f"{mse_mean:.4f},{sim_mean:.4f}")
            
            f.write(','.join(row) + '\n')

    def configure_optimizers(self):
        optimizer = torch.optim.Adam(self.parameters(), lr=2e-4, weight_decay=0)
        
        warmup_scheduler = torch.optim.lr_scheduler.LinearLR(
            optimizer, start_factor=0.01, end_factor=1.0, total_iters=1
        )
        cosine_scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=self.args.trainer_max_epochs - 1, eta_min=1e-6
        )
        scheduler = torch.optim.lr_scheduler.SequentialLR(
            optimizer, schedulers=[warmup_scheduler, cosine_scheduler], milestones=[1]
        )
        
        return {'optimizer': optimizer, 'lr_scheduler': {'scheduler': scheduler, 'interval': 'epoch'}}

    def get_dataset(self, args, mode):
        if mode == "train":
            dataset = ChromosomeDataset(
                chr_list=args.train_chrom,
                cell_line_list=args.cell_type,
                depth_list=args.train_depth,
                data_root=args.data_root,
                genome=args.genome,
                use_augmentation=args.use_augmentation,
                centromere_filter=args.centromere_filter,
                verbose=False
            )
        elif mode == 'val':
            dataset = ChromosomeDataset(
                chr_list=args.valid_chrom,
                cell_line_list=args.cell_type_valid,
                depth_list=args.valid_depth,
                data_root=args.data_root,
                genome=args.genome,
                use_augmentation=False,
                centromere_filter=args.centromere_filter,
                verbose=False
            )
        
        return dataset

    def get_dataloader(self, args, mode):
        dataset = self.get_dataset(args, mode)
        shuffle = True
        batch_size = args.dataloader_batch_size
        num_workers = args.dataloader_num_workers

        if not args.dataloader_ddp_disabled:
            gpus = args.trainer_num_gpu
            batch_size = int(args.dataloader_batch_size / gpus)
            num_workers = int(args.dataloader_num_workers / gpus)

        dataloader = torch.utils.data.DataLoader(
            dataset,
            shuffle=shuffle,
            batch_size=batch_size,
            num_workers=num_workers,
            pin_memory=True,
            prefetch_factor=2,
            persistent_workers=True
        )
        return dataloader

    def get_model(self, args):
        encoder_split = EncoderSplit()
        model = ConvTransModel(
            encoder_split=encoder_split,
            num_classes=self.num_cell_types,
            mid_hidden=256,
            record_attn=False
        )
        if args.trainer_num_gpu > 1:
            model = torch.nn.SyncBatchNorm.convert_sync_batchnorm(model)

        
        if args.ckpt_path is not None:
            checkpoint = torch.load(args.ckpt_path, map_location='cpu')
            model_weights = checkpoint.get('state_dict', checkpoint)
            model_weights = {key.replace('model.', ''): value for key, value in model_weights.items()}
            
            # check discriminator fit
            current_state = model.state_dict()
            incompatible_keys = []
            for key in ['fc.weight', 'fc.bias', 'resnet.fc.weight', 'resnet.fc.bias']:
                if key in model_weights and key in current_state:
                    if model_weights[key].shape != current_state[key].shape:
                        incompatible_keys.append(key)
            
            if incompatible_keys:
                print(f"\n{'='*70}")
                print(f"⚠ Discriminator shape mismatch: {incompatible_keys}")
                print(f"Removing all discriminator weights for clean initialization")
                print(f"{'='*70}\n")
                model_weights = {k: v for k, v in model_weights.items() 
                               if not any(k.startswith(p) for p in ['resnet.', 'fc.', 'pos_embed.'])}
            
            model.load_state_dict(model_weights, strict=False)
            print('Loaded checkpoint!')
        
        return model


if __name__ == '__main__':
    main()