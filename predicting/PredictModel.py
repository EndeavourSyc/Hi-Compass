import torch
import torch.nn as nn

from . import blocks



class PredictModel(nn.Module):
    def __init__(self, encoder_split=blocks.EncoderSplit(), num_classes=30, mid_hidden=256, record_attn=False):
        super(PredictModel, self).__init__()
        print('Initializing PredictModel')
        self.encoder_split = encoder_split
        self.attn = blocks.AttnModule(hidden=mid_hidden, record_attn=record_attn)
        self.decoder = blocks.Decoder(mid_hidden*2)
        self.record_attn = record_attn
        self.adpool = nn.AdaptiveAvgPool2d((209, 209))
        

    def forward(self, seq, atac, real_depth, ctcf,):
        seq = self.move_feature_forward(seq).float()
        fused_feat = self.encoder_split(seq=seq, atac=atac, real_depth=real_depth, ctcf=ctcf)
        fused_feat = self.move_feature_forward(fused_feat)
        if self.record_attn:
            fused_feat, attn_weights = self.attn(fused_feat)
        else:
            fused_feat = self.attn(fused_feat)
        fused_feat = self.move_feature_forward(fused_feat)
        fused_feat = self.diagonalize(fused_feat)
        x = self.decoder(fused_feat).squeeze(1)
        x = self.adpool(x)
        
        # Ignore diagonal
        diag_mask = torch.ones_like(x)
        diag_mask = torch.tril(diag_mask, 0) * torch.triu(diag_mask, -0)
        diag_mask = torch.where(diag_mask == 0, 1, 0) 
        x = diag_mask * x
        x = torch.clamp(x,0,10,out=None) + torch.eye(x.shape[1], dtype=x.dtype, device=x.device) * (1e-6)

        if self.record_attn:
            return x, attn_weights
        else:
            return x

    def move_feature_forward(self, x):
        return x.transpose(1, 2).contiguous()

    def diagonalize(self, x):
        x_i = x.unsqueeze(2).repeat(1, 1, 256, 1)
        x_j = x.unsqueeze(3).repeat(1, 1, 1, 256)
        input_map = torch.cat([x_i, x_j], dim=1)
        return input_map