import torch
import torch.nn as nn

from . import blocks


# class ConvModel(nn.Module):
#     def __init__(self, num_genomic_features, mid_hidden=256):
#         super(ConvModel, self).__init__()
#         print('Initializing ConvModel')
#         self.encoder = blocks.EncoderSplit(num_genomic_features, output_size=mid_hidden, num_blocks=12)
#         self.decoder = blocks.Decoder(mid_hidden * 2)
#         self.adpool = nn.AdaptiveAvgPool2d((209, 209))

#     def forward(self, x):
#         '''
#         Input feature:
#         batch_size, length * res, feature_dim
#         '''
#         x = self.move_feature_forward(x).float()
#         x = self.encoder(x)
#         x = self.diagonalize(x)
#         x = self.decoder(x).squeeze(1)
#         x = self.adpool(x)
#         return x

#     def move_feature_forward(self, x):
#         '''
#         input dim:
#         bs, img_len, feat
#         to: 
#         bs, feat, img_len
#         '''
#         return x.transpose(1, 2).contiguous()

#     def diagonalize(self, x):
#         # print("x::", x.shape)
#         x_i = x.unsqueeze(2).repeat(1, 1, 256, 1)
#         # print("x_i::", x_i.shape)
#         x_j = x.unsqueeze(3).repeat(1, 1, 1, 256)
#         # print("x_j::", x_j.shape)
#         input_map = torch.cat([x_i, x_j], dim=1)
#         # print("input_map::", input_map.shape)
#         return input_map


class ConvTransModel(nn.Module):
    def __init__(self, encoder_split=blocks.EncoderSplit(), num_classes=23, mid_hidden=256, record_attn=False):
        super(ConvTransModel, self).__init__()
        print('Initializing Discriminator ConvTransModel')
        self.encoder_split = encoder_split
        self.attn = blocks.AttnModule(hidden=mid_hidden, record_attn=record_attn)
        self.decoder = blocks.Decoder(mid_hidden*2)
        self.record_attn = record_attn
        self.adpool = nn.AdaptiveAvgPool2d((209, 209))
        
        # Discriminator branch
        self.pos_embed = nn.Linear(3, 12)
        self.resnet = blocks.ResNet(blocks.ResBlock, [2, 2, 2, 2], num_classes=num_classes)
        self.fc = nn.Linear(12 + num_classes, num_classes)

        # InsulationScoreModule
        self.ismodule = blocks.InsulationScoreModule()

    def forward(self, seq, atac, real_depth, ctcf, pos_start, pos_end, chrom_num):
        seq = self.move_feature_forward(seq).float()
        fused_feat = self.encoder_split(seq=seq, atac=atac, real_depth=real_depth, ctcf=ctcf)
        # print(fused_feat.shape)
        fused_feat = self.move_feature_forward(fused_feat)
        if self.record_attn:
            fused_feat, attn_weights = self.attn(fused_feat)
        else:
            fused_feat = self.attn(fused_feat)
        fused_feat = self.move_feature_forward(fused_feat)
        fused_feat = self.diagonalize(fused_feat)
        # print(fused_feat.shape)
        x = self.decoder(fused_feat).squeeze(1)
        x = self.adpool(x)
        
        # Ignore diagonal
        diag_mask = torch.ones_like(x)
        diag_mask = torch.tril(diag_mask, 0) * torch.triu(diag_mask, -0)
        diag_mask = torch.where(diag_mask == 0, 1, 0) 
        x = diag_mask * x
        x = torch.clamp(x,0,9.99,out=None) + torch.eye(x.shape[1], dtype=x.dtype, device=x.device) * (1e-6)

        # InsulationScoreModule
        insulation = self.ismodule(x)

        # Discriminator
        mat = x.unsqueeze(1)
        start_pos = pos_start.unsqueeze(-1)
        end_pos = pos_end.unsqueeze(-1)
        chrom_num = chrom_num.unsqueeze(-1)
        pos_info = torch.cat([start_pos, end_pos, chrom_num], dim=1)
        pos_info = pos_info.float()
        # print(pos_info.shape)
        pos_feat = self.pos_embed(pos_info)
        mat_feat = self.resnet(mat)
        com_feat = torch.cat([pos_feat, mat_feat], dim=1)
        pred_cls = self.fc(com_feat)
        if self.record_attn:
            return x, pred_cls, insulation, attn_weights
        else:
            return x, pred_cls, insulation

    def move_feature_forward(self, x):
        return x.transpose(1, 2).contiguous()

    def diagonalize(self, x):
        x_i = x.unsqueeze(2).repeat(1, 1, 256, 1)
        x_j = x.unsqueeze(3).repeat(1, 1, 1, 256)
        input_map = torch.cat([x_i, x_j], dim=1)
        return input_map