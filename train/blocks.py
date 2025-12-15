import torch
import torch.nn as nn
import numpy as np
import copy


class ConvBlock(nn.Module):
    def __init__(self, size, stride=2, hidden_in=64, hidden=64):
        super(ConvBlock, self).__init__()
        pad_len = int(size / 2)
        self.scale = nn.Sequential(
            nn.Conv1d(hidden_in, hidden, size, stride, pad_len),
            nn.BatchNorm1d(hidden),
            nn.ReLU(),
        )
        self.res = nn.Sequential(
            nn.Conv1d(hidden, hidden, size, padding=pad_len),
            nn.BatchNorm1d(hidden),
            nn.ReLU(),
            nn.Conv1d(hidden, hidden, size, padding=pad_len),
            nn.BatchNorm1d(hidden),
        )
        self.relu = nn.ReLU()

    def forward(self, x):
        scaled = self.scale(x)
        identity = scaled
        res_out = self.res(scaled)
        out = self.relu(res_out + identity)
        return out
    


class ATACDepthEncoder(nn.Module):
    def __init__(self, in_channel=1, output_size=128, filter_size=5, num_blocks=12, depth_ranges=None):
        super(ATACDepthEncoder, self).__init__()
        self.filter_size = filter_size
        self.conv_start = nn.Sequential(
            nn.Conv1d(in_channel, 32, 3, 2, 1, bias=False).float(),
            nn.BatchNorm1d(32),
            nn.ReLU(),
        )
        self.conv_start_1d = nn.Conv1d(in_channel, 32, 3, 2, 1)
        self.bn1d = nn.BatchNorm1d(32)
        self.activate = nn.ReLU()
        sizes_3 = [3] * num_blocks
        hiddens = [32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256, 256]
        hidden_ins = [32, 32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256]
        self.res_blocks = self.get_res_blocks(num_blocks, sizes_3, hidden_ins, hiddens)

        self.conv_end = nn.Conv1d(256, output_size, 1)
        
        # 添加深度嵌入层
        if depth_ranges is None:
            depth_ranges = [500000, 800000, 1200000, 3000000, 100000000, 200000000]
        self.depth_ranges = depth_ranges
        self.depth_embed = nn.Embedding(len(depth_ranges) + 1, 1)

    def forward(self, atac, depth):
        # 将深度值映射到区间索引
        depth_idx = torch.bucketize(depth, torch.tensor(self.depth_ranges, device=depth.device))
        
        # 将区间索引输入到嵌入层
        depth_feat = self.depth_embed(depth_idx)
        # depth_feat = depth_feat.unsqueeze(2).expand(-1, -1, atac.size(1))
        
        # 将ATAC序列与深度特征连接
        # print('atac shape: ', atac.shape)
        # print('depth_feat shape: ', depth_feat.shape)

        # x =(atac * depth_feat).unsqueeze(1).to(torch.float32) # (atac * depth_feat).unsqueeze(1).double()
        x = torch.mul(atac, depth_feat).unsqueeze(1).to(torch.float32) 
        # x =atac.unsqueeze(1).to(torch.float32)
        # print(type(x),x.dtype)
        # print('start', x.shape)
        x = self.conv_start_1d(x) 
        # print('conv_start_1d', x.shape)
        x = self.bn1d(x)
        # print('bn1d', x.shape)
        x = self.activate(x)
        x = self.res_blocks(x)
        # print('res_blocks', x.shape)
        out = self.conv_end(x)
        # print('conv_end', x.shape)
        return out

    def get_res_blocks(self, n, sizes, his, hs):
        blocks = []
        for i, filter_size, h, hi in zip(range(n), sizes, hs, his):
            blocks.append(ConvBlock(filter_size, hidden_in=hi, hidden=h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks


class DNAEncoder(nn.Module):
    def __init__(self, in_channel=5, output_size=64, filter_size=5, num_blocks=12):
        super(DNAEncoder, self).__init__()
        self.filter_size = filter_size
        self.conv_start = nn.Sequential(
            nn.Conv1d(in_channel, 32, 3, 2, 1),
            nn.BatchNorm1d(32),
            nn.ReLU(),
        )
        sizes_3 = [3] * num_blocks
        
        hiddens = [32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256, 256]
        hidden_ins = [32, 32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256]
        self.res_blocks = self.get_res_blocks(num_blocks, sizes_3, hidden_ins, hiddens)
        self.conv_end = nn.Conv1d(256, output_size, 1)
        
        # 添加特征权重
        self.feature_weights = nn.Parameter(torch.ones(1, in_channel, 1))

    def forward(self, x):
        # 应用特征权重
        # print(x.shape)
        # x = x * self.feature_weights
        x = torch.mul(x, self.feature_weights)
        x = self.conv_start(x)
        
        x = self.res_blocks(x)
        out = self.conv_end(x)
        return out
    def get_res_blocks(self, n, sizes, his, hs):
        blocks = []
        for i, filter_size, h, hi in zip(range(n), sizes, hs, his):
            blocks.append(ConvBlock(filter_size, hidden_in=hi, hidden=h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks



class CTCFEncoder(nn.Module):
    def __init__(self, in_channel=1, output_size=64, filter_size=5, num_blocks=12):
        super(CTCFEncoder, self).__init__()
        self.filter_size = filter_size
        self.conv_start = nn.Sequential(
            nn.Conv1d(in_channel, 32, 3, 2, 1).float(),
            nn.BatchNorm1d(32),
            nn.ReLU(),
        )
        self.conv_start_1d = nn.Conv1d(in_channel, 32, 3, 2, 1)
        self.bn1d = nn.BatchNorm1d(32)
        self.activate = nn.ReLU()
        sizes_3 = [3] * num_blocks

        hiddens = [32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256, 256]
        hidden_ins = [32, 32, 32, 32, 32, 64, 64, 128, 128, 128, 128, 256]
        self.res_blocks = self.get_res_blocks(num_blocks, sizes_3, hidden_ins, hiddens)

        self.conv_end = nn.Conv1d(256, output_size, 1)

    def forward(self, x):
        x = x.unsqueeze(1).to(torch.float32)
        x = self.conv_start_1d(x) 
        x = self.bn1d(x)
        x = self.activate(x)
        x = self.res_blocks(x)
        out = self.conv_end(x)
        return out

    def get_res_blocks(self, n, sizes, his, hs):
        blocks = []
        for i, filter_size, h, hi in zip(range(n), sizes, hs, his):
            blocks.append(ConvBlock(filter_size, hidden_in=hi, hidden=h))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks

# seq, atac, real_depth, ctcf, mat
class EncoderSplit(nn.Module):
    def __init__(self, atac_depth_encoder=ATACDepthEncoder(), dna_encoder=DNAEncoder(), ctcf_encoder=CTCFEncoder(), fusion_method='attention'):
        super(EncoderSplit, self).__init__()
        self.atac_depth_encoder = atac_depth_encoder
        self.dna_encoder = dna_encoder
        self.ctcf_encoder = ctcf_encoder
        self.fusion_method = fusion_method

        if fusion_method == 'attention':
            self.attention = AttnModule(hidden=256, )
        elif fusion_method == 'gating':
            self.gate = nn.Sequential(
                nn.Linear(256, 256),
                nn.Sigmoid()
            )

    def forward(self, seq, atac, real_depth, ctcf):
        atac_depth_out = self.atac_depth_encoder(atac, real_depth)
        dna_out = self.dna_encoder(seq)
        ctcf_out = self.ctcf_encoder(ctcf)

        if self.fusion_method == 'attention':
            # print(atac_depth_out.shape, dna_out.shape, ctcf_out.shape)
            concat_out = torch.cat([atac_depth_out, dna_out, ctcf_out], dim=1)
            # print(concat_out.shape)
            concat_out = concat_out.permute(0, 2, 1)
            fused_out = self.attention(concat_out)
            fused_out = fused_out.permute(0, 2, 1)
            # print(fused_out.shape)
            # fused_out = fused_out.mean(dim=2)
        elif self.fusion_method == 'gating':
            concat_out = torch.cat([atac_depth_out, dna_out, ctcf_out], dim=1)
            gate_weights = self.gate(concat_out)
            fused_out = gate_weights * concat_out
        else:
            fused_out = torch.cat([atac_depth_out, dna_out, ctcf_out], dim=1)

        return fused_out


class ResBlockDilated(nn.Module):
    def __init__(self, size, hidden=64, stride=1, dil=2):
        super(ResBlockDilated, self).__init__()
        pad_len = dil
        self.res = nn.Sequential(
            nn.Conv2d(hidden, hidden, size, padding=pad_len,
                      dilation=dil),
            nn.BatchNorm2d(hidden),
            nn.ReLU(),
            nn.Conv2d(hidden, hidden, size, padding=pad_len,
                      dilation=dil),
            nn.BatchNorm2d(hidden),
        )
        self.relu = nn.ReLU()

    def forward(self, x):
        identity = x
        res_out = self.res(x)
        out = self.relu(res_out + identity)
        return out


class Decoder(nn.Module):
    def __init__(self, in_channel, hidden=256, filter_size=3, num_blocks=5):
        super(Decoder, self).__init__()
        self.filter_size = filter_size

        self.conv_start = nn.Sequential(
            nn.Conv2d(in_channel, hidden, 3, 1, 1),
            nn.BatchNorm2d(hidden),
            nn.ReLU(),
        )
        self.res_blocks = self.get_res_blocks(num_blocks, hidden)
        self.conv_end = nn.Conv2d(hidden, 1, 1)

    def forward(self, x):
        x = self.conv_start(x)
        x = self.res_blocks(x)
        out = self.conv_end(x)
        return out

    def get_res_blocks(self, n, hidden):
        blocks = []
        for i in range(n):
            dilation = 2 ** (i + 1)
            blocks.append(ResBlockDilated(self.filter_size, hidden=hidden, dil=dilation))
        res_blocks = nn.Sequential(*blocks)
        return res_blocks


class TransformerLayer(torch.nn.TransformerEncoderLayer):
    # Pre-LN structure

    def forward(self, src, src_mask=None, src_key_padding_mask=None):
        # MHA section
        src_norm = self.norm1(src)
        src_side, attn_weights = self.self_attn(src_norm, src_norm, src_norm,
                                                attn_mask=src_mask,
                                                key_padding_mask=src_key_padding_mask)
        src = src + self.dropout1(src_side)

        # MLP section
        src_norm = self.norm2(src)
        src_side = self.linear2(self.dropout(self.activation(self.linear1(src_norm))))
        src = src + self.dropout2(src_side)
        return src, attn_weights


class TransformerEncoder(torch.nn.TransformerEncoder):

    def __init__(self, encoder_layer, num_layers, norm=None, record_attn=False):
        super(TransformerEncoder, self).__init__(encoder_layer, num_layers)
        self.layers = self._get_clones(encoder_layer, num_layers)
        self.num_layers = num_layers
        self.norm = norm
        self.record_attn = record_attn

    def forward(self, src, mask=None, src_key_padding_mask=None):
        r"""Pass the input through the encoder layers in turn.

        Args:
            src: the sequence to the encoder (required).
            mask: the mask for the src sequence (optional).
            src_key_padding_mask: the mask for the src keys per batch (optional).

        Shape:
            see the docs in Transformer class.
        """
        output = src

        attn_weight_list = []

        for mod in self.layers:
            output, attn_weights = mod(output, src_mask=mask, src_key_padding_mask=src_key_padding_mask)
            attn_weight_list.append(attn_weights.unsqueeze(0).detach())
        if self.norm is not None:
            output = self.norm(output)

        if self.record_attn:
            return output, torch.cat(attn_weight_list)
        else:
            return output

    def _get_clones(self, module, N):
        return torch.nn.modules.ModuleList([copy.deepcopy(module) for i in range(N)])


class PositionalEncoding(nn.Module):

    def __init__(self, hidden, dropout=0.1, max_len=256):
        super().__init__()
        self.dropout = nn.Dropout(p=dropout)
        position = torch.arange(max_len).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, hidden, 2) * (-np.log(10000.0) / hidden))
        pe = torch.zeros(max_len, 1, hidden)
        pe[:, 0, 0::2] = torch.sin(position * div_term)
        pe[:, 0, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        """
        Args:
            x: Tensor, shape [seq_len, batch_size, embedding_dim]
        """
        x = x + self.pe[:x.size(0)]
        return self.dropout(x)


class AttnModule(nn.Module):
    def __init__(self, hidden=256, layers=8, record_attn=False, inpu_dim=256):
        super(AttnModule, self).__init__()

        self.record_attn = record_attn
        self.pos_encoder = PositionalEncoding(hidden, dropout=0.1)
        encoder_layers = TransformerLayer(hidden,
                                          nhead=8,
                                          dropout=0.1,
                                          dim_feedforward=512,
                                          batch_first=True)
        self.module = TransformerEncoder(encoder_layers,
                                         layers,
                                         record_attn=record_attn)

    def forward(self, x):
        x = self.pos_encoder(x)
        output = self.module(x)
        return output

    def inference(self, x):
        return self.module(x)


class ResBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1):
        super(ResBlock, self).__init__()
        self.conv1 = nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.relu = nn.ReLU(inplace=True)
        self.conv2 = nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn2 = nn.BatchNorm2d(out_channels)

        self.shortcut = nn.Sequential()
        if stride != 1 or in_channels != out_channels:
            self.shortcut = nn.Sequential(
                nn.Conv2d(in_channels, out_channels, kernel_size=1, stride=stride, bias=False),
                nn.BatchNorm2d(out_channels)
            )

    def forward(self, x):
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.bn2(self.conv2(out))
        out += self.shortcut(x)
        out = self.relu(out)
        return out


class ResNet(nn.Module):
    def __init__(self, block, num_blocks, num_classes=2):
        super(ResNet, self).__init__()
        self.in_channels = 64

        self.conv1 = nn.Conv2d(1, 64, kernel_size=3, stride=1, padding=1, bias=False)
        self.bn1 = nn.BatchNorm2d(64)
        self.relu = nn.ReLU(inplace=True)
        self.layer1 = self._make_layer(block, 64, num_blocks[0], stride=1)
        self.layer2 = self._make_layer(block, 128, num_blocks[1], stride=2)
        self.layer3 = self._make_layer(block, 256, num_blocks[2], stride=2)
        self.layer4 = self._make_layer(block, 512, num_blocks[3], stride=2)
        self.avg_pool = nn.AdaptiveAvgPool2d((1, 1))
        self.fc = nn.Linear(512, num_classes)

    def _make_layer(self, block, out_channels, num_blocks, stride):
        strides = [stride] + [1] * (num_blocks - 1)
        layers = []
        for stride in strides:
            layers.append(block(self.in_channels, out_channels, stride))
            self.in_channels = out_channels
        return nn.Sequential(*layers)

    def forward(self, x):
        out = self.relu(self.bn1(self.conv1(x)))
        out = self.layer1(out)
        out = self.layer2(out)
        out = self.layer3(out)
        out = self.layer4(out)
        out = self.avg_pool(out)
        out = out.view(out.size(0), -1)
        out = self.fc(out)
        return out

class InsulationScoreModule(nn.Module):
    def __init__(self, window_size=10):
        super().__init__()
        self.window_size = window_size
        self.avg_pool = nn.AvgPool2d((self.window_size, self.window_size), stride=(1, 1),padding=(0, 0))
        # self.ad_pool = nn.AdaptiveAvgPool1d(209)
        

    def forward(self, x):
        x = self.avg_pool(x)
        # x = x.squeeze(dim=-1)
        x = torch.diagonal(x, offset=self.window_size + 1, dim1=-2, dim2=-1)
        # x = self.ad_pool(x)
        return x
# class InsulationScoreModule(nn.Module):
#     def __init__(self, window_size=10, ignor_diag=2, hidden_dim=256):
#         super(InsulationScoreModule, self).__init__()
#         self.conv1 = nn.Conv2d(1, hidden_dim, kernel_size=3, stride=1, padding=1)
#         self.window_size = window_size
#         self.ignor_diag = ignor_diag - 1
#         self.relu = nn.ReLU(inplace=True)
#         self.conv2 = nn.Conv2d(hidden_dim, 1, kernel_size=(209, 1), stride=1)
#         # self.avg_pool = nn.AdaptiveAvgPool1d(209)

#     def forward(self, x):
#         mask = torch.ones_like(x)
#         mask = torch.tril(mask, self.window_size) * torch.triu(mask, -self.window_size)  # window outlier mask
#         mask_reverse = torch.tril(mask, self.ignor_diag) * torch.triu(mask, -self.ignor_diag)
#         mask_reverse = torch.where(mask_reverse == 0, 1, 0)  # ignore diag mask
#         x = x * mask * mask_reverse
#         # print(x)
#         x = x.unsqueeze(1)  # (batch_size, 1, 256, 256)
#         x = self.conv1(x)
#         x = self.relu(x)  # (batch_size, hidden_dim, 256, 256)
#         x = self.conv2(x)  # (batch_size, 1, 1, 256)
#         x = x.squeeze(1).squeeze(1)
#         # x = self.avg_pool(x)
#         return x