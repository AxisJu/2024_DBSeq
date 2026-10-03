# Project adaptations of DOLPHIN: https://github.com/mcgilldinglab/DOLPHIN
# MIT License. Copyright (c) 2024 McGill Ding Lab.
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

def compute(inputs):
    """Compute using supplied in-memory input mappings, helper functions and options."""
    import torch
    import torch.nn as nn
    import pyro
    import pyro.distributions as dist
    from pyro.infer import SVI, Trace_ELBO
    from pyro.optim import Adam
    from torch_geometric.nn import GATConv
    import pyro.poutine as poutine
    pyro.distributions.enable_validation(False)
    pyro.set_rng_seed(0)

    class Gra_Encoder(nn.Module):

        def __init__(self, in_node_fea, gat_channel, nhead, gat_dropout, concat, in_fea, hidden_dim, z_dim, p_dropout):
            super().__init__()
            self.hidden_channels = gat_channel
            self.convs = nn.ModuleList()
            for i, h_dim in enumerate(self.hidden_channels):
                if i == 0 or concat == False:
                    self.convs.append(GATConv(in_node_fea, h_dim, heads=nhead, dropout=gat_dropout, concat=concat))
                else:
                    self.convs.append(GATConv(in_node_fea * nhead, h_dim, heads=nhead, dropout=gat_dropout, concat=concat))
                in_node_fea = h_dim
            self.act = nn.ReLU(inplace=True)
            modules = []
            in_dim = in_fea * gat_channel[-1] * (nhead if concat else 1)
            modules.append(nn.Sequential(nn.Linear(in_dim, hidden_dim[0]), nn.BatchNorm1d(hidden_dim[0], momentum=0.01, eps=0.001), nn.LayerNorm(hidden_dim[0], elementwise_affine=False), nn.ReLU(inplace=True), nn.Dropout(p=p_dropout)))
            if len(hidden_dim) > 1:
                for h_dim in hidden_dim[1:]:
                    modules.append(nn.Sequential(nn.Linear(hidden_dim[0], h_dim), nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), nn.LayerNorm(h_dim, elementwise_affine=False), nn.ReLU(inplace=True), nn.Dropout(p=p_dropout)))
            self.encoder = nn.Sequential(*modules)
            self.fc_mu = nn.Linear(hidden_dim[-1], z_dim)
            self.fc_var = nn.Linear(hidden_dim[-1], z_dim)

        def forward(self, data, batch):
            x = data.x
            for conv in self.convs:
                x = self.act(conv(x, data.edge_index, data.edge_attr))
            x_gat_conv = x.reshape(batch, -1)
            fea_out = self.encoder(x_gat_conv)
            z_loc = self.fc_mu(fea_out)
            z_scale = torch.exp(self.fc_var(fea_out)) + 0.0001
            return (z_loc, z_scale)

    class Fea_Decoder(nn.Module):

        def __init__(self, z_dim, hidden_dim, out_dim):
            super().__init__()
            hidden_dim_decoder = hidden_dim[::-1]
            modules = []
            curr_dim = z_dim
            for h_dim in hidden_dim_decoder:
                modules.append(nn.Sequential(nn.Linear(curr_dim, h_dim), nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), nn.LayerNorm(h_dim, elementwise_affine=False), nn.ReLU(inplace=True), nn.Dropout(p=0)))
                curr_dim = h_dim
            self.decoder = nn.Sequential(*modules)
            self.final = nn.Linear(hidden_dim_decoder[-1], out_dim)

        def forward(self, z):
            hidden = self.decoder(z)
            return torch.sigmoid(self.final(hidden))

    class Adj_Decoder(nn.Module):

        def __init__(self, z_dim, hidden_dim, out_dim):
            super().__init__()
            hidden_dim_decoder = hidden_dim[::-1]
            modules = []
            curr_dim = z_dim
            for h_dim in hidden_dim_decoder:
                modules.append(nn.Sequential(nn.Linear(curr_dim, h_dim), nn.BatchNorm1d(h_dim, momentum=0.01, eps=0.001), nn.LayerNorm(h_dim, elementwise_affine=False), nn.ReLU(inplace=True), nn.Dropout(p=0)))
                curr_dim = h_dim
            self.decoder = nn.Sequential(*modules)
            self.final = nn.Linear(hidden_dim_decoder[-1], out_dim)

        def forward(self, z):
            hidden = self.decoder(z)
            return torch.sigmoid(self.final(hidden))

    class VAE(nn.Module):

        def __init__(self, in_node_fea, gat_channel, nhead, gat_dropout, concat, in_fea, list_gra_enc_hid, z_dim, gra_p_dropout, list_fea_dec_hid, in_adj, list_adj_dec_hid, kl_beta, fea_lambda, adj_lambda):
            super().__init__()
            self.z_dim = z_dim
            self.in_fea = in_fea
            self.in_adj = in_adj
            self.kl_beta = kl_beta
            self.fea_lambda = fea_lambda
            self.adj_lambda = adj_lambda
            self.gra_encoder = Gra_Encoder(in_node_fea, gat_channel, nhead, gat_dropout, concat, in_fea, list_gra_enc_hid, z_dim, gra_p_dropout)
            self.fea_decoder = Fea_Decoder(z_dim, list_fea_dec_hid, in_fea)
            self.fea_log_theta = nn.Parameter(torch.randn(in_fea))
            self.fea_gate_logits = nn.Parameter(torch.randn(in_fea))
            self.adj_decoder = Adj_Decoder(z_dim, list_adj_dec_hid, in_adj)
            self.adj_log_theta = nn.Parameter(torch.randn(in_adj))
            self.adj_gate_logits = nn.Parameter(torch.randn(in_adj))

        def model(self, x_gra):
            pyro.module('fea_decoder', self)
            pyro.module('adj_decoder', self)
            batch = x_gra.y.shape[0]
            with pyro.plate('data', batch):
                x_fea = x_gra.x_fea
                x_adj = x_gra.x_adj
                z_loc = x_fea.new_zeros((batch, self.z_dim))
                z_scale = x_fea.new_ones((batch, self.z_dim))
                with poutine.scale(scale=self.kl_beta):
                    z = pyro.sample('latent', dist.Normal(z_loc, z_scale).to_event(1))
                fea_px_scale = self.fea_decoder(z)
                fea_theta = torch.exp(self.fea_log_theta)
                fea_nb_logits = (fea_px_scale + 0.0001).log() - (fea_theta + 0.0001).log()
                fea_x_dist = dist.ZeroInflatedNegativeBinomial(total_count=fea_theta, logits=fea_nb_logits, gate_logits=self.fea_gate_logits)
                with poutine.scale(scale=self.fea_lambda):
                    fea_rx = pyro.sample('obs_fea', fea_x_dist.to_event(1), obs=x_fea)
                adj_px_scale = self.adj_decoder(z)
                adj_theta = torch.exp(self.adj_log_theta)
                adj_nb_logits = (adj_px_scale + 0.0001).log() - (adj_theta + 0.0001).log()
                adj_x_dist = dist.ZeroInflatedNegativeBinomial(total_count=adj_theta, logits=adj_nb_logits, gate_logits=self.adj_gate_logits)
                with poutine.scale(scale=self.adj_lambda):
                    adj_rx = pyro.sample('obs2', adj_x_dist.to_event(1), obs=x_adj)
                return (fea_rx, adj_rx)

        def guide(self, x_gra):
            pyro.module('gra_encoder', self)
            batch = x_gra.y.shape[0]
            with pyro.plate('data', batch):
                qz_m, qz_v = self.gra_encoder(x_gra, batch)
                with poutine.scale(scale=self.kl_beta):
                    rz = pyro.sample('latent', dist.Normal(qz_m, qz_v.sqrt()).to_event(1))
                return rz

        def getZ(self, x_gra):
            batch = x_gra.y.shape[0]
            z_mu, z_var = self.gra_encoder(x_gra, batch)
            return (z_mu, z_mu + z_var)

    def define_svi(in_node_fea, in_fea, in_adj, params, device):
        vae = VAE(in_node_fea, params['gat_channel'], params['nhead'], params['gat_dropout'], params['concat'], in_fea, params['list_gra_enc_hid'], params['z_dim'], params['gra_p_dropout'], params['list_fea_dec_hid'], in_adj, params['list_adj_dec_hid'], params['kl_beta'], params['fea_lambda'], params['adj_lambda'])
        vae.to(device)
        vae.train()
        adam_args = {'lr': params['lr'], 'foreach': False}
        optimizer = Adam(adam_args)
        elbo = Trace_ELBO(strict_enumeration_warning=False)
        svi = SVI(vae.model, vae.guide, optimizer, loss=elbo)
        return (vae, svi)
    return locals()
