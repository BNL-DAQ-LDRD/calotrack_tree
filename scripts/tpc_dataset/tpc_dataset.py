from pathlib import Path
import yaml
import pandas as pd
import uproot

import torch
from torch.utils.data import Dataset
from torch_geometric.data import Data


class TPCDataset(Dataset):
    def __init__(self, data_root):
        data_root = Path(data_root)

        # load meta
        with open(data_root/'meta.yaml', 'r') as fh:
            self.meta = yaml.safe_load(fh)

        # load data
        self.data = uproot.open(data_root/'data.root')["T;1"]
        self.num_events = len(self.data['nRecoClusters'].array(library='np'))

        print(f'load {self.num_events} from {str(data_root)}')

        self.reco_cols = ['E', 'x', 'y', 'z']
        self.particle_reg_cols = ['px', 'py', 'pz', 'vtx_x', 'vtx_y', 'vtx_z', 'energy']
        self.particle_seg_col = 'track_id'

    def __len__(self):
        return self.num_events

    def __getitem__(self, index):
        reco_cluster = self.load_branch('reco_cluster', index)
        particle = self.load_branch('particle', index)
        g4hit = self.load_branch('track_g4hit', index)

        # reconstructed clusters in TPC
        tpc_reco_cluster = reco_cluster[ reco_cluster.detid == 2 ]
        if len(tpc_reco_cluster) == 0:
            return {'features': Data(x=torch.randn(0, len(self.reco_cols))),
                    'seg_target': Data(x=torch.randn(0)),
                    'reg_target': Data(x=torch.randn(0, len(self.particle_reg_cols)))}

        # match g4hits
        tpc_reco_cluster_with_track_id = pd.merge(tpc_reco_cluster,
                                                  g4hit[['trparticle_track_id', 'id']],
                                                  left_on='g4hit_id',
                                                  right_on='id',
                                                  how='left')

        missing_g4hit = tpc_reco_cluster_with_track_id.trparticle_track_id.isna().any()
        assert not missing_g4hit, \
            'There exist reconstructed clusters with missing particle track'

        # match track_id
        tpc_reco_cluster_with_particle = pd.merge(tpc_reco_cluster_with_track_id,
                                                  particle,
                                                  left_on ='trparticle_track_id',
                                                  right_on='track_id',
                                                  how='left')

        features = torch.tensor(tpc_reco_cluster_with_particle[self.reco_cols].values)
        seg_target = torch.tensor(tpc_reco_cluster_with_particle[self.particle_seg_col].values)
        reg_target = torch.tensor(tpc_reco_cluster_with_particle[self.particle_reg_cols].values)

        return {'features': Data(x=features),
                'seg_target': Data(x=seg_target),
                'reg_target': Data(x=reg_target)}

    def load_branch(self, branch, event_id):

        # construct the key name for the count
        count_key = 'n' + ''.join(branch.replace('_', ' ').title().split()) + 's'

        # load data
        keys = self.meta[branch]
        dataframe_data = {}
        for key in keys:
            if key == count_key:
                num_records = self.data[key].array(library='np')[event_id]
            else:
                if key.startswith(branch):
                    column_name = key.removeprefix(f'{branch}_')
                elif key.startswith(branch.title()):
                    column_name = key.removeprefix(f'{branch.title()}_')
                else:
                    raise KeyError('unkonwn key {key}')

                dataframe_data[column_name] = self.data[key].array(library='np')[event_id]

        # construct data frame and make sure
        # the number of record mataches
        dataframe = pd.DataFrame(data=dataframe_data)

        assert len(dataframe) == num_records, \
            f"unmatched number of records ({num_records} != {len(dataframe)}) "

        print(f'Loaded {num_records} from branch {branch} of event {event_id}.')

        return dataframe
