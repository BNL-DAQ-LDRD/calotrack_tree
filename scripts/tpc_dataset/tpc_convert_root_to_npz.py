"""
Convert root file to NumPy data
"""
from collections import defaultdict
from pathlib import Path
import yaml
import numpy as np
import pandas as pd
import uproot
from tqdm import tqdm

from torch.utils.data import Dataset


class TPCDataset(Dataset):
    """
    Load reconstructed clusters from TPC with 4 features: energy, x, y, z.
    Match each reco. cluster to a particle via the g4hit id.
    The segmentation target for each reco. cluster is the track id of the
    particle that leaves the trajectory. The regression target is the
    initial momentum of the particle (vx, vy, vz), vertex of the particle
    (vtx_x, vtx_y, vtx_z), and the energy of the particle.
    """
    def __init__(self, data_root):

        super().__init__()

        data_root = Path(data_root)

        # load data
        self.data = uproot.open(data_root/'data.root')["T;1"]
        self.num_events = len(self.data['nRecoClusters'].array(library='np'))

        print(f'load {self.num_events} from {str(data_root)}')

        self.reco_keys = ['E', 'x', 'y', 'z']
        self.particle_reg_keys = ['px', 'py', 'pz', 'vtx_x', 'vtx_y', 'vtx_z', 'energy']
        self.particle_seg_key = 'track_id'

        print('loading reconstructed clusters... ', end='')
        self.reco_cluster = self.load_branch(
            'reco_cluster',
            ['E', 'x', 'y', 'z', 'g4hit_id', 'detid']
        )
        print('Done')
        print('loading G4 hits... ', end='')
        self.g4hit = self.load_branch(
            'track_g4hit',
            ['id', 'trparticle_track_id']
        )
        print('Done')
        print('loading particles... ', end='')
        self.particle = self.load_branch(
            'particle',
            ['track_id', 'px', 'py', 'pz', 'vtx_x', 'vtx_y', 'vtx_z', 'energy']
        )
        print('Done')

    def __len__(self):
        return self.num_events

    def __getitem__(self, index):
        event = self.load_event(index)

        if isinstance(event, str):
            return event

        # NOTE: you may return the features and targets directly
        # without using the Data function from torch_geometric.
        # This function is just to make running GNN easier.
        return {'features': event[self.reco_keys].values,
                'seg_target': event[self.particle_seg_key].values,
                'reg_target': event[self.particle_reg_keys].values}

    def load_branch(self, branch, keys):
        """
        Load a branch (calorimeter hits, tracking reco. cluster, etc) with a set of keys
        """
        return {key: self.data[f'{branch}_{key}'].array(library='np') for key in keys}

    @staticmethod
    def __load_event(branch_map, event_id):
        event_data = {key: val[event_id] for key, val in branch_map.items()}
        return pd.DataFrame(data=event_data)

    def load_event(self, event_id):
        """
        Load an event
        """
        reco_cluster = self.__load_event(self.reco_cluster, event_id)
        tpc_reco_cluster = reco_cluster[ reco_cluster.detid == 2 ]
        if len(tpc_reco_cluster) == 0:
            return 'empty_tpc'

        reco_bad = tpc_reco_cluster[self.reco_keys].isnull().values.any()
        if reco_bad:
            return 'missing_reco_data'

        g4hit = self.__load_event(self.g4hit, event_id)
        particle = self.__load_event(self.particle, event_id)

        tpc_reco_cluster_with_track_id = pd.merge(tpc_reco_cluster,
                                                  g4hit,
                                                  left_on='g4hit_id',
                                                  right_on='id',
                                                  how='left')

        missing_g4hit = tpc_reco_cluster_with_track_id.trparticle_track_id.isna().sum()
        if missing_g4hit > 0:
            return 'missing_g4hits'

        # match track_id
        tpc_reco_cluster_with_particle = pd.merge(tpc_reco_cluster_with_track_id,
                                                  particle,
                                                  left_on ='trparticle_track_id',
                                                  right_on='track_id',
                                                  how='left')

        particle_bad = tpc_reco_cluster_with_particle[self.particle_reg_keys].isnull().values.any()
        if particle_bad:
            return 'missing_particle_data'

        return tpc_reco_cluster_with_particle


if __name__ == '__main__':

    # Set parameters here.
    DATA_ROOT = '/home/sphenix_fm/100kevt/'
    SAVE_ROOT = '/home/sphenix_fm/100kevt_npy'
    TRAIN_FRACTION = .8

    # Create path to save the NumPy files.
    SAVE_ROOT = Path(SAVE_ROOT)
    if not SAVE_ROOT.exists():
        SAVE_ROOT.mkdir(parents=True)

    # Create dataset, load events,
    # filter out problematic ones, and save the good ones.
    dataset = TPCDataset(DATA_ROOT)

    problematic_data = defaultdict(list)

    pbar = tqdm(range(len(dataset)), total=len(dataset))
    for event_id in pbar:

        data = dataset[event_id]

        if isinstance(data, str):
            problematic_data[data].append(event_id)
            continue

        np.savez_compressed(SAVE_ROOT/f'event_{event_id}', **data)

        pbar.set_postfix({key: len(val)
                          for key, val in problematic_data.items()})
        pbar.update(1)

    with open(SAVE_ROOT/'problematic_events.yaml', 'w', encoding='UTF-8') as handle:
        yaml.dump(dict(problematic_data), handle)

    # Create train-test split and save the lists of train/test events.
    fnames = list(SAVE_ROOT.glob('*npz'))
    event_ids = sorted([int(fname.stem.split('_')[-1]) for fname in fnames])
    print(f'Numbers of good events: {len(event_ids)}')

    num_train = int(len(event_ids) * TRAIN_FRACTION)

    np.random.shuffle(event_ids)

    train_event_ids = sorted(event_ids[:num_train])
    test_event_ids = sorted(event_ids[num_train:])

    np.savetxt(SAVE_ROOT/'train.txt', train_event_ids, fmt='%d', newline='\n')
    np.savetxt(SAVE_ROOT/'test.txt', test_event_ids, fmt='%d', newline='\n')
