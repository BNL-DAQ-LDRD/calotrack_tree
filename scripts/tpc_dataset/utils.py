import uproot
import numpy as np
import pandas as pd
import yaml

DATAROOT = "/home/yhuang2/PROJs/calotrack_tree_data/pp100evt_2.root"
META = '/home/yhuang2/PROJs/calotrack_tree/scripts/meta.yaml'
PID = '/home/yhuang2/PROJs/calotrack_tree/scripts/pid.yaml'
with open(PID, 'r') as fh:
    PID_MAP = yaml.safe_load(fh)


def get_all_branches():
    with open(META, 'r') as fh:
        meta = yaml.safe_load(fh)
        return meta.keys()

def load_branch(branch, event_id):
    # load keys in a branch from the meta data
    with open(META, 'r') as fh:
        meta = yaml.safe_load(fh)
    keys = meta[branch]

    # construct the key name for the count
    count_key = 'n' + ''.join(branch.replace('_', ' ').title().split()) + 's'

    # load data
    data = uproot.open(DATAROOT)["T;1"]
    dataframe_data = {}
    for key in keys:
        if key == count_key:
            num_records = data[key].array(library='np')[event_id]
        else:
            if key.startswith(branch):
                column_name = key.removeprefix(f'{branch}_')
            elif key.startswith(branch.title()):
                column_name = key.removeprefix(f'{branch.title()}_')
            else:
                raise KeyError('unkonwn key {key}')
            
            dataframe_data[column_name] = data[key].array(library='np')[event_id]

    # construct data frame and make sure 
    # the number of record mataches
    dataframe = pd.DataFrame(data=dataframe_data)
    
    assert len(dataframe) == num_records, \
        f"unmatched number of records ({num_records} != {len(dataframe)}) "

    print(f'Loaded {num_records} from branch {branch} of event {event_id}.')
    
    return dataframe