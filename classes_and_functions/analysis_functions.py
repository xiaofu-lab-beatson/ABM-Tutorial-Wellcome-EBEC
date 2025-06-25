"""_summary_

This Python script contains a function to detect tumours using the DBSCAN algorithm 

"""

import pandas as pd
from sklearn.cluster import DBSCAN
from typing import Tuple 

def get_tumour_sizes(
    tumour_t: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """_summary_
    
    This function detects tumours using the DBSCAN algorithm.

    Args:
        tumour_t (pd.DataFrame): a DataFrame containing information about site_id, x, y, site_type, cell_id, adjacent_site_ids_str, zonation_type

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: A Tuple of DataFrame objects including
            tumour_t_sizes - a DataFrame object containing information about DBSCAN cluster ids and sizes
            tumour_t_labelled - a DataFrame object containing information as in tumour_t and further including annotated DBSCAN cluster ids that cells belong to
    """
    
    
    dbs = DBSCAN(eps=1.05, min_samples=1) # spacing = 1 in simulation
    dbs.fit(tumour_t[['x','y']].values)
    tumour_t_labelled = tumour_t.copy()
    tumour_t_labelled['label'] = dbs.labels_
    tumour_t_sizes = tumour_t_labelled.groupby('label', as_index=False).agg({'cell_id':'count'})
    tumour_t_sizes.rename(columns={'cell_id': 'size'}, inplace=True)
    
    return tumour_t_sizes, tumour_t_labelled
    