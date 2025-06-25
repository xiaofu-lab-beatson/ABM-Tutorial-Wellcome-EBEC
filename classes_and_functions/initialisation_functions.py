"""_summary_

This Python script contains function to initialise CancerCell objects nad update lattice accordingly

"""

from classes_and_functions.cell_classes import CancerCell, Hepatocyte
import numpy as np
import pandas as pd

from typing import Dict, Tuple

#### This function without docstrings is NOT used for the tutorial ####
def init_lattice_in_simulation(lattice: pd.DataFrame, lattice_settings: Dict) -> pd.DataFrame:
    lattice_in_simulation = lattice[
        ["site_id", "x", "y", "site_type"]
    ].copy()

    lattice_in_simulation["cell_id"] = lattice_in_simulation["site_id"].copy()

    lattice_in_simulation["cell_id"] = lattice_in_simulation["cell_id"].astype(np.uint32)
    lattice_in_simulation["site_id"] = lattice_in_simulation["site_id"].astype(np.uint16)
    lattice_in_simulation["site_type"] = lattice_in_simulation["site_type"].astype(np.uint8)
    
    # get adjacent site ids for all lattice sites - it feels that this procedure can be optimised
    def get_adjacent_site_xys(xy_array: np.ndarray, spacing: float) -> np.ndarray:

        radians = np.linspace(0, 2*np.pi *5/6, 6)
        adjacent_xy_arrays = np.array([
                (xy_array[0] + spacing * np.cos(radian), xy_array[1] + spacing * np.sin(radian))
                for radian in radians
        ])

        return adjacent_xy_arrays

    def get_adjacent_site_ids(adjacent_xy_arrays, lattice, threshold):
        
        adjacent_site_ids = []
        for xy_array in adjacent_xy_arrays:
            
            lattice_xy_arrays = lattice[['x','y']].values
            lattice_site_ids  = lattice.site_id.values
            distances = np.linalg.norm(lattice_xy_arrays-xy_array, axis=1)
            which_ones = np.where(distances<0.1*lattice_settings['spacing'])[0]
            if which_ones:
                adjacent_site_ids.append(str(lattice_site_ids[which_ones[0]]))
            
        if len(adjacent_site_ids):
            adjacent_site_ids_str = ','.join(adjacent_site_ids)
        else:
            adjacent_site_ids_str = ""
        
        return adjacent_site_ids_str

    def annotate_lattice_with_adjacent_sites(lattice):
        
        list_of_adjacent_site_ids_strs = []
        for site_id, x, y, site_type in lattice[['site_id', 'x', 'y', 'site_type']].values:
            
            if site_id % 10000 == 0 or site_id == lattice.shape[0]-1:
                print(f"> {site_id/lattice.shape[0]:.2%} processed...")
            
            if site_type==2:
                adjacent_xy_arrays = get_adjacent_site_xys(
                    xy_array=np.array((x,y)), spacing=lattice_settings['spacing']
                )
                adjacent_site_ids_str = get_adjacent_site_ids(
                    adjacent_xy_arrays=adjacent_xy_arrays,
                    lattice=lattice.loc[lattice.site_type.isin([0,1,2])],
                    threshold=0.1*lattice_settings['spacing']
                )
                list_of_adjacent_site_ids_strs.append(adjacent_site_ids_str)
                
            else:
                list_of_adjacent_site_ids_strs.append('n/a')
                
        lattice["adjacent_site_ids_str"] = list_of_adjacent_site_ids_strs
        
        return lattice

    def annotate_lattice_with_zonation(lattice):
        
        cv_xy_arrays = lattice.loc[lattice.site_type==0, ['x', 'y']].values
        list_of_zonation_types = []
        
        for site_id, x, y, site_type in lattice[['site_id', 'x', 'y', 'site_type']].values:
            
            if site_type==2:
                xy_array=np.array((x,y))
                distances = np.linalg.norm(cv_xy_arrays-xy_array, axis=1)
                if distances.min() <= 5**lattice_settings['spacing']:
                    zonation_type = "peri-central"
                else:
                    zonation_type = "other"
            else:
                zonation_type = "n/a"
                
            list_of_zonation_types.append(zonation_type)
            
        lattice["zonation_type"] = list_of_zonation_types
        
        return lattice

    print("> annotate adjacent sites ...")
    lattice_in_simulation = annotate_lattice_with_adjacent_sites(lattice=lattice_in_simulation)
    print("> annotate zonation type ...")
    lattice_in_simulation = annotate_lattice_with_zonation(lattice=lattice_in_simulation)
    
    return lattice_in_simulation

def init_cell_dictionaries(
    lattice: pd.DataFrame, 
    n_cancer_cells_init: int, 
    CancerCell: CancerCell, 
    Hepatocyte: Hepatocyte
) -> Tuple[Dict, pd.DataFrame]:
    """_summary_
    
    This function initialises cancer cells in the lattice.

    Args:
        lattice (pd.DataFrame): a DataFrame containing information about site_id, x, y, site_type, cell_id, adjacent_site_ids_str, zonation_type
        n_cancer_cells_init (int): the number of cancer cells to initialise in the lattice
        CancerCell (CancerCell): the CancerCell class used for creating new CancerCell objects
        Hepatocyte (Hepatocyte): the Hepatocyte class used for creating new Hepatocyte objects

    Returns:
        Tuple[Dict, pd.DataFrame]: a Tuple of objects including
            cell_dictionaries - a Dict of Dict containing the CancerCell and Hepatocyte objects
            lattice - a DataFrame containing information following initialisation of CancerCell objects
    """
    
    # initial configuration of hepatocytes
    dict_of_hepatocytes  = {} # id : Hepatocyte()
    for site_id, x, y, site_type, cell_id, adjacent_site_ids_str, zonation_type in lattice.values:
        
        if site_type == 2:
        
            cell_attributes = {
                "cell_id": cell_id, "site_id": site_id,
                "cell_position": (x,y),
                "cell_state": 0 
            }
            hep = Hepatocyte(cell_attributes=cell_attributes)
            
            dict_of_hepatocytes.update({site_id: hep})

    # introduce the first cancer cell
    dict_of_cancer_cells = {} # id : CancerCell()

    print("BEFORE: total number of hepatocytes: %d " % len(dict_of_hepatocytes))

    cancer_cell_id = int(lattice.cell_id.max() + 1)
    site_ids_to_sample = [
        site_id
        for site_id in dict_of_hepatocytes.keys()
        # if site_id < 100
    ]
    # print(site_ids_to_sample)

    cancer_cell_site_ids = np.random.choice(site_ids_to_sample, size=n_cancer_cells_init, replace=False)
    print(f"> selecting {cancer_cell_site_ids.size} sites to create the first CancerCell objects ")

    for cancer_cell_site_id in cancer_cell_site_ids:

        cancer_cell_xy = dict_of_hepatocytes[cancer_cell_site_id].attributes['cell_position']

        cancer_cell_attributes = {
            "cell_id": cancer_cell_id, "site_id": cancer_cell_site_id,
            "cell_position": cancer_cell_xy,
            "cell_state": 1 # proliferative 
        }

        ## [1] create cancer cell object
        cancer_cell = CancerCell(cell_attributes=cancer_cell_attributes)
        dict_of_cancer_cells.update({cancer_cell_id: cancer_cell})

        ## [2] delete hepatocyte
        del dict_of_hepatocytes[cancer_cell_site_id]

        ## [3] update lattice information
        # print(lattice.loc[lattice.site_id==cancer_cell_site_id])
        lattice.loc[
            lattice.site_id==cancer_cell_site_id,
            ["site_type", "cell_id"]
        ] = (4, cancer_cell_id) # double check cell type corresponds to cancer cell
        # print(lattice.loc[lattice.site_id==cancer_cell_site_id])

        cancer_cell_id += 1
        
    print("AFTER : total number of hepatocytes : %d " % len(dict_of_hepatocytes))
    print("AFTER : total number of cancer cells: %d " % len(dict_of_cancer_cells))

    cell_dictionaries = {
        "CancerCell": dict_of_cancer_cells,
        "Hepatocyte": dict_of_hepatocytes
    }

    return cell_dictionaries, lattice