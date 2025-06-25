"""_summary_

This Python script contains functions to simulate tumour growth and liver damage.  

"""

import pandas as pd
from classes_and_functions.cell_classes import CancerCell, Hepatocyte
import numpy as np

from typing import Dict, Tuple

def update_cell_states(
    cell_dictionaries: Dict[str, Dict],
    lattice: pd.DataFrame,
    parameters: Dict[str, float],
    CancerCell: CancerCell, 
    Hepatocyte: Hepatocyte,
    model_type: str="model_1"
) -> Tuple[Dict, pd.DataFrame]:
    """_summary_
    
    This function simulates cancer cell proliferation and migration, hepatocyte death, with cell_dictionaries and lattice updated accordingly.

    Args:
        cell_dictionaries (Dict[str, Dict]): a Dict of Dict containing the CancerCell and Hepatocyte objects
        lattice (pd.DataFrame): a DataFrame containing information following initialisation of CancerCell objects
        parameters (Dict[str, float]): parameters to be used in the simulation (see settings.py)
        CancerCell (CancerCell): the CancerCell class used for creating new CancerCell objects
        Hepatocyte (Hepatocyte): the Hepatocyte class used for creating new Hepatocyte objects
        model_type (str, optional): model type to implement in the simulation. Defaults to "model_1".

    Returns:
        Tuple[Dict, pd.DataFrame]: an updated Dict of Dict containing the CancerCell and Hepatocyte objects (with updated attributes; newly created and deleted objects)
    """
    
    
    p_cc_grow = parameters['P_CC_GROW']
    p_hep_damaged = parameters['P_HEP_DAMAGED']
    p_hep_cleared = parameters['P_HEP_CLEARED']

    dict_of_cancer_cells, dict_of_hepatocytes = \
        cell_dictionaries['CancerCell'], cell_dictionaries['Hepatocyte']
        
    new_dict_of_cancer_cells = dict_of_cancer_cells.copy()
    new_dict_of_hepatocytes  = dict_of_hepatocytes.copy()
        
    # ASSUME for now: only those hepatocytes adjacent to proliferative cancer cells need to be processed
    list_of_hep_ids_to_process = []
        
    # update states of cancer cells  
    for cancer_cell_id, cancer_cell in dict_of_cancer_cells.items():
        # cancer_cell_attributes = cancer_cell.attributes
        cancer_cell_attributes = cancer_cell.get_attributes()
        
        # cancer_cell_id = cancer_cell_attributes['cell_id']
        cancer_cell_site_id = cancer_cell_attributes['site_id']
        cancer_cell_state = cancer_cell_attributes['cell_state']
        cancer_cell_position = cancer_cell_attributes['cell_position']
        
        # get adjacent site ids
        adjacent_site_ids_str = lattice.loc[
            lattice.site_id==cancer_cell_site_id,
            "adjacent_site_ids_str"
        ].values[0]
        adjacent_site_ids = [
            int(str) for str in adjacent_site_ids_str.split(',')
        ]
        
        # proliferative cancer cells grow
        if cancer_cell_state == 1:
            
            # interate over adjacent sites
            for adjacent_site_id in adjacent_site_ids:
                adjacent_site_type, adjacent_site_x, adjacent_site_y, adjacent_site_cell_id = lattice.loc[
                    lattice.site_id==adjacent_site_id, 
                    ["site_type", 'x', 'y', 'cell_id']
                ].values[0]
                
                # grow into NO "Not Occupied" adjacent site
                if adjacent_site_type == 3: # site type = "NO"
                    if np.random.random() < p_cc_grow: # grow to the adjacent site
                        
                        # add a new CancerCell
                        new_cancer_cell_id = int(lattice.cell_id.max() + 1)
                        new_cancer_cell_site_id = adjacent_site_id
                        
                        new_cancer_cell_xy = (adjacent_site_x, adjacent_site_y)

                        new_cancer_cell_attributes = {
                            "cell_id": new_cancer_cell_id, "site_id": new_cancer_cell_site_id,
                            "cell_position": new_cancer_cell_xy,
                            "cell_state": 1 # proliferative 
                        }

                        ## [1] create cancer cell object
                        new_cancer_cell = CancerCell(cell_attributes=new_cancer_cell_attributes)
                        
                        new_dict_of_cancer_cells |= {new_cancer_cell_id: new_cancer_cell}
                        
                        ## [2] update lattice
                        lattice.loc[
                            lattice.site_id==adjacent_site_id, 
                            ["site_type", "cell_id"]
                        ] = (4, new_cancer_cell_id) # sitetype = cancer cell
                        
                    else:
                        if model_type=="model_4": # move to the adjacent site
                            
                            # update the attributes
                            # ... site id
                            # ... xy position
                            cancer_cell_site_id_new = adjacent_site_id
                            cancer_cell_xy_new = (adjacent_site_x, adjacent_site_y)
                            updated_cancer_cell_attributes = {
                                "cell_id": cancer_cell_id, "site_id": cancer_cell_site_id_new,
                                "cell_position": cancer_cell_xy_new,
                                "cell_state": 1 # proliferative 
                            }
                            new_dict_of_cancer_cells[cancer_cell_id].set_attributes(updated_cancer_cell_attributes)
                            
                            # update the lattice site
                            # ... previous site to be emptied 
                            # ... new site to be filled 
                            lattice.loc[
                                lattice.site_id==cancer_cell_site_id, 
                                ["site_type", "cell_id"]
                            ] = (3, np.nan) # sitetype = not occupied 
                            lattice.loc[
                                lattice.site_id==cancer_cell_site_id_new, 
                                ["site_type", "cell_id"]
                            ] = (4, cancer_cell_id) # sitetype = cancer cell
                        
                elif adjacent_site_type == 2: # site type = "HEP"
                    list_of_hep_ids_to_process.append(adjacent_site_cell_id)
                
                # (more conditions...)
        
    # update states of hepatocytes
    # list_of_hep_ids_to_process = list(dict_of_hepatocytes.keys())
    
    total_number_of_hepatocytes_to_process = len(list_of_hep_ids_to_process)
    number_processed = 0
    # for hep_id, hep in dict_of_hepatocytes.items():
    for hep_id in list_of_hep_ids_to_process:
        
        if hep_id not in dict_of_hepatocytes.keys() or hep_id not in new_dict_of_hepatocytes.keys():
            continue
        
        hep = dict_of_hepatocytes[hep_id]
        
        number_processed += 1
        # if number_processed % (total_number_of_hepatocytes_to_process / 20) == 0:
        #     print(f"> {number_processed / total_number_of_hepatocytes_to_process:.2%} hepatocytes processed...") 
        
        # hep_attributes = hep.attributes
        hep_attributes = hep.get_attributes()
        
        hep_xy = hep_attributes['cell_position']
        hep_id = hep_attributes['cell_id']
        hep_state = hep_attributes['cell_state']
        hep_site_id = hep_attributes['site_id']
        
        # get adjacent site ids
        adjacent_site_ids_str = lattice.loc[
            lattice.site_id==hep_site_id,
            "adjacent_site_ids_str"
        ].values[0]
        adjacent_site_ids = [
            int(str) for str in adjacent_site_ids_str.split(',')
        ]
        
        # quiescent hepatocytes change cell states
        if hep_state == 0:
        
            # interate over adjacent sites
            for adjacent_site_id in adjacent_site_ids:
                adjacent_site_type, adjacent_site_x, adjacent_site_y = lattice.loc[
                    lattice.site_id==adjacent_site_id, 
                    ["site_type", 'x', 'y']
                ].values[0]
                
                # turn into apoptotic state if a proliferative cancer cell is adjacent
                if adjacent_site_type == 4: # site type = "CC"
                    if np.random.random() < p_hep_damaged:
                        updated_hep_attributes = hep_attributes.copy()
                        updated_hep_attributes['cell_state'] = 2
                        hep.attributes = updated_hep_attributes
                        new_dict_of_hepatocytes[hep_id] = hep
                    
        # apoptotic hepatocytes get cleared
        elif hep_state == 2:
            if np.random.random() < p_hep_cleared: # get cleared
                
                # [1] delete this hepatocyte 
                del new_dict_of_hepatocytes[hep_id]
                
                # [2] update lattice site type
                lattice.loc[
                    lattice.site_id==hep_site_id, 
                    ["site_type", "cell_id"]
                ] = (3, np.nan) # change to Not Occupied
                
            else: # not get cleared
                
                if model_type=="model_1": # fibrosis not considered in this model
                    pass
                
                elif model_type=="model_2": # fibrosis is considered; for simplicity, for apoptotic hepatocytes not cleared, they turn ECM deposited
                    
                    if lattice.loc[
                        lattice.site_id==hep_site_id
                    ].zonation_type.values[0]=="peri-central": # check if this hepatocyte is located in peri-zonal zonation
                    
                        # [1] delete this hepatocyte 
                        del new_dict_of_hepatocytes[hep_id]
                        
                        # [2] update lattice site type
                        lattice.loc[
                            lattice.site_id==hep_site_id, 
                            ["site_type", "cell_id"]
                        ] = (5, np.nan) # change to ECM
                    
                    # (to be considered) whether or not to introduce ECM as a class

            # (more actions)

        # (more conditions)
                
    new_cell_dictionaries = {
        "CancerCell": new_dict_of_cancer_cells, "Hepatocyte": new_dict_of_hepatocytes
    }
    
    return new_cell_dictionaries, lattice

def implicit_immune_predation(
    cell_dictionaries: Dict[str, Dict],
    lattice: pd.DataFrame,
    parameters: Dict[str, float],
    model_type: str="model_3"
) -> Tuple[Dict, pd.DataFrame]:
    """_summary_
    
    This function simulates cancer cell death (implicitly killed by cytotoxic immune cells), with cell_dictionaries and lattice updated accordingly

    Args:
        cell_dictionaries (Dict[str, Dict]): a Dict of Dict containing the CancerCell and Hepatocyte objects
        lattice (pd.DataFrame): a DataFrame containing information following initialisation of CancerCell objects
        parameters (Dict[str, float]): parameters to be used in the simulation (see settings.py)
        model_type (str, optional): model type to implement in the simulation. Defaults to "model_3".

    Returns:
        Tuple[Dict, pd.DataFrame]: an updated Dict of Dict containing the CancerCell and Hepatocyte objects (with updated attributes; newly created and deleted objects)
    """
    
    if model_type not in ["model_3"]:
        print("model type is wrong! this function shouldn't be called!")
        return None
    
    p_cc_killed = parameters['P_CC_KILLED']
    
    dict_of_cancer_cells, dict_of_hepatocytes = \
        cell_dictionaries['CancerCell'], cell_dictionaries['Hepatocyte']
    
    new_dict_of_cancer_cells = dict_of_cancer_cells.copy()
    new_dict_of_hepatocytes  = dict_of_hepatocytes.copy()
    
    n_lattice_sites = lattice.shape[0]
    n_lattice_sites_by_tumour = lattice.loc[lattice.site_type==4].shape[0]
    lattice_tumour = lattice.loc[lattice.site_type==4].copy()
    
    # randomly sample K sites as being immune infiltrated/attacked
    lattice_site_ids_immune_attack = np.random.choice(lattice.site_id.values, size=n_lattice_sites_by_tumour, replace=False)
    
    # get the list of cell ids under attack that are tumour
    cancer_cell_ids_to_be_killed = lattice_tumour.loc[
        lattice_tumour.site_id.isin(lattice_site_ids_immune_attack)
    ].cell_id.values
    
    # kill cancer cells
    for cancer_cell_id in cancer_cell_ids_to_be_killed:
        cancer_cell = dict_of_cancer_cells[cancer_cell_id]
        cancer_cell_attributes = cancer_cell.attributes
        cancer_cell_site_id = cancer_cell_attributes['site_id']
        cancer_cell_state = cancer_cell_attributes['cell_state']
        cancer_cell_position = cancer_cell_attributes['cell_position']
        
        # killed with a probability
        if np.random.random() < p_cc_killed:
            
            # [1] delete this hepatocyte 
            del new_dict_of_cancer_cells[cancer_cell_id]
            
            # [2] update lattice site type
            lattice.loc[
                lattice.site_id==cancer_cell_site_id, 
                ["site_type", "cell_id"]
            ] = (3, np.nan) # change to Not Occupied
              
    new_cell_dictionaries = {
        "CancerCell": new_dict_of_cancer_cells, "Hepatocyte": new_dict_of_hepatocytes
    }
    
    return new_cell_dictionaries, lattice