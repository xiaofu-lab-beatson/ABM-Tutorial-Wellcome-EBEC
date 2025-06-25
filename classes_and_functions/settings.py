"""_summary_

This Python script contains functions to get configurations and parameters required for setting up a simulation. 

"""

from typing import Tuple, Dict

def get_cell_configurations() -> Tuple[Dict]:
    """_summary_

    Returns:
        Tuple[Dict]: A Tuple of Dict variables, including
            site_types: Dict[int, str] maps integer to string names of site types
            sites_states: Dict[int, str] maps integer to string names of site states
            color_map: Dict[str, str] contains the colors corresponding to site types 
            markersize_map: Dict[str, float] contains the marker sizes corresponding to site types 
    """
    
    # configuration for cells
    site_types = {
        0: "CV", # central vein
        1: "PT", # portal triad 
        2: "HEP", # hepatocyte
        3: "NO", # not occupied
        4: "CC", # cancer cell 
        5: "ECM" 
    }

    sites_states = {
        0: "quiescent", 
        1: "proliferative", 
        2: "apoptotic", 
        # 3: "migratory"
    }

    color_map = {
        "CV" : "blue", 
        "PT": "red", 
        "HEP": "lightgreen",
        "NO" : "#EEEEEE",
        "CC" : "#525100",
        "ECM": "magenta"
    }
    markersize_map = {
        "CV": 2, "PT": 2, 
        "HEP": 0.75,
        "NO": 0.75,
        "CC": 1.25,
        "ECM": 1.25
    }
    
    return (site_types, sites_states, color_map, markersize_map)


def get_simulation_parameters(model_type="model_1") -> Dict[str, float]:
    """_summary_

    Args:
        model_type (str, optional): model type to implement in the simulation. Defaults to "model_1".

    Returns:
        Dict[str, float]: parameters to be used in the simulation, including 
            P_CC_GROW - probability of a cancer cell growing to an adjacent site
            P_HEP_DAMAGED - prbability of a healthy hepatocyte damaged by cancer cells to become apoptotic
            P_HEP_CLEARED - probability of an apoptotic hepatocyte becoming cleared
            P_CC_KILLED - probability of a cancer cell killed by an (implicitly implemented) cytotoxic T cell
    """
    
    # simulate tumour growth
    P_CC_GROW = 1           # probability of a cancer cell growing to an adjacent site
    P_HEP_DAMAGED = 0.5     # probability of a healthy hepatocyte damaged by cancer cells to become apoptotic
    P_HEP_CLEARED = 0.5     # probability of an apoptotic hepatocyte becoming cleared

    parameters = {}
    parameters['P_CC_GROW'] = P_CC_GROW
    parameters['P_HEP_DAMAGED'] = P_HEP_DAMAGED
    parameters['P_HEP_CLEARED'] = P_HEP_CLEARED
    
    if model_type=="model_1" or model_type=="model_2": # model_2 = model_1 + fibrosis at peri-central regions
        
        return parameters
    
    if model_type=="model_3": # model_3 = model_1 + implicit immune predation 
        P_CC_KILLED = 0.5       # probability of cancer cells being killed, by implicit immune predation 
        parameters['P_CC_KILLED'] = P_CC_KILLED
        
        return parameters
    
    if model_type=="model_4": # model_4 = model_1 + move or grow
        P_CC_GROW = 0.5       # probability of cancer cells growing; if not growing then move
        parameters['P_CC_GROW'] = P_CC_GROW
        
        return parameters
        