import pathlib as pl
import pcraster as pcr
import numpy as np
import pandas as pd


def calculate_delta_catchment_connections(deltas: pcr.Field,
                                          catchments: pcr.Field) -> pd.DataFrame:
    
    catchments_array = pcr.pcr2numpy(catchments, 0)
    deltas_array = pcr.pcr2numpy(deltas, 0)
    
    connections = pd.DataFrame({"delta": deltas_array.flatten(),
                                "catchment": catchments_array.flatten()})
    
    connections = connections.loc[connections["delta"] != 0]
    connections = connections.loc[connections["catchment"] != 0]
    connections = connections.drop_duplicates()
    
    catchments_unique, catchments_counts = np.unique(catchments_array,
                                                     return_counts = True)
    catchment_sizes = pd.DataFrame({"catchment": catchments_unique,
                                    "catchment_size": catchments_counts})
    
    connections = connections.merge(catchment_sizes,
                                    on = "catchment",
                                    how = "left")
    
    deltas_unique, deltas_counts = np.unique(deltas_array,
                                             return_counts = True)
    delta_sizes = pd.DataFrame({"delta": deltas_unique,
                                "delta_size": deltas_counts})
    
    connections = connections.merge(delta_sizes,
                                    on = "delta",
                                    how = "left")
    
    return connections


def calculate_delta_catchment_replacements(connections: pd.DataFrame,) -> pd.DataFrame:
    
    delta_connections = connections.copy()
    
    catchment_ids = connections["catchment"].unique()
    
    replacements = {"from" : [],
                    "to" : []}
    
    catchment_id = catchment_ids[-1]
    for catchment_id in catchment_ids:
        
        # Find all delta ids that are connected to the current catchment
        catchment_sel = delta_connections["catchment"] == catchment_id
        catchment_delta_ids = delta_connections.loc[catchment_sel, "delta"].unique()
        if catchment_delta_ids.size == 0:
            continue
        
        # Find all catchment ids that are connected to the current, and all subsequent, delta ids
        all_catchment_ids = np.array([catchment_id])
        all_catchment_delta_ids = catchment_delta_ids
        prev_all_catchment_ids = np.array([])
        prev_all_catchment_delta_ids = np.array([])
        
        # Loop until no new delta ids are found
        while prev_all_catchment_delta_ids.size != all_catchment_delta_ids.size and prev_all_catchment_ids.size != all_catchment_ids.size:
            prev_all_catchment_ids = all_catchment_ids.copy()
            prev_all_catchment_delta_ids = all_catchment_delta_ids.copy()
            
            all_catchment_sel = delta_connections["delta"].isin(all_catchment_delta_ids)
            all_catchment_ids = delta_connections.loc[all_catchment_sel, "catchment"].unique()
            all_catchment_delta_sel = delta_connections["catchment"].isin(all_catchment_ids)
            all_catchment_delta_ids = delta_connections.loc[all_catchment_delta_sel, "delta"].unique()

        # Assign catchment id to all delta-connected catchments
        replacements["from"] += all_catchment_ids.tolist()
        replacements["to"] += [catchment_id] * all_catchment_ids.size
        
        # Drop all delta-catchment connections that have been processed
        drop_sel = np.isin(delta_connections["catchment"], all_catchment_ids)
        drop_indices = delta_connections.loc[drop_sel].index
        delta_connections = delta_connections.drop(drop_indices)
        
    if delta_connections.size > 0:
        raise ValueError("Not all clumps have been processed")
        
    replacements = pd.DataFrame(replacements)
    return replacements


save_dir = pl.Path("./saves")
out_dir = pl.Path("./saves")

resolutions = [dir.stem for dir in save_dir.iterdir() if dir.is_dir()]

resolution = resolutions[0]
for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    save_resolution_dir = save_dir / resolution
    out_resolution_dir = out_dir / resolution
        
    delta_catchments_file = out_resolution_dir / "delta_catchments_global.map"
    if delta_catchments_file.exists():
        print(f"> exists")
        continue
    
    catchments_file = save_resolution_dir / "catchments_global.map"
    deltas_file = save_resolution_dir / "deltas_global.map"
    if not catchments_file.exists() or not deltas_file.exists():
        print(f"> skipping")
        continue
    
    pcr.setclone(str(catchments_file))
    catchments = pcr.readmap(str(catchments_file))
    deltas = pcr.readmap(str(deltas_file))
        
    # Map delta-catchment connections
    out_connections_file = out_resolution_dir / "delta_catchment_connections.csv"
    if out_connections_file.exists():
        connections = pd.read_csv(out_connections_file)
    else:
        connections = calculate_delta_catchment_connections(deltas, catchments)
        out_connections_file.parent.mkdir(parents = True, exist_ok = True)
        connections.to_csv(out_connections_file, index = False)
    
    # Map delta-catchment replacements
    out_replacements_file = out_resolution_dir / "delta_catchment_replacements.csv"
    if out_replacements_file.exists():
        replacements = pd.read_csv(out_replacements_file)
    else:
        replacements = calculate_delta_catchment_replacements(connections)
        out_replacements_file.parent.mkdir(parents = True, exist_ok = True)
        replacements.to_csv(out_replacements_file, index = False)   
    
    # Assign catchment id to all delta-connected catchments
    catchments_array = pcr.pcr2numpy(catchments, 0)
    delta_catchments_array = catchments_array.copy()
    
    replacements_remaining = replacements.copy()
    catchment_ids = replacements_remaining["to"].unique()
    for catchment_index, catchment_id in enumerate(catchment_ids):
        print(f"processing catchment {catchment_id} ({catchment_index}/{catchment_ids.size})")
        
        catchment_sel = replacements_remaining["to"] == catchment_id
        if catchment_sel.sum() == 0:
            continue
        
        all_catchment_ids = replacements_remaining.loc[catchment_sel, "from"]
        all_catchment_mask = np.isin(catchments_array, all_catchment_ids)
        delta_catchments_array[all_catchment_mask] = catchment_id
        
        drop_indices = replacements_remaining.loc[catchment_sel].index
        replacements_remaining = replacements_remaining.drop(drop_indices)
    
    delta_catchments = pcr.numpy2pcr(pcr.Nominal, delta_catchments_array, 0)
    # pcr.aguila(delta_catchments)
    # pcr.aguila(delta_catchments != catchments)
    
    delta_catchments_file.parent.mkdir(parents = True, exist_ok = True)
    pcr.report(delta_catchments, str(delta_catchments_file))