import pathlib as pl
import numpy as np
import pandas as pd
import pcraster as pcr

def calculate_country_catchment_connections(countries: pcr.Field,
                                            catchments: pcr.Field,
                                            ldd: pcr.Field) -> pd.DataFrame:
    
    pit = pcr.pit(ldd) != 0
        
    catchments_array = pcr.pcr2numpy(catchments, 0)
    countries_array = pcr.pcr2numpy(countries, 0)
    pit_array = pcr.pcr2numpy(pit, 0)
    
    connections = pd.DataFrame({"country": countries_array[pit_array == 1].flatten(),
                                "catchment": catchments_array[pit_array == 1].flatten()})
    
    connections = connections.groupby(["catchment", "country"]).aggregate({"country": "size"})
    connections = connections.rename(columns = {"country": "pit_size"}).reset_index()
    max_pit_size_indices = connections.groupby(["catchment", "country"]).idxmax()["pit_size"]
    connections = connections.loc[max_pit_size_indices]
    
    catchments_unique, catchments_counts = np.unique(catchments_array,
                                                     return_counts = True)
    catchment_sizes = pd.DataFrame({"catchment": catchments_unique,
                                    "catchment_size": catchments_counts})
    
    connections = connections.merge(catchment_sizes,
                                    on = "catchment",
                                    how = "left")
    
    countries_unique, countries_counts = np.unique(countries_array,
                                                   return_counts = True)
    country_sizes = pd.DataFrame({"country": countries_unique,
                                  "country_size": countries_counts})
    
    connections = connections.merge(country_sizes,
                                    on = "country",
                                    how = "left")
    
    return connections


def calculate_country_catchment_replacements(connections: pd.DataFrame,
                                             threshold: int) -> pd.DataFrame:
    
    country_connections = connections.copy()
    
    # Loop over countries based on their total catchment size (smallest first)
    country_sizes = country_connections.groupby("country").aggregate({"catchment_size": "sum"})
    country_ids = country_sizes.sort_values("catchment_size").index
    
    replacements = {"from": [],
                    "to": []}
    
    if country_ids.size == 0:
        return pd.DataFrame(replacements)
    
    country_id = country_ids[-1]
    for country_id in country_ids:
        
        country_sel = country_connections["country"] == country_id
        if country_sel.sum() == 0:
            continue
        
        all_country_connections = country_connections.loc[country_sel]
        all_country_connections = all_country_connections.sort_values("catchment_size")
        all_country_connections_cumsum = all_country_connections["catchment_size"].cumsum()
        
        threshold_sel = all_country_connections_cumsum < threshold            
        if threshold_sel.sum() == 0:
            # Drop all delta-catchment connections that have been processed
            drop_sel = np.isin(country_connections["country"], country_id)
            drop_indices = country_connections.loc[drop_sel].index
            country_connections = country_connections.drop(drop_indices)
            continue
        
        all_catchment_ids = all_country_connections.loc[threshold_sel, "catchment"].to_numpy()
        catchment_id = all_catchment_ids[0]
        
        # Assign catchment id to all delta-connected catchments
        replacements["from"] += all_catchment_ids.tolist()
        replacements["to"] += [catchment_id] * all_catchment_ids.size
        
        # Drop all delta-catchment connections that have been processed
        drop_sel = np.isin(country_connections["country"], country_id)
        drop_indices = country_connections.loc[drop_sel].index
        country_connections = country_connections.drop(drop_indices)
    
    if country_connections.size > 0:
        print(country_connections)
        raise ValueError("Not all countrys have been processed")
    
    replacements = pd.DataFrame(replacements)
    return replacements


save_dir = pl.Path("./saves")
out_dir = pl.Path("./saves")

resolutions = [dir.stem for dir in save_dir.iterdir() if dir.is_dir()]
# resolutions = ["03min"]

for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    save_resolution_dir = save_dir / resolution
    out_resolution_dir = out_dir / resolution

    regions = [dir.stem for dir in save_resolution_dir.iterdir() if dir.is_dir()]
    # regions = ["Europe"]
    
    for region in regions:
        print(f"processing region {region}")
        
        save_region_dir = save_resolution_dir / region
        out_region_dir = out_resolution_dir / region
        
        country_catchments_file = out_region_dir / f"country_catchments_{region}.map"
        if country_catchments_file.exists():
            print("> exists")
            continue
        
        catchments_file = save_region_dir / f"delta_catchments_{region}.map"
        countries_file = save_region_dir / f"countries_{region}.map"
        ldd_file = save_region_dir / f"ldd_{region}.map"
        if not catchments_file.exists() or not countries_file.exists() and ldd_file.exists():
            print("> skipping")
            continue
        
        pcr.setclone(str(catchments_file))
        catchments = pcr.readmap(str(catchments_file))
        countries = pcr.readmap(str(countries_file))
        ldd = pcr.readmap(str(ldd_file))
        
        # Map country-catchment connections
        out_connections_file = out_region_dir / "country_catchment_connections.csv"
        if out_connections_file.exists():
            connections = pd.read_csv(out_connections_file)
        else:
            connections = calculate_country_catchment_connections(countries=countries,
                                                                catchments=catchments,
                                                                ldd=ldd)
            out_connections_file.parent.mkdir(parents = True, exist_ok = True)
            connections.to_csv(out_connections_file, index = False)
        
        threshold = connections["catchment_size"].max()
        
        # Map country-catchment replacements
        out_replacements_file = out_region_dir / "country_catchment_replacements.csv"
        if out_replacements_file.exists():
            replacements = pd.read_csv(out_replacements_file)
        else:
            replacements = calculate_country_catchment_replacements(connections=connections,
                                                                    threshold=threshold)
            out_replacements_file.parent.mkdir(parents = True, exist_ok = True)
            replacements.to_csv(out_replacements_file, index = False)
        
        # Assign catchment id to all delta-connected catchments
        catchments_array = pcr.pcr2numpy(catchments, 0)
        country_catchments_array = catchments_array.copy()
        
        replacements_remaining = replacements.copy()
        catchment_ids = replacements_remaining["to"].unique()
        for catchment_index, catchment_id in enumerate(catchment_ids):
            print(f"processing catchment {catchment_id} ({catchment_index}/{catchment_ids.size})")
            
            catchment_sel = replacements_remaining["to"] == catchment_id
            if catchment_sel.sum() == 0:
                continue
            
            all_catchment_ids = replacements_remaining.loc[catchment_sel, "from"]
            all_catchment_mask = np.isin(country_catchments_array, all_catchment_ids)
            country_catchments_array[all_catchment_mask] = catchment_id
            
            drop_indices = replacements_remaining.loc[catchment_sel].index
            replacements_remaining = replacements_remaining.drop(drop_indices)
        
        country_catchments = pcr.numpy2pcr(pcr.Nominal, country_catchments_array, 0)
        country_catchments = pcr.ifthenelse(catchments == 0, pcr.nominal(0), country_catchments)
        # pcr.aguila(country_catchments)
        # pcr.aguila(country_catchments != catchments)
        
        country_catchments_file.parent.mkdir(parents = True, exist_ok = True)
        pcr.report(country_catchments, str(country_catchments_file))
        