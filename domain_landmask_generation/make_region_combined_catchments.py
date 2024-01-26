from typing import Any
import pathlib as pl
import numpy as np
import pandas as pd
import pcraster as pcr


def aquire_map_attributes(calc_lonlat: bool = False) -> dict[str, float]:
    resolution = pcr.clone().cellSize()
    north = pcr.clone().north()
    west = pcr.clone().west()
    ncols = pcr.clone().nrCols()
    nrows = pcr.clone().nrRows()
    
    attributes = {"resolution": resolution,
                  "north": north,
                  "west": west,
                  "ncols": ncols,
                  "nrows": nrows}
    
    if calc_lonlat:
        lons = np.arange(west + resolution / 2, west + resolution * ncols, resolution)
        lats = np.arange(north - resolution / 2, north - resolution * nrows, -resolution)
        attributes["lons"] = lons
        attributes["lats"] = lats
    
    return attributes


def calculate_catchment_details(catchments: pcr.Field,
                              ldd: pcr.Field) -> pd.DataFrame:
    
    attributes = aquire_map_attributes(calc_lonlat = True)
    
    pit = pcr.pit(ldd) != 0
    
    catchments_array = pcr.pcr2numpy(catchments, 0)
    pit_array = pcr.pcr2numpy(pit, 0)
    lons_array = np.tile(attributes["lons"], (attributes["nrows"], 1))
    lats_array = np.tile(attributes["lats"], (attributes["ncols"], 1)).T
    
    catchments_unique, catchments_counts = np.unique(catchments_array,
                                                     return_counts = True)
    catchment_sizes = pd.DataFrame({"catchment": catchments_unique,
                                    "size": catchments_counts})
    catchment_sizes = catchment_sizes.loc[catchment_sizes["catchment"] != 0]
    
    catchment_pits = pd.DataFrame({"catchment": catchments_array[pit_array == 1],
                                   "lon": lons_array[pit_array == 1],
                                   "lat": lats_array[pit_array == 1]})
    catchment_pits = catchment_pits.loc[catchment_pits["catchment"] != 0]
    
    catchment_details = catchment_sizes.merge(catchment_pits,
                                            on = "catchment",
                                            how = "left")
    
    return catchment_details


def calculate_combined_catchment_replacements(catchment_details: pd.DataFrame,
                                              threshold: int,
                                              expansion_distance: float) -> pd.DataFrame:
    
    catchment_details = catchment_details[catchment_details["size"] < threshold]
    catchment_details = catchment_details.sort_values(by = "size")
    catchment_ids = catchment_details.index.unique()
    
    catchment_details_remaining = catchment_details.copy()
    
    replacements = {"from": [],
                    "to": []}
    
    if catchment_ids.size == 0:
        replacements = pd.DataFrame(replacements)
        return replacements
    
    catchment_id = catchment_ids[0]
    for catchment_id in catchment_ids:
        catchment_detail = catchment_details_remaining.loc[catchment_id]
        if catchment_detail.index.size == 0:
            continue
        
        size = catchment_detail["size"]
        lon = catchment_detail["lon"]
        lat = catchment_detail["lat"]
        if size > threshold:
            continue
        
        closest = catchment_details_remaining.loc[catchment_details_remaining.index != catchment_id].copy()
        
        distances = np.sqrt((closest["lon"] - lon)**2 + (closest["lat"] - lat)**2)
        closest["distance"] = distances
        closest = closest.loc[closest["distance"] < expansion_distance]
        if closest.index.size == 0:
            continue
        
        closest = closest.sort_values(by = "distance")
        
        new_sizes = closest["size"] + size
        closest["new_size"] = new_sizes
        closest = closest.loc[closest["new_size"] < threshold]
        if closest.index.size == 0:
            continue
        
        closest = closest.iloc[0]
        
        replacements["from"].append(catchment_id)
        replacements["to"].append(closest.name)
        
        catchment_details_remaining.at[closest.name, "size"] += size
        catchment_details_remaining = catchment_details_remaining.drop(catchment_id)
        
    replacements = pd.DataFrame(replacements)
    return replacements


save_dir = pl.Path("./saves")
out_dir = pl.Path("./saves")
expansion_distance = 10.0 # in degrees

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
        
        combined_catchments_file = out_region_dir / f"combined_catchments_{region}.map"
        if combined_catchments_file.exists():
            print("> exists")
            #continue
        
        catchments_file = save_region_dir / f"country_catchments_{region}.map"
        ldd_file = save_region_dir / f"ldd_{region}.map"
        if not catchments_file.exists() or not ldd_file.exists():
            print("> skipping")
            continue
        
        pcr.setclone(str(catchments_file))
        catchments = pcr.readmap(str(catchments_file))
        ldd = pcr.readmap(str(ldd_file))
        
        catchment_details = calculate_catchment_details(catchments=catchments,
                                                        ldd=ldd)
        catchment_details = catchment_details.groupby(["catchment"]).mean()
        threshold = catchment_details["size"].max()
        
        replacements = calculate_combined_catchment_replacements(catchment_details=catchment_details,
                                                                 threshold=threshold,
                                                                 expansion_distance = 30.0)
        
        # Assign catchment id to all delta-connected catchments
        catchments_array = pcr.pcr2numpy(catchments, 0)
        combined_catchments_array = catchments_array.copy()
        
        for index, replacement in replacements.iterrows():
            from_catchment = replacement["from"]
            to_catchment = replacement["to"]
            from_catchment_mask = combined_catchments_array == from_catchment
            combined_catchments_array[from_catchment_mask] = to_catchment
        
        combined_catchments = pcr.numpy2pcr(pcr.Nominal, combined_catchments_array, 0)
        combined_catchments = pcr.ifthenelse(catchments == 0, pcr.nominal(0), combined_catchments)
        # pcr.aguila(catchments)
        # pcr.aguila(combined_catchments)
        # pcr.aguila(combined_catchments != catchments)
        
        combined_catchments_file.parent.mkdir(parents = True, exist_ok = True)
        pcr.report(combined_catchments, str(combined_catchments_file))
            
