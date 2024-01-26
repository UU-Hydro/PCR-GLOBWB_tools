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


def aquire_bounding_box_attributes(pcrmap: pcr.Field) -> dict[str, float]:
    attributes = aquire_map_attributes(calc_lonlat = True)

    array = pcr.pcr2numpy(pcrmap, False).astype(np.bool_)
    lat_sel = np.max(array, axis = 1)
    lon_sel = np.max(array, axis = 0)

    lons = attributes["lons"][lon_sel]
    lats = attributes["lats"][lat_sel]

    attributes["origional"] = attributes.copy()
    attributes["lons"] = lons
    attributes["lats"] = lats
    attributes["west"] = min(lons) - attributes["resolution"] / 2
    attributes["north"] = max(lats) + attributes["resolution"] / 2
    attributes["ncols"] = len(lons)
    attributes["nrows"] = len(lats)
    attributes["lon_sel"] = lon_sel
    attributes["lat_sel"] = lat_sel
    
    return attributes


def report_bounding_box(pcrmap: pcr.Field,
                        output_file: pl.Path,
                        input_pcrtype = pcr.Nominal,
                        output_pcrtype = pcr.Boolean) -> None:
    
    if output_pcrtype == pcr.Boolean:
        pcrmap = pcrmap != 0
    pcrmap_array = pcr.pcr2numpy(pcr.scalar(pcrmap), -1)
    
    if np.max(pcrmap_array) <= 0:
        raise ValueError("No landmask found")
    
    active_pcrmap_mask = pcr.ifthenelse(pcrmap != 0, pcr.boolean(1), pcr.boolean(0))
    attributes = aquire_bounding_box_attributes(pcrmap=active_pcrmap_mask)
    
    active_indices = np.ix_(attributes["lat_sel"], attributes["lon_sel"])
    active_pcrmap_array = pcrmap_array[active_indices]
    
    pcr.setclone(attributes["nrows"],
                 attributes["ncols"],
                 attributes["resolution"],
                 attributes["west"],
                 attributes["north"])
    
    active_pcrmap = pcr.numpy2pcr(output_pcrtype, active_pcrmap_array, -1)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    pcr.report(active_pcrmap, str(output_file))
    
    pcr.setclone(attributes["origional"]["nrows"],
                 attributes["origional"]["ncols"],
                 attributes["origional"]["resolution"],
                 attributes["origional"]["west"],
                 attributes["origional"]["north"])


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
        
        catchments_file = save_region_dir / f"combined_catchments_{region}.map"
        if not catchments_file.exists():
            print("> skipping")
            continue
        
        pcr.setclone(str(catchments_file))
        catchments = pcr.readmap(str(catchments_file))
        
        catchments_array = pcr.pcr2numpy(catchments, 0)
        catchment_ids = np.unique(catchments_array)
        catchment_ids = catchment_ids[catchment_ids != 0]
        
        domains_array = np.zeros_like(catchments_array)
        for index, catchment_id in enumerate(catchment_ids):
            domain_index = index + 1
            catchment_sel = catchments_array == catchment_id
            domains_array[catchment_sel] = domain_index
        
        domains = pcr.numpy2pcr(pcr.Nominal, domains_array, 0)
        # pcr.aguila(domains)
        
        domains_file = out_region_dir / f"domains_{region}.map"
        domains_file.parent.mkdir(parents = True, exist_ok = True)
        pcr.report(domains, str(domains_file))
        
        if catchment_ids.size == 1:
            catchment_id = catchment_ids[0]
            landsurface = catchments == int(catchment_id)
            # pcr.aguila(landsurface)
            
            landsurface_file = out_region_dir / f"landmask_{region}.map"
            landsurface_file.parent.mkdir(parents = True, exist_ok = True)
            pcr.report(landsurface, str(landsurface_file))
        
        else:
            for index, catchment_id in enumerate(catchment_ids):
                domain_index = index + 1
                out_domain_dir = out_region_dir / f"{domain_index}"
                
                landsurface = catchments == int(catchment_id)
                # pcr.aguila(landsurface)
                
                landsurface_file = out_domain_dir / f"landmask_{domain_index}.map"
                landsurface_file.parent.mkdir(parents = True, exist_ok = True)
                report_bounding_box(pcrmap=landsurface,
                                    output_file=landsurface_file,
                                    input_pcrtype=pcr.Boolean,
                                    output_pcrtype=pcr.Boolean)
                
