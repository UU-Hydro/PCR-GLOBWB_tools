from typing import Optional 
import pathlib as pl
import numpy as np
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


def derive_catchment_from_point(catchments: pcr.Field,
                                lon: float,
                                lat: float) -> pcr.Field:
    
    attributes = aquire_map_attributes(calc_lonlat = True)
        
    lon_diff = np.abs(attributes["lons"] - lon)
    lat_diff = np.abs(attributes["lats"] - lat)
    if np.min(lon_diff) > 0.5 * attributes["resolution"]:
        raise ValueError("Point longitude not found in map")
    if np.min(lat_diff) > 0.5 * attributes["resolution"]:
        raise ValueError("Point latitude not found in map")    

    lon_index = np.argmin(lon_diff)
    lat_index = np.argmin(lat_diff)

    catchments_array = pcr.pcr2numpy(catchments, 0)    
    catchment_id = int(catchments_array[lat_index, lon_index])
    
    catchments = pcr.ifthenelse(catchments == catchment_id, catchments, pcr.nominal(0))
    return catchments


def derive_catchments_from_extent(catchments: pcr.Field,
                                  ldd: Optional[pcr.Field],
                                  lons: tuple[float, float],
                                  lats: tuple[float, float],) -> list[int]:
        
    lon_min = min(lons)
    lat_min = min(lats)
    lon_max = max(lons)
    lat_max = max(lats)

    attributes = aquire_map_attributes(calc_lonlat = True)

    lons_sel = np.logical_and(attributes["lons"] >= lon_min,
                              attributes["lons"] <= lon_max)
    lats_sel = np.logical_and(attributes["lats"] >= lat_min,
                              attributes["lats"] <= lat_max)
    
    catchments_array = pcr.pcr2numpy(catchments, 0)
    
    pit = pcr.pit(ldd)
    pit_array = pcr.pcr2numpy(pit, 0)
    
    extent_indices = np.ix_(lats_sel, lons_sel)
    extent_catchments_array = catchments_array[extent_indices]
    extent_pit_array = pit_array[extent_indices]
    catchment_ids = extent_catchments_array[extent_pit_array > 0]
    
    catchment_mask_array = (np.isin(catchments_array, catchment_ids)).astype(np.int8)
    catchment_mask_array[catchments_array == 0] = -1
    catchment_mask = pcr.numpy2pcr(pcr.Boolean, catchment_mask_array, -1)
    
    catchments = pcr.ifthenelse(catchment_mask, catchments, pcr.nominal(0))
    return catchments


def report_bounding_box(pcrmap: pcr.Field,
                        output_file: pl.Path,
                        input_pcrtype = pcr.Nominal,
                        output_pcrtype = pcr.Boolean) -> None:
    
    if output_pcrtype == pcr.Boolean:
        pcrmap = pcrmap != 0
    pcrmap_array = pcr.pcr2numpy(pcr.scalar(pcrmap), -1)
    
    if np.max(pcrmap_array) <= 0:
        raise ValueError("No landmask found")
    
    active_pcrmap_mask = pcr.defined(pcrmap)
    if input_pcrtype == pcr.Nominal:
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
    

input_dir = pl.Path("./input")
save_dir = pl.Path("./saves")
out_dir = pl.Path("./saves")

points = {"Rhine": {"lon": 8.465353254198007, 
                    "lat": 49.2819248583303},
          "Po": {"lon": 10.667410622484812,
                 "lat": 45.039166232422225},
          "Tugela": {"lon": 30.169456045054208,
                     "lat": -28.750632363733498}}
extents = {"Europe": {"lons": (-11, 42),
                      "lats": (33, 73)}}

resolutions = [dir.stem for dir in save_dir.iterdir() if dir.is_dir()]

for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    input_resolution_dir = input_dir / resolution
    save_resolution_dir = save_dir / resolution
    out_resolution_dir = out_dir / resolution

    catchments_file = save_resolution_dir / "delta_catchments_global.map"
    pcr.setclone(str(catchments_file))
    catchments = pcr.readmap(str(catchments_file))

    ldd_file = input_resolution_dir / "ldd.map"
    ldd = pcr.readmap(str(ldd_file))

    region, point = list(points.items())[0]
    for region, point in points.items():
        print(f"processing region {region}")
        
        out_region_dir = out_resolution_dir / region
        
        lon = point["lon"]
        lat = point["lat"]
        
        catchments_file = out_region_dir / f"delta_catchments_{region}.map"
        ldd_file = out_region_dir / f"ldd_{region}.map"
        if catchments_file.exists() and ldd_file.exists():
            print("> exists")
            continue
        
        point_catchments = derive_catchment_from_point(catchments = catchments,
                                                    lon = lon,
                                                    lat = lat)
        # pcr.aguila(point_catchments)
        
        point_ldd = pcr.ifthen(point_catchments != 0, ldd)
        # pcr.aguila(point_ldd)
        
        catchments_file.parent.mkdir(parents = True, exist_ok = True)
        report_bounding_box(pcrmap = point_catchments,
                            output_file = catchments_file,
                            output_pcrtype = pcr.Nominal)
        
        ldd_file.parent.mkdir(parents = True, exist_ok = True)
        report_bounding_box(pcrmap = point_ldd,
                            output_file = ldd_file,
                            input_pcrtype=pcr.Ldd,
                            output_pcrtype = pcr.Ldd)

    region, extent = list(extents.items())[0]
    for region, extent in extents.items():
        print(f"processing region {region}")
        
        out_region_dir = out_resolution_dir / region
        
        lons = extent["lons"]
        lats = extent["lats"]
        
        catchments_file = out_region_dir / f"delta_catchments_{region}.map"
        ldd_file = out_region_dir / f"ldd_{region}.map"
        if catchments_file.exists() and ldd_file.exists():
            print("> exists")
            continue
        
        extent_catchments = derive_catchments_from_extent(catchments = catchments,
                                                        ldd=ldd,
                                                        lons = lons,
                                                        lats = lats)
        # pcr.aguila(extent_catchments)
        
        extent_ldd = pcr.ifthen(extent_catchments != 0, ldd)
        # pcr.aguila(extent_ldd)
        
        catchments_file.parent.mkdir(parents = True, exist_ok = True)
        report_bounding_box(pcrmap = extent_catchments,
                            output_file = catchments_file,
                            output_pcrtype = pcr.Nominal)
        
        ldd_file.parent.mkdir(parents = True, exist_ok = True)
        report_bounding_box(pcrmap = extent_ldd,
                            output_file = ldd_file,
                            input_pcrtype=pcr.Ldd,
                            output_pcrtype = pcr.Ldd)