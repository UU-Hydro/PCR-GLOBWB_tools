from typing import Optional
import pathlib as pl
import numpy as np
import numba as nb
import pcraster as pcr
import geopandas as gpd
from rasterio import features
from rasterio import transform
from rasterio.enums import MergeAlg
import copy


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


@nb.njit
def expand_fill(array: np.ndarray,
                fill_value: bool,
                distance: int,) -> np.ndarray:
    
    array_indices = np.where(array)
    for index, (y_index, x_index) in enumerate(zip(*array_indices)):
        if index % 1000 == 0:
            print(f"processing index {index} of {array_indices[0].size}")
        
        y_index_min = max(y_index - distance, 0)
        y_index_max = min(y_index + distance + 1, array.shape[0])
        x_index_min = max(x_index - distance, 0)
        x_index_max = min(x_index + distance + 1, array.shape[1])
        array[y_index_min:y_index_max, x_index_min:x_index_max] = fill_value
    return array


def find_nearest_neighbour(array: np.ndarray,
                           fill_value,
                           distance: int,
                           mask: Optional[np.ndarray] = None,
                           simple_square: bool = False) -> np.ndarray:
    
    if mask is None:
        mask = array == fill_value
        
    if not simple_square:
        window_shape = distance * 2
        if window_shape % 2 == 0:
            window_shape += 1
            
        shape_x_indices = np.arange(window_shape) - (window_shape - 1) / 2
        shape_y_indices = np.arange(window_shape) - (window_shape - 1) / 2
        shape_x_indices = np.tile(shape_x_indices, (window_shape, 1))
        shape_y_indices = np.tile(shape_y_indices, (window_shape, 1)).T
        window_distances = np.sqrt(shape_x_indices ** 2 + shape_y_indices ** 2)
        window_sel = window_distances <= distance + 1e-6
    
    mask_indices = np.where(mask)
    for index, (y_index, x_index) in enumerate(zip(*mask_indices)):
        if index % 1000 == 0:
            print(f"processing index {index} of {mask_indices[0].size}")
        
        y_index_min = y_index - distance
        y_index_max = y_index + distance + 1
        x_index_min = x_index - distance
        x_index_max = x_index + distance + 1
        
        if not simple_square:
            window_sel_used = window_sel
            if y_index_min < 0 or y_index_max > array.shape[0] or x_index_min < 0 or x_index_max > array.shape[1]:    
                window_sel_used = copy.deepcopy(window_sel_used)        
                if y_index_min < 0:
                    window_sel_used = window_sel_used[-y_index_min:, :]
                if y_index_max > array.shape[0]:
                    window_sel_used = window_sel_used[:-(y_index_max - array.shape[0]), :]
                if x_index_min < 0:
                    window_sel_used = window_sel_used[:, -x_index_min:]
                if x_index_max > array.shape[1]:
                    window_sel_used = window_sel_used[:, :-(x_index_max - array.shape[1])]
        
        y_index_min = max(y_index_min, 0)
        y_index_max = min(y_index_max, array.shape[0])
        x_index_min = max(x_index_min, 0)
        x_index_max = min(x_index_max, array.shape[1])
        
        neighbours = array[y_index_min:y_index_max, x_index_min:x_index_max]
        if not simple_square:
            neighbours = neighbours[window_sel_used]
            
        neighbours = neighbours[neighbours != fill_value]
        if neighbours.size == 0:
            continue
        
        neighbour_unique, neighbour_count = np.unique(neighbours, return_counts = True)
        neighbour = neighbour_unique[np.argmax(neighbour_count)]
        array[y_index, x_index] = neighbour
    
    return array


input_dir = pl.Path("./input")
save_dir = pl.Path("./saves")
out_dir = pl.Path("./saves")
expansion_resolution = 0.25 # in degrees

deltas_file = input_dir / "global_delta_map" / "global_map.shp"
deltas_gdf = gpd.read_file(deltas_file)
deltas_gdf = deltas_gdf.loc[deltas_gdf["DeltaID"] > 0]
deltas_geoms = [(geom, int(id)) for geom, id in zip(deltas_gdf["geometry"], deltas_gdf["DeltaID"])]

resolutions = [dir.stem for dir in input_dir.iterdir() if dir.is_dir() if "sec" in dir.stem or "min" in dir.stem]
resolution = resolutions[-1]
for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    input_resolution_dir = input_dir / resolution
    save_resolution_dir = save_dir / resolution
    out_resolution_dir = out_dir / resolution
    
    out_deltas_file = out_resolution_dir / "deltas_global.map"
    if out_deltas_file.exists():
        print("> exists")
        #continue
        
    ldd_file = input_resolution_dir / "ldd.map"
    if not ldd_file.exists():
        print("> skipping")
        continue
    
    ldd = pcr.readmap(str(ldd_file))
    attributes = aquire_map_attributes()
    
    resolution_transform = transform.from_origin(west = attributes["west"],
                                                north = attributes["north"],
                                                xsize = attributes["resolution"],
                                                ysize = attributes["resolution"])
    
    deltas_array = features.rasterize(deltas_geoms,
                                    out_shape = (attributes["nrows"], attributes["ncols"]),
                                    transform = resolution_transform,
                                    fill = 0,
                                    all_touched = True,
                                    merge_alg = MergeAlg.replace,
                                    dtype = np.int32)
    
    deltas = pcr.numpy2pcr(pcr.Nominal, deltas_array, 0)
    # pcr.aguila(deltas)
    
    # Iterate a couple of times to get the downstream areas of deltas as well
    expanded_deltas_array = pcr.pcr2numpy(deltas, 0)
    ldd_array = pcr.pcr2numpy(ldd, 0)
    
    expansion_distance = int(expansion_resolution / attributes["resolution"])
    
    missing_array = expanded_deltas_array != 0
    missing_array = expand_fill(array = missing_array,
                                fill_value=True,
                                distance = expansion_distance)
    missing_array = np.logical_and(missing_array, ldd_array != 0)
    missing_array = np.logical_and(missing_array, expanded_deltas_array == 0)
    # missing = pcr.numpy2pcr(pcr.Boolean, missing_array, 0)
    # pcr.aguila(missing)
    
    expanded_deltas_array = find_nearest_neighbour(mask = missing_array,
                                                   array = expanded_deltas_array,
                                                   fill_value=0,
                                                   distance = expansion_distance)
    expanded_deltas = pcr.numpy2pcr(pcr.Nominal, expanded_deltas_array, 0)
    # pcr.aguila(expanded_deltas)
    
    out_deltas_file.parent.mkdir(parents = True, exist_ok = True)
    pcr.report(expanded_deltas, str(out_deltas_file))