from typing import Optional
import pathlib as pl
import numpy as np
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
expansion_distance = 1.0 # in degrees

countries_file = input_dir / "global_country_map" / "geoBoundariesCGAZ_ADM0.shp"
countries_gdf = gpd.read_file(countries_file)
countries_gdf = countries_gdf.loc[countries_gdf["shapeType"] == "ADM0"]
countries_gdf.loc[:, "CountryID"] = countries_gdf.groupby(['shapeName']).ngroup() + 1
countries_geoms = [(geom, int(id)) for geom, id in zip(countries_gdf["geometry"], countries_gdf["CountryID"])]

resolutions = [dir.stem for dir in save_dir.iterdir() if dir.is_dir()]
# resolutions = ["03min"]

resolution = resolutions[0]
for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    save_resolution_dir = save_dir / resolution
    out_resolution_dir = out_dir / resolution

    regions = [dir.stem for dir in save_resolution_dir.iterdir() if dir.is_dir()]
    # regions = ["Europe"]
    
    region = regions[0]
    for region in regions:
        print(f"processing region {region}")
        
        save_region_dir = save_resolution_dir / region
        out_region_dir = out_resolution_dir / region
        
        out_countries_file = out_region_dir / f"countries_{region}.map"
        if out_countries_file.exists():
            print("> exists")
            continue
        
        ldd_file = save_region_dir / f"ldd_{region}.map"
        if not ldd_file.exists():
            print("> skipping")
            continue
        
        pcr.setclone(str(ldd_file))
        ldd = pcr.readmap(str(ldd_file))
        attributes = aquire_map_attributes()
        
        resolution_transform = transform.from_origin(west = attributes["west"],
                                                    north = attributes["north"],
                                                    xsize = attributes["resolution"],
                                                    ysize = attributes["resolution"])
        
        expanded_countries_array = features.rasterize(countries_geoms,
                                            out_shape = (attributes["nrows"], attributes["ncols"]),
                                            transform = resolution_transform,
                                            fill = 0,
                                            all_touched = True,
                                            merge_alg = MergeAlg.replace,
                                            dtype = np.int32)
        
        countries = pcr.numpy2pcr(pcr.Nominal, expanded_countries_array, 0)
        countries = pcr.clump(countries) # Remove small islands, we combine these back to the nearest country later
        # pcr.aguila(clump_countries)
        
        expanded_countries_array = pcr.pcr2numpy(countries, 0)
        ldd_array = pcr.pcr2numpy(ldd, 0)
        
        # Iterate a couple of times to get the nearby areas of countries as well
        initial_window = 2
        niters = int(expansion_distance / attributes["resolution"] / initial_window)
        for expand_iter in range(niters):
            print(f"expanding countries iter {expand_iter + 1}/{niters}")
            missing_array = np.logical_and(expanded_countries_array == 0,
                                           ldd_array != 0)
            if missing_array.sum() == 0:
                break
            expanded_countries_array = find_nearest_neighbour(mask = missing_array,
                                                              array = expanded_countries_array,
                                                              fill_value=0,
                                                              distance = initial_window**(expand_iter + 1))
        expanded_countries = pcr.numpy2pcr(pcr.Nominal, expanded_countries_array, 0)
        # pcr.aguila(expanded_countries)
        
        missing_array = np.logical_and(expanded_countries_array == 0,
                                ldd_array != 0)
        if missing_array.sum() != 0:
            raise ValueError("Warning: not all countries are defined")
            
        out_countries_file.parent.mkdir(parents = True, exist_ok = True)
        pcr.report(expanded_countries, str(out_countries_file))