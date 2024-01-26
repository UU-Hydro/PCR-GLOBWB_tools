import pathlib as pl
import pcraster as pcr

input_dir = pl.Path("./input")
out_dir = pl.Path("./saves")

resolutions = [dir.stem for dir in input_dir.iterdir() if dir.is_dir() and "sec" in dir.stem or "min" in dir.stem]

resolution = resolutions[0]
for resolution in resolutions:
    print(f"processing resolution {resolution}")
    
    input_resolution_dir = input_dir / resolution
    out_resolution_dir = out_dir / resolution
    
    catchments_file = out_resolution_dir / "catchments_global.map"
    if catchments_file.exists():
        print("> exists")
        continue
    
    ldd_file = input_resolution_dir / "ldd.map"

    ldd = pcr.readmap(str(ldd_file))
    pits = pcr.pit(ldd)
    catchments = pcr.catchment(ldd, pits)

    catchments_file.parent.mkdir(parents = True, exist_ok = True)
    pcr.report(catchments, str(catchments_file))