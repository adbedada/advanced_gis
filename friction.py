'''

This script contains functions to process geospatial data and compute cost distance'
cli:
    python friction.py single_friction_cost -i "input_file.tiff" -o "output_friction.tiff"
'''


import rasterio
import gdal, ogr
import numpy as np
import gdalconst
import rasterio as rio
import rasterio.merge, rasterio.mask
from skimage import graph
import argparse


pixelWidth = 250
pixelHeight = -250
_rows = 120
_cols = 148

def resample(inputfile, outputfile):

    src = gdal.Open(inputfile, gdal.GA_ReadOnly)
    src_trfm = src.GetGeoTransform()
    src_proj = src.GetProjection()
    cols = _cols
    rows = _rows
    originX = src_trfm[0]
    originY = src_trfm[3]

    dst = gdal.GetDriverByName('GTiff').Create(outputfile, cols, rows, 1, gdal.GDT_Float32)
    dst.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    dst.SetProjection(src_proj)
    ref_proj = dst.GetProjection()
    gdal.ReprojectImage(src, dst, src_proj, ref_proj, gdalconst.GRA_Mode)

    del dst #flush


def single_friction_cost(rasterfile, frictionfile):
    img = gdal.Open(rasterfile, gdal.GA_ReadOnly)
    arr = img.ReadAsArray()
    value = np.where(arr == 3, 3.45, 15)

    src = gdal.Open(inputfile, gdal.GA_ReadOnly)
    src_trfm = src.GetGeoTransform()
    src_proj = src.GetProjection()
    cols = _cols
    rows = _rows
    originX = src_trfm[0]
    originY = src_trfm[3]

    dst = gdal.GetDriverByName('GTiff').Create(frictionfile, cols, rows, 1, gdal.GDT_Float32)
    dst.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    dst.SetProjection(src_proj)

    outband = dst.GetRasterBand(1)
    outband.WriteArray(value)
    outband.FlushCache()


def rasterize(shapefile, exampletiff, outputfile, options=["ATTRIBUTE=val"]):
    src = gdal.Open(exampletiff, gdal.GA_ReadOnly)
    src_trfm = src.GetGeoTransform()
    src_proj = src.GetProjection()
    cols = _cols
    rows = _rows
    originX = src_trfm[0]
    originY = src_trfm[3]

    shp = ogr.Open(shapefile)
    layer = shp.GetLayer()

    dst = gdal.GetDriverByName('GTiff').Create(outputfile, cols, rows, 1, gdal.GDT_Float32)
    dst.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    dst.SetProjection(src_proj)
    band = dst.GetRasterBand(1)
    NoData_value = -9999
    band.SetNoDataValue(NoData_value)
    band.FlushCache()
    gdal.RasterizeLayer(dst, [1], layer, options=options)


def merge_rasters(inputfile1,inputfile2, outputfile):
    src1 = rio.open(inputfile1)
    src2 = rio.open(inputfile2)
    transform = src1.transform
    crs = src1.crs
    list_files = (src1, src2)
    merged_files, road_affine = rasterio.merge.merge(list_files)
    dtype = merged_files.dtype

    with rasterio.open(outputfile,
                       'w', driver='GTiff',
                       height=src1.shape[0],
                       width=src1.shape[1],
                       count=1,
                       dtype=dtype,
                       crs=crs,
                       transform=transform) as dst:
        dst.write(merged_files[0], 1)


def multiple_friction_costs(inputfile, outputfile):
    img = gdal.Open(inputfile, gdal.GA_ReadOnly)
    arr = img.ReadAsArray()
    val1 = np.where(arr == 65, 0.23, arr)
    val2 = np.where(val1 == 55, 0.26, val1)
    val3 = np.where(val2 == 35, 0.43, val2)
    val4 = np.where(val3 == 25, 0.6, val3)

    src = gdal.Open(inputfile, gdal.GA_ReadOnly)
    src_trfm = src.GetGeoTransform()
    src_proj = src.GetProjection()
    cols = _cols
    rows = _rows
    originX = src_trfm[0]
    originY = src_trfm[3]

    dst = gdal.GetDriverByName('GTiff').Create(outputfile, cols, rows, 1, gdal.GDT_Float32)
    dst.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    dst.SetProjection(src_proj)

    outband = dst.GetRasterBand(1)
    outband.WriteArray(val4)
    outband.FlushCache()


def replace_null(inputfile1, inputfile2, outputfile):
    src1 = rio.open(inputfile1)
    arr1 = src1.read()
    src2 = rio.open(inputfile2)
    arr2 = src2.read()

    val = arr1 == -9999
    arr1[val] = arr2[val]
    dtype = arr1.dtype
    transform = src1.transform
    crs = src1.crs
    with rasterio.open(outputfile,
                       'w', driver='GTiff',
                       height=src1.shape[0],
                       width=src1.shape[1],
                       count=1,
                       dtype=dtype,
                       crs=crs,
                       transform=transform) as dst:
        dst.write(arr1[0], 1)


def compute_cost(inputfile, outputfile):
    src1 = rio.open(inputfile)
    arr1 = src1.read()[0]
    lg = graph.MCP_Geometric(arr1)
    lcd = lg.find_costs(starts =[(41,3)])[0]
    dtype = lcd.dtype
    transform = src1.transform
    crs = src1.crs
    with rasterio.open(outputfile,
                       'w', driver='GTiff',
                       height=src1.shape[0],
                       width=src1.shape[1],
                       count=1,
                       dtype=dtype,
                       crs=crs,
                       transform=transform) as dst:
        dst.write(lcd, 1)

############ PRODUCE COST FUNCTIONS ##################

#resample("landuse/hdr.adf","u-resampled.tiff")
#single_friction_cost("lu-resampled.tiff", "output/lu-friction.tiff")
#rasterize("PavedRoads.shp","output/lu-friction.tiff", "output/PavedRoads.tiff", options=["ATTRIBUTE=SpeedLimit"])
#merge_rasters("output/PavedRoads.tiff", "output/UnpavedRoads.tiff", "output/merged_roads.tiff")
#multiple_friction_cost("output/merged_roads.tiff", "output/roads_friction.tiff")
#replace_null("output/roads_friction.tiff", "output/lu-friction.tiff","output/all_friction.tiff")
#compute_cost("output/all_friction.tiff", "output/friction_cost.tiff")


parser = argparse.ArgumentParser(description='Compute Cost Distance.')
subparsers = parser.add_subparsers()

parser_resample = subparsers.add_parser('resample')
parser_resample.add_argument('-i', '--input', dest='inputfile', help="Raster to resample")
parser_resample.add_argument('-o', '--output', dest='outputfile', help="Name of resampled output")

parser_one_friction = subparsers.add_parser('single_friction_cost')
parser_one_friction.add_argument('-i', '--input', dest='rasterfile', help="N=Input raster to compute friction cost")
parser_one_friction.add_argument('-o', '--output', dest='frictionfile', help="Name for computed friction cost file")

parser_rasterize = subparsers.add_parser('rasterize')
parser_rasterize.add_argument('-s', '--shapefile', dest='shapefile', help="Input raster to compute friction cost")
parser_rasterize.add_argument('-e', '--example', dest='exampletiff', help="Name for computed friction cost file")
parser_rasterize.add_argument('-o', '--output', dest='output', help="Name for computed friction cost file")
parser_rasterize.set_defaults(func=rasterize)

parser_merge_raster = subparsers.add_parser('merge_rasters')
parser_merge_raster.add_argument('-i1', '--input1', dest='inputfile1', help="Nnput raster to compute friction cost")
parser_merge_raster.add_argument('-i2', '--input2', dest='inputfile2', help="Nnput raster to compute friction cost")
parser_merge_raster.add_argument('-o', '--output', dest='outputfile', help="Name for computed friction cost file")

parser_replace_null = subparsers.add_parser('replace_null')
parser_replace_null.add_argument('-i1', '--input1', dest='inputfile1', help="Input raster to compute friction cost")
parser_replace_null.add_argument('-i2', '--input2', dest='inputfile2', help="Input raster to compute friction cost")
parser_replace_null.add_argument('-o', '--output', dest='outputfile', help="Name for computed friction cost file")

parser_compute_cost = subparsers.add_parser('compute_cost')
parser_compute_cost.add_argument('-i1', '--input1', dest='inputfile1', help="Input raster to compute friction cost")
parser_compute_cost.add_argument('-i2', '--input2', dest='inputfile2', help="Input raster to compute friction cost")
parser_compute_cost.add_argument('-o', '--output', dest='outputfile', help="Name for computed friction cost file")

parser.parse_args(['--resample', 'single_friction_costew', '--baz', 'Z'])
args = parser.parse_args()
inputfile = str(args.inputfile)
outputfile = str(args.outputfile)
#
if __name__ == "__main__":
    resample(inputfile, outputfile)




