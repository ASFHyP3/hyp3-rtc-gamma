#!/usr/bin/python
import logging
import os, sys
import argparse
import shutil
import zipfile
import math
import osr
import glob
import datetime
import numpy as np
from osgeo import gdal
from lxml import etree

import saa_func_lib as saa
from byteSigmaScale import byteSigmaScale
from createAmp import createAmp
from getSubSwath import get_bounding_box_file
from execute import execute
from getParameter import getParameter
from makeAsfBrowse import makeAsfBrowse
from ingest_S1_granule import ingest_S1_granule
from rtc2color import rtc2color
from asf_geometry import geometry_geo2proj
from getBursts import getBursts
from make_arc_thumb import pngtothumb
from geocode_sentinel import geocode_sentinel
from rtc_sentinel import rtc_sentinel_gamma

def create_dem_par_files(infiles):
    dpar_files = []
    for fi in infiles:
        basename = fi.replace(".SAFE","")
        dataType = "float"
        pixel_size = 30
        lat_max,lat_min,lon_max,lon_min = get_bounding_box_file(fi)
        post = 30
        create_dem_par(basename,dataType,pixel_size,lat_max,lat_min,lon_max,lon_min,post)
        dpar_files.append(basename+".par")

def create_dem_par(basename,dataType,pixel_size,lat_max,lat_min,lon_max,lon_min,post):
    demParIn = "{}_dem_par.in".format(basename)
    zone, false_north, y_min, y_max, x_min, x_max = geometry_geo2proj(lat_max,lat_min,lon_max,lon_min)
    logging.debug("Original Output Coordinates: {} {} {} {}".format(y_min, y_max, x_min, x_max))
    if post is not None:
        shift = 0
        x_max = math.ceil(x_max/post)*post+shift
        x_min = math.floor(x_min/post)*post-shift
        y_max = math.ceil(y_max/post)*post+shift
        y_min = math.floor(y_min/post)*post-shift
        logging.debug("Snapped Output Coordinates: {} {} {} {}".format(y_min, y_max, x_min, x_max))

    f = open(demParIn,"w")
    f.write("UTM\n")
    f.write("WGS84\n")
    f.write("1\n")
    f.write("{}\n".format(zone))
    f.write("{}\n".format(false_north))
    f.write("{}\n".format(basename))
    if "float" in dataType:
        f.write("REAL*4\n")
    elif "int16" in dataType:
        f.write("INTEGER*2\n")
    f.write("0.0\n")
    f.write("1.0\n")

    xsize = np.floor(abs((x_max-x_min)/pixel_size))
    ysize = np.floor(abs((y_max-y_min)/pixel_size))

    f.write("{}\n".format(int(xsize)))
    f.write("{}\n".format(int(ysize)))
    f.write("{} {}\n".format(-1.0*pixel_size,pixel_size))
    f.write("{} {}\n".format(y_max,x_min))
    f.close()

    return demParIn

def ingest_raw_stack(infiles,res):
  logging.info("Infile list: {}".format(infiles))
  look_fact = res / 10.0
  mgrd_list = []
  for x in xrange(len(infiles)):
    infile = infiles[x]
    if not os.path.exists(infile):
        logging.error("ERROR: Input file {} does not exist".format(infile))
        exit(1)
    if "zip" in infile:
        zip_ref = zipfile.ZipFile(infile, 'r')
        zip_ref.extractall(".")
        zip_ref.close()    
        infiles[x] = infile.replace(".zip",".SAFE")

    # Get list of files to process
    vvlist = glob.glob("{}/*/*vv*.tiff".format(infile))
    hhlist = glob.glob("{}/*/*hh*.tiff".format(infile))

    # Try to get the precision state vectors
#    try:
#        cmd = "get_orb.py {}".format(infile)
#        logging.info("Getting precision orbit information")
#        execute(cmd,uselogging=True)
#    except:
#        logging.warning("Unable to fetch precision state vectors... continuing")

    for fi in vvlist:
        outf= fi.split("/")[-1]
        outfile = outf.split(".")[0]
        outfile = outfile + ".mgrd"
        if os.path.isfile(outfile):
            logging.info("Output file {} already exists. Skipping.".format(outfile))
            mgrd_list.append(outfile)
        else:
            logging.info("Creating output file {}".format(outfile))
            ingest_S1_granule(infile,"VV",look_fact,outfile)
            mgrd_list.append(outfile)
         
    for fi in hhlist:
        outf= fi.split("/")[-1]
        outfile = outf.split(".")[0]
        outfile = outfile + ".mgrd"
        logging.info("Creating output file {}".format(outfile))
        ingest_S1_granule(infile,"HH",look_fact,outfile)
	mgrd_list.append(outfile)
 
  return(mgrd_list,infiles)
	
def create_diff_par_in():
    cfgdir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "config"))
    f = open("{}/diff_par.in".format(cfgdir),"r")
    outfile = "diff_par.in"
    g = open(outfile,"w")
    for line in f:
        g.write("{}".format(line))
    f.close()
    g.close()
    return(outfile)

def get_dates(file_list,sep="-"):
    date_and_file = []
    for fi in file_list:
        string = fi.split(sep)[4]
	date = string.split(sep)[0].upper()
	m = [fi,date]
	date_and_file.append(m)
    date_and_file.sort(key = lambda row: row[1])
    files = []
    dates = []
    for x in xrange(len(date_and_file)):
        files.append(date_and_file[x][0])
        dates.append(date_and_file[x][1])  
    return(files,dates)


def match_raw_stack(infiles,res):
    mgrd_files,infiles = ingest_raw_stack(infiles,res)
    dpar_files = create_dem_par_files(infiles)
    logging.info("mgrd file {}".format(mgrd_files))
    sorted_mgrd,dates = get_dates(mgrd_files)
    sorted_files,tmp = get_dates(infiles,sep="_")
    for x in xrange(len(dates)):
        if dates[x] != tmp[x]:
            logging.error("Mismatch between SAFE files and data files!")
            logging.error("dates {}".format(dates))
            logging.error("tmp {}".format(tmp))
            exit(-1)
    diff_par_in = create_diff_par_in()
    
    for x in xrange(len(sorted_mgrd)-1):
        diff_par = "{}_{}.par".format(dates[x],dates[x+1])
        if os.path.exists(diff_par):
            logging.info("Diff par file {} exists; skipping".format(diff_par))
        else:
            par1 = sorted_mgrd[x]+".par"
            par2 = sorted_mgrd[x+1]+".par"
            offs = sorted_mgrd[x].replace("mgrd","offs")
            coffs = sorted_mgrd[x].replace("mgrd","coffs")
            ccp = sorted_mgrd[x].replace("mgrd","ccp")
            cmd =  "create_diff_par {} {} {} 1 < {}".format(par1,par2,diff_par,diff_par_in)
            logging.info("Running command {}".format(cmd))
            os.system(cmd)
       
            try: 
                cmd = "init_offsetm {} {} {}".format(sorted_mgrd[x],sorted_mgrd[x+1],diff_par)
                logging.info("Running command {}".format(cmd))
                os.system(cmd)
            except:
                logging.info("Warning: could not determine initial offset for a single patch")
                logging.info("Warning: using multiple patches")
                cmd = "offset_pwrm {} {} {} {} {}".format(sorted_mgrd[x],sorted_mgrd[x+1],diff_par,offs,ccp) 
                logging.info("Running command {}".format(cmd))
                os.system(cmd)
                try:
                    cmd = "offset_fitm {} {} {} {}".format(offs,ccp,diff_par,coffs)
                except:
                    logging.error("ERROR: no patches had sufficient SNR for match for initial offset")
                    exit(-1)
    
            polyra = getParameter(diff_par,"range_offset_polynomial",uselogging=True)
            polyaz = getParameter(diff_par,"azimuth_offset_polynomial",uselogging=True)
       
            logging.info("Initial Range Polynomial: {}".format(polyra))
            logging.info("Initial Azimuth Polynomial: {}".format(polyaz))

            cmd = "offset_pwrm {} {} {} {} {}".format(sorted_mgrd[x],sorted_mgrd[x+1],diff_par,offs,ccp) 
            logging.info("Running command {}".format(cmd))
            os.system(cmd)
            try:
                cmd = "offset_fitm {} {} {} {}".format(offs,ccp,diff_par,coffs)
            except:
                logging.error("ERROR: failed to coregister scene!")
                exit(-1)

    prod_dir = os.path.join(os.getcwd(),"RTC_PRODUCTS")
    if not os.path.exists(prod_dir):
        os.mkdir(prod_dir)
    rtc_granule(sorted_files[0],prod_dir,res)
    polyaz = 0 
    polyra = 0
    for x in xrange(1,len(sorted_mgrd)):
        diff_par = "{}_{}.par".format(dates[x-1],dates[x])
        diff_par = os.path.join(os.getcwd(),diff_par)
        polyaz,polyra = get_poly(polyaz,polyra,diff_par) 
        new_par = fix_diff_par(polyaz,polyra,diff_par)
        rtc_granule(sorted_files[x],prod_dir,res,stack=new_par)


def rtc_granule(fi,prod_dir,res,stack=None):
    back = os.getcwd()
    mydir = os.path.basename(fi)
    mydir = mydir.replace(".SAFE","")
    mydir = mydir[17:32]
    if os.path.exists(mydir):
        logging.info("Old {} directory found; deleting".format(mydir))
        shutil.rmtree(mydir)
    os.mkdir(mydir)
    os.chdir(mydir)
    os.symlink("../{}".format(fi),fi)

    if res == 10.0:
        look_fact = 1
    else: 
        look_fact = int(res/10.0) * 2

    rtc_sentinel_gamma(fi,matchFlag=True,deadFlag=True,gammaFlag=True,res=res,pwrFlag=True,looks=look_fact,
                       terms=1,stack=stack,noCrossPol=True)

    for tmp in glob.glob("PRODUCT/*_V*.tif"):
        shutil.copy(tmp,prod_dir)
    for tmp in glob.glob("PRODUCT/*_H*.tif"):
        shutil.copy(tmp,prod_dir)

    os.chdir(back)
   
def fix_diff_par(az,ra,diff_par):
    f = open(diff_par,"r")
    outfile = diff_par.replace(".par","") + "_accum.par"
    g = open(outfile,"w")
    for line in f:
        if "range_offset_polynomial" in line:
            g.write("range_offset_polynomial: {} 0.0 0.0 0.0 0.0 0.0\n".format(ra))
        elif "azimuth_offset_polynomial" in line:
            g.write("azimuth_offset_polynomial: {} 0.0 0.0 0.0 0.0 0.0\n".format(az))
        else:
            g.write("{}".format(line))
    f.close()
    g.close()

def get_poly(az,ra,diff_par):
     polyra = getParameter(diff_par,"range_offset_polynomial",uselogging=True)
     polyaz = getParameter(diff_par,"azimuth_offset_polynomial",uselogging=True)
     ra_new = float(polyra.split()[0])
     az_new = float(polyaz.split()[0])
     return az_new+az,ra_new+ra


if __name__ == '__main__':

  parser = argparse.ArgumentParser(prog='ingest_stack',description="ingest a stack of S1 GRD images")
  parser.add_argument("infile",nargs="+",help="zip files or SAFE directories")
  parser.add_argument("-r","--res",type=float,default=30,help="Desired output pixel spacing for stack (default 30m)")
  args = parser.parse_args()

  logFile = "{}_{}_log.txt".format("ingest_stack",os.getpid())
  logging.basicConfig(filename=logFile,format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',level=logging.DEBUG)
  logging.getLogger().addHandler(logging.StreamHandler())
  logging.info("Starting run")

  match_raw_stack(args.infile,args.res)
