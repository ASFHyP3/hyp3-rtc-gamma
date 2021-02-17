import numpy as np
import rasterio as rio
import xarray as xr
import rioxarray
import logging
import os
import argparse
import re
import pycrs
from osgeo import gdal
from datetime import datetime
from rasterio.enums import Resampling


#
#     The code assumes auxiliary files are in the same location as the tiff files
# 

#
# example: S1A_IW_RT30_20180727T161143_G_gpn_VV
#


def parse_asf_rtc_name(infile):
    data = infile.split('_')
    parsed = {}
    try:
        parsed['platform'] = data[0]
        parsed['beam_mode'] = data[1]
        parsed['pixel_spacing'] = data[2][2:-1]
        parsed['start_time'] = data[3]
    except:
        raise Exception(f'ERROR: Unable to parse filename {infile}')

    if 'G' in data[4]:
        parsed['package'] = 'Gamma'
    elif 'S' in data[4]:
        parsed['package'] = 'S1TBX'

    try:
        if data[5][0] == 'g':
            parsed['radiometry'] = 'gamma0'
        elif data[5][0] == 's':
            parsed['radiometry'] = 'sigma0'
        elif data[5][0] == 'b':
            parsed['radiometry'] = 'beta0'
    except:
        logging.error(f'ERROR: Unable to determine radiometry from string {infile} letter {data[5][0]}')
        raise Exception(f'ERROR: Unable to determine radiometry from string {infile} letter {data[5][0]}')

    try:
        if data[5][1] == 'p':
            parsed['scale'] = 'power'
        elif data[5][1] == 'a':
            parsed['scale'] = 'amplitude'
    except:
        logging.error(f'ERROR: Unable to determine scaling from string {infile} letter {data[5][1]}')
        raise Exception(f'ERROR: Unable to determine scaling from string {infile} letter {data[5][1]}')

    try:
        if data[5][2] == 'n':
            parsed['filtered'] = 'no'
        elif data[5][2] == 'f':
            parsed['filtered'] = 'yes'
    except:
        logging.error(f'ERROR: Unable to determine filtering from string {infile} letter {data[5][2]}')
        raise Exception(f'ERROR: Unable to determine filtering from string {infile} letter {data[5][2]}')

    parsed['polarization'] = data[6]

    return parsed


def get_dataset(infile):
    logging.debug('Input file: {}'.format(infile))
    dataset = rio.open(infile)
    image_dts = os.path.basename(infile)[12:27]
    granule = get_granule(infile, image_dts)
    proc_dt = get_processing_datetime(infile, granule)
    ulx, uly = dataset.transform * (0, 0)
    lrx, lry = dataset.transform * (dataset.width, dataset.height)
    x_extent = [ulx, lrx]
    y_extent = [uly, lry]
    return dataset, image_dts, proc_dt, x_extent, y_extent, granule


def get_granule(infile, image_dts):
    lines, logfile = read_log_file(infile)
    search_string = 'S1[AB]_.._...._1S[DS][HV]_' + '{}'.format(
        image_dts) + '_\d{8}T\d{6}_\d+_([0-9A-Fa-f]+)_([0-9A-Fa-f]+)'
    granule = re.search(search_string, lines)
    if not granule:
        raise Exception('ERROR: No granule name found in {}'.format(logfile))
    granule = granule.group(0)
    logging.info('Found granule {}'.format(granule))
    return granule


def get_processing_datetime(infile, granule):
    lines, logfile = read_log_file(infile)
    date_string = re.search('\d\d/\d\d/\d\d\d\d \d\d:\d\d:\d\d [A,P]M - INFO -\s*Input name\s*: {}'.format(granule),
                            lines)
    if not date_string:
        raise Exception('ERROR: No date_string found in {}'.format(infile))
    date_string = date_string.group(0)
    proc_date = re.search('\d\d/\d\d/\d{4} \d\d:\d\d:\d\d [A,P]M', date_string)
    if not proc_date:
        raise Exception('ERROR: No processing date found in {}'.format(date_string))
    proc_date = proc_date.group(0)
    date, time, period = proc_date.split()
    hour, minute, second = time.split(':')
    if 'PM' in period and int(hour) != 12:
        hour = int(hour) + 12
    proc_date = datetime(int(date[6:10]), int(date[0:2]), int(date[3:5]), int(hour), int(minute), int(second))
    logging.info('Processing_date {}'.format(proc_date))
    return proc_date


def read_log_file(infile):
    _, basename = os.path.split(infile)
    pol = get_pol(basename)
    if pol == 'VH':
        pol = 'VV'
    elif pol == 'HV':
        pol = 'HH'
    file_name = infile.replace('_{}.tif'.format(pol), '.log')
    if not os.path.exists(file_name):
        raise Exception('ERROR: Unable to find file {}'.format(file_name))
    logging.debug('Reading file {}'.format(file_name))
    with open(file_name, 'r') as f:
        lines = f.read()
    return lines, file_name


def get_pol(infile):
    if 'VV' in infile:
        pol = 'VV'
    elif 'VH' in infile:
        pol = 'VH'
    elif 'HH' in infile:
        pol = 'HH'
    elif 'HV' in infile:
        pol = 'HV'
    else:
        raise Exception('Could not determine polarization of file ' + infile)
    return pol


def fill_cfg(crs_wkt, prod_type, granule, proc_dt, hyp3_ver, gamma_ver, pix_x, pix_y):
    logging.info('Adding metadata')
    cfg = {'title': 'SAR RTC', 
           'institution': 'Alaska Satellite Facility (ASF)', 
           'mission': f'Sentinel-1{granule[2]}', 
           'crs_wkt': crs_wkt,
           'x_spacing': pix_x, 
           'y_spacing': pix_y, 
           'rtc_processing_date': datetime.strftime(proc_dt, '%Y%m%dT%H%M%S'),
           'source': f"ASF DAAC HyP3 {datetime.now().strftime('%Y')} using hyp3_gamma "
                     f'v{hyp3_ver} running GAMMA release {gamma_ver}. '
                     f'Contains modified Copernicus Sentinel data {granule[17:21]}, processed by ESA',
           'Conventions': 'CF - 1.6', 
           'references': 'asf.alaska.edu',
           'comment': 'This is an early prototype.'}

    #
    # we'll want this line for stacks:
    #    cfg['product_type'] = prod_type + ' stack'
    #
    #    cfg['feature_type'] = '????????'
    #
    #    Audit trail of process chain with time stamps
    #    cfg['history'] = f'Data acquired {granule[17:32]}; Processed to RTC at ASF on ' \
    #                      f'{datetime.strftime(proc_dt, '%Y%m%dT%H%M%S')} using Hyp3 v{hyp3_ver}; ' \
    #                      f'netCDF created {datetime.now().strftime('%Y%m%dT%H%M%S')} using ' \
    #                      f'hyp3_stacking {stacking_ver}'
    #
    #   The compliance checker says that I need the field _NCProperties.  But, when I try to use it I get:
    #   file_creation_stack = f'rioxarray={rioxarray.__version__}, xarray={xr.__version__}, rasterio={rio.__version__}'
    #   cfg['_NCProperties'] = file_creation_stack
    #       AttributeError: NetCDF: String match to name in use
    #

    cfg_x = {'axis': 'X', 'units': 'm', 'standard_name': 'projection_x_coordinate', 'long_name': 'Easting'}
    cfg_y = {'axis': 'Y', 'units': 'm', 'standard_name': 'projection_y_coordinate', 'long_name': 'Northing'}

    return cfg, cfg_x, cfg_y


def get_science_code(infile):
    readme = os.path.join(os.path.dirname(infile), 'README.txt')
    with open(readme, 'r') as f:
        lines = f.readlines()
    sc_name = 'UNKNOWN'
    hyp3_ver = 'UNKNOWN'
    gamma_ver = 'UNKNOWN'
    for line in lines:
        if 'RTC' in line and 'GAMMA' in line:
            sc_name = 'RTC GAMMA'
        elif 'RTC' in line and 'S1TBX' in line:
            sc_name = 'RTC S1TBX'

        if 'HYP3' in line and 'software version' in line:
            obj = re.search('\d+\.\d+\.*\d*', line)
            hyp3_ver = obj.group(0)

        if 'release' in line and '{}'.format(sc_name[0]) in line:
            obj = re.search('\d{8}', line)
            gamma_ver = obj.group(0)

    return hyp3_ver, gamma_ver


def create_dB(data, input_type):
    if 'amp' in input_type:
        power_data = data * data
        decibel_data = 10 * np.log(power_data)
    elif 'power' in input_type:
        decibel_data = 10 * np.log(data)
    else:
        raise Exception(f'Unknown input_type: {input_type}')
    return decibel_data


def scale_data(backscatter, scene, output_scale):
    power_data = backscatter
    if 'amp' in scene['scale']:
        power_data = backscatter * backscatter

    if 'power' in output_scale:
        scaled_data = power_data
    elif 'amp' in output_scale:
        scaled_data = backscatter
    elif 'db' in output_scale.lower():
        scaled_data = 10 * np.log(power_data)
    else:
        raise Exception(f'Unknown output scale {output_scale}')
    return scaled_data


def check_for_all_zeros(data):
    is_all_zero = np.all((data == 0))
    if is_all_zero:
        logging.warning('Data array is all zeros!')


def initialize_metadata(data_array, cfg, cfg_x, cfg_y):
    for key in cfg:
        data_array.attrs[key] = cfg[key]
    for key in cfg_x:
        data_array.x.attrs[key] = cfg_x[key]
    for key in cfg_y:
        data_array.y.attrs[key] = cfg_y[key]


#    for key in cfg_p:
#        data_array.variables[crs_tmp].attrs[key] = cfg_p[key]


def single2netcdf(prod_type, outfile, infile, output_scale=None, resolution=None):
    logging.info('single2netcdf: {} {} {} {}'.format(prod_type, outfile, infile, output_scale))
    dataset, image_dts, proc_dt, x_extent, y_extent, granule = get_dataset(infile)

    # get initial pixel size for dataset
    pix_x = dataset.transform[0]
    pix_y = dataset.transform[4]

    # resample if necessary    
    if resolution:
        if pix_x < resolution:
            res = int(resolution)
            root, unused = os.path.splitext(os.path.basename(infile))
            resampled_file = f'{root}_{res}.tif'
            logging.info(f'Resampling {infile} to file {resampled_file}')
            gdal.Translate(resampled_file, infile, xRes=resolution, yRes=resolution, resampleAlg='cubic')
            pix_x = resolution
            pix_y = -1 * resolution
            dataset.close()
            dataset = rio.open(resampled_file)
        elif pix_x >= resolution:
            logging.warning(f'Desired output resolution less than original  ({resolution} vs {pix_x})')
            logging.warning('No resampling performed')
    else:
        logging.info('Skipping resample step')

    no_data_value = dataset.nodata
    epsg_code = dataset.crs.to_epsg()
    crs_wkt = dataset.crs.wkt
    crs = pycrs.parse.from_ogc_wkt(crs_wkt)
    crs_name = crs.proj.name.ogc_wkt.lower()
    hyp3_ver, gamma_ver = get_science_code(infile)
    cfg, cfg_x, cfg_y = fill_cfg(crs_wkt, prod_type, granule, proc_dt, hyp3_ver, gamma_ver, pix_x, pix_y)

    logging.debug('Config is {}'.format(cfg))
    logging.info('Processing file {}'.format(infile))

    scene = parse_asf_rtc_name(os.path.basename(infile))

    backscatter = dataset.read(1)
    backscatter = scale_data(backscatter, scene, output_scale)
    backscatter = np.ma.masked_invalid(backscatter, copy=True)
    check_for_all_zeros(backscatter)

    x_coords = np.arange(x_extent[0], x_extent[1], pix_x)
    y_coords = np.arange(y_extent[0], y_extent[1], pix_y)

    # Deal with inconsistent rounding in arange!
    logging.info(f'Dataset.width = {dataset.width}')
    logging.info(f'Dataset.height = {dataset.height}')

    if len(x_coords) > dataset.width:
        x_coords = x_coords[0:dataset.width]
    if len(y_coords) > dataset.width:
        y_coords = y_coords[0:dataset.height]

    data_array = xr.Dataset({
        'y': y_coords,
        'x': x_coords,
        crs_name: crs_name,
        'backscatter': (('y', 'x'), backscatter.filled(0.0), {
            '_FillValue': no_data_value,
            'grid_mapping': crs_name,
            'long_name': 'SAR RTC',
            'radiometry': scene['radiometry'],
            'scaling': scene['scale'],
            'standard_name': 'SAR',
            'product_name': os.path.basename(infile),
            'granule': granule,
            'acquisition_start_time': granule[17:32],
            'acquisition_end_time': granule[33:48],
            'polarization': get_pol(infile),
            'radiation_frequency': 3.0 / 0.555,
            'radiation_frequency_unit': 'GHz',
            'radiation_wavelength': 3.0 / ((3.0 / 0.555) * 10),
            'radiation_wavelength_unit': 'm',
            'sensor_band_identifier': 'C'
        })
    })

    # Set the metadata based upon the config dictionaries
    initialize_metadata(data_array, cfg, cfg_x, cfg_y)

    # Set the CRS for this dataset
    data_array = data_array.rio.set_crs(epsg_code)
    dsp = data_array.rio.write_crs(epsg_code)

    # Fix the CRS for this dataset
    dsp.variables[crs_name].attrs['spatial_ref'] = crs_wkt
    dsp.variables[crs_name].attrs['crs_wkt'] = crs_wkt

    # FIXME: Want to delete the backscatter coordinates, but not sure how?
    #    del data_array.variables['backscatter'].attrs['coordinates']
    logging.info('Writing file {}'.format(outfile))
    dsp.to_netcdf(outfile, encoding={
        'x': {'_FillValue': None},
        'y': {'_FillValue': None},
    })
    logging.info('Successful Completion!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='single2netcdf.py',
                                     description='Convert an RTC stack from .tif format into netCDF4',
                                     epilog='The log and README files must be in the same directory as the .tif')

    parser.add_argument('infile', help='Name of input geotiff file')
    parser.add_argument('outfile', help='Name of output netcdf file')
    parser.add_argument('product_type', help='Type of data being stacked (RTC or INSAR)', metavar='InputType',
                        choices=['INSAR', 'RTC'])
    parser.add_argument('-o', '--output_scale', help='Output scale type\n', choices=['power', 'amp', 'dB'],
                        default='power')
    parser.add_argument('-r', '--resolution', help='Desired output resolution', type=float)
    args = parser.parse_args()

    logFile = 'single2netcdf_{}.log'.format(os.getpid())
    logging.basicConfig(filename=logFile, format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.DEBUG)
    logging.getLogger().addHandler(logging.StreamHandler())
    logging.info('Starting run')

    single2netcdf(args.product_type, args.outfile, args.infile, args.output_scale, args.resolution)
