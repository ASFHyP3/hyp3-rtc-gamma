"""Process Sentinel-1 data into interferograms using GAMMA"""

import argparse
import glob
import logging
import os
import re
import shutil
import sys
from datetime import datetime
from secrets import token_hex

from hyp3lib import GranuleError
from hyp3lib.SLC_copy_S1_fullSW import SLC_copy_S1_fullSW
from hyp3lib.execute import execute
from hyp3lib.getParameter import getParameter
from hyp3lib.get_orb import downloadSentinelOrbitFile
from hyp3lib.makeAsfBrowse import makeAsfBrowse
from hyp3lib.par_s1_slc_single import par_s1_slc_single
from lxml import etree

from hyp3_gamma.insar.create_metadata_insar_gamma import create_readme_file
from hyp3_gamma.insar.getDemFileGamma import get_dem_file_gamma
from hyp3_gamma.insar.interf_pwr_s1_lt_tops_proc import interf_pwr_s1_lt_tops_proc
from hyp3_gamma.insar.unwrapping_geocoding import unwrapping_geocoding

log = logging.getLogger(__name__)


def get_bursts(mydir, name):
    back = os.getcwd()
    os.chdir(os.path.join(mydir, "annotation"))

    total_bursts = None
    time = []
    for myfile in os.listdir("."):
        if name in myfile:
            root = etree.parse(myfile)
            for coord in root.iter('azimuthAnxTime'):
                time.append(float(coord.text))
            for count in root.iter('burstList'):
                total_bursts = int(count.attrib['count'])

    os.chdir(back)

    return time, total_bursts


def get_burst_overlaps(reference_dir, secondary_dir):
    log.info("Calculating burst overlaps; in directory {}".format(os.getcwd()))
    burst_tab1 = "%s_burst_tab" % reference_dir[17:25]
    burst_tab2 = "%s_burst_tab" % secondary_dir[17:25]

    with open(burst_tab1, "w") as f1:
        with open(burst_tab2, "w") as f2:
            for name in ['001.xml', '002.xml', '003.xml']:
                time1, total_bursts1 = get_bursts(reference_dir, name)
                log.info("total_bursts1, time1 {} {}".format(total_bursts1, time1))
                time2, total_bursts2 = get_bursts(secondary_dir, name)
                log.info("total_bursts2, time2 {} {}".format(total_bursts2, time2))
                cnt = 1
                start1 = 0
                start2 = 0
                found = 0
                x = time1[0]
                for y in time2:
                    if abs(x - y) < 0.20:
                        log.info("Found burst match at 1 %s" % cnt)
                        found = 1
                        start1 = 1
                        start2 = cnt
                    cnt += 1

                if found == 0:
                    y = time2[0]
                    cnt = 1
                    for x in time1:
                        if abs(x - y) < 0.20:
                            log.info("Found burst match at %s 1" % cnt)
                            start1 = cnt
                            start2 = 1
                        cnt += 1

                size1 = total_bursts1 - start1 + 1
                size2 = total_bursts2 - start2 + 1

                if size1 > size2:
                    size = size2
                else:
                    size = size1

                f1.write("%s %s\n" % (start1, start1 + size - 1))
                f2.write("%s %s\n" % (start2, start2 + size - 1))

    return burst_tab1, burst_tab2


def get_copol(granule_name):
    polarization = granule_name[14:16]
    if polarization in ['SV', 'DV']:
        return 'vv'
    if polarization in ['SH', 'DH']:
        return 'hh'
    raise GranuleError(f'Cannot determine co-polarization of granule {granule_name}')


def least_precise_orbit_of(orbits):
    if any([orb is None for orb in orbits]):
        return 'O'
    if any(['RESORB' in orb for orb in orbits]):
        return 'R'
    return 'P'


def timedetla_in_days(delta):
    seconds_in_a_day = 60 * 60 * 24
    total_seconds = abs(delta.total_seconds())
    return round(total_seconds/seconds_in_a_day)


def get_product_name(reference_name, secondary_name, orbit_files, pixel_spacing=80, water_mask=False):
    plat1 = reference_name[2]
    plat2 = secondary_name[2]

    datetime1 = reference_name[17:32]
    datetime2 = secondary_name[17:32]

    ref_datetime = datetime.strptime(datetime1, '%Y%m%dT%H%M%S')
    sec_datetime = datetime.strptime(datetime2, '%Y%m%dT%H%M%S')
    days = timedetla_in_days(ref_datetime - sec_datetime)

    pol1 = reference_name[15:16]
    pol2 = secondary_name[15:16]
    orb = least_precise_orbit_of(orbit_files)
    mask = 'w' if water_mask else 'u'
    product_id = token_hex(2).upper()

    return f'S1{plat1}{plat2}_{datetime1}_{datetime2}_{pol1}{pol2}{orb}{days:03}_INT{pixel_spacing}_G_{mask}eF_' \
           f'{product_id}'


def move_output_files(output, reference, prod_dir, long_output, include_los_displacement, include_look_vectors,
                      include_wrapped_phase, include_inc_map, water_masking, include_dem):
    inName = "{}.mli.geo.tif".format(reference)
    outName = "{}_amp.tif".format(os.path.join(prod_dir, long_output))
    shutil.copy(inName, outName)

    inName = "{}.cc.geo.tif".format(output)
    outName = "{}_corr.tif".format(os.path.join(prod_dir, long_output))
    if os.path.isfile(inName):
        shutil.copy(inName, outName)

    inName = "{}.vert.disp.geo.org.tif".format(output)
    outName = "{}_vert_disp.tif".format(os.path.join(prod_dir, long_output))
    shutil.copy(inName, outName)

    inName = "{}.adf.unw.geo.tif".format(output)
    outName = "{}_unw_phase.tif".format(os.path.join(prod_dir, long_output))
    shutil.copy(inName, outName)

    if include_wrapped_phase:
        inName = "{}.diff0.man.adf.geo.tif".format(output)
        outName = "{}_wrapped_phase.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    if include_dem:
        inName = "{}.dem.tif".format(output)
        outName = "{}_dem.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    if include_los_displacement:
        inName = "{}.los.disp.geo.org.tif".format(output)
        outName = "{}_los_disp.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    if include_inc_map:
        inName = "{}.inc.tif".format(output)
        outName = "{}_inc_map.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    if water_masking:
        inName = "final_mask.tif"
        outName = "{}_water_mask.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    if include_look_vectors:
        inName = "{}.lv_theta.tif".format(output)
        outName = "{}_lv_theta.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)
        inName = "{}.lv_phi.tif".format(output)
        outName = "{}_lv_phi.tif".format(os.path.join(prod_dir, long_output))
        shutil.copy(inName, outName)

    makeAsfBrowse("{}.diff0.man.adf.bmp.geo.tif".format(output),
                  "{}_color_phase".format(os.path.join(prod_dir, long_output)))

    makeAsfBrowse("{}.adf.unw.geo.bmp.tif".format(output),
                  "{}_unw_phase".format(os.path.join(prod_dir, long_output)))


def make_parameter_file(mydir, parameter_file_name, alooks, rlooks, dem_source):
    res = 20 * int(alooks)

    reference_date = mydir[:15]
    secondary_date = mydir[17:]
    reference_date_short = reference_date[:8]

    log.info("In directory {} looking for file with date {}".format(os.getcwd(), reference_date_short))
    reference_file = glob.glob("*%s*.SAFE" % reference_date)[0]
    secondary_file = glob.glob("*%s*.SAFE" % secondary_date)[0]

    parfile = f'{reference_date_short}.mli.par'
    erad_nadir = getParameter(parfile, 'earth_radius_below_sensor')
    erad_nadir = erad_nadir.split()[0]
    sar_to_earth_center = getParameter(parfile, 'sar_to_earth_center')
    sar_to_earth_center = sar_to_earth_center.split()[0]
    height = float(sar_to_earth_center) - float(erad_nadir)

    with open("baseline.log") as f:
        for line in f:
            if "estimated baseline perpendicular component" in line:
                # FIXME: RE is overly complicated here. this is two simple string splits
                t = re.split(":", line)
                s = re.split(r'\s+', t[1])
                baseline = float(s[1])

    back = os.getcwd()
    os.chdir(os.path.join(reference_file, "annotation"))

    utctime = None
    for myfile in os.listdir("."):
        if "001.xml" in myfile:
            root = etree.parse(myfile)
            for coord in root.iter('productFirstLineUtcTime'):
                utc = coord.text
                log.info("Found utc time {}".format(utc))
                t = utc.split("T")
                log.info("{}".format(t))
                s = t[1].split(":")
                log.info("{}".format(s))
                utctime = ((int(s[0]) * 60 + int(s[1])) * 60) + float(s[2])
    os.chdir(back)

    heading = None
    name = f'{reference_date[:8]}.mli.par'
    with open(name, "r") as f:
        for line in f:
            if "heading" in line:
                t = re.split(":", line)
                # FIXME: RE is overly complicated here. this is two simple string splits
                s = re.split(r'\s+', t[1])
                heading = float(s[1])

    reference_file = reference_file.replace(".SAFE", "")
    secondary_file = secondary_file.replace(".SAFE", "")

    with open(parameter_file_name, 'w') as f:
        f.write('Reference Granule: %s\n' % reference_file)
        f.write('Secondary Granule: %s\n' % secondary_file)
        f.write('Baseline: %s\n' % baseline)
        f.write('UTCtime: %s\n' % utctime)
        f.write('Heading: %s\n' % heading)
        f.write('Spacecraft height: %s\n' % height)
        f.write('Earth radius at nadir: %s\n' % erad_nadir)
        f.write('Range looks: %s\n' % rlooks)
        f.write('Azimuth looks: %s\n' % alooks)
        f.write('INSAR phase filter:  adf\n')
        f.write('Phase filter parameter: 0.6\n')
        f.write('Resolution of output (m): %s\n' % res)
        f.write('Range bandpass filter: no\n')
        f.write('Azimuth bandpass filter: no\n')
        f.write('DEM source: %s\n' % dem_source)
        f.write('DEM resolution (m): %s\n' % (res * 2))
        f.write('Unwrapping type: mcf\n')
        f.write('Unwrapping threshold: none\n')
        f.write('Speckle filtering: off\n')


def insar_sentinel_gamma(reference_file, secondary_file, rlooks=20, alooks=4, include_look_vectors=False,
                         include_los_displacement=False, include_wrapped_phase=False, include_inc_map=False,
                         water_masking=False, include_dem=False):
    log.info("\n\nSentinel-1 differential interferogram creation program\n")

    wrk = os.getcwd()
    reference_date = reference_file[17:32]
    reference_date_short = reference_file[17:25]
    secondary_date = secondary_file[17:32]
    secondary_date_short = secondary_file[17:25]
    igramName = "{}_{}".format(reference_date, secondary_date)

    if "IW_SLC__" not in reference_file:
        raise GranuleError(f'Reference file {reference_file} is not of type IW_SLC!')
    if "IW_SLC__" not in secondary_file:
        raise GranuleError(f'Secondary file {secondary_file} is not of type IW_SLC!')

    pol = get_copol(reference_file)
    log.info("Processing the {} polarization".format(pol))

    #  Ingest the data files into gamma format
    log.info("Starting par_S1_SLC")
    orbit_files = []
    for granule in (reference_file, secondary_file):
        orbit_file, _ = downloadSentinelOrbitFile(granule)
        par_s1_slc_single(granule, pol, os.path.abspath(orbit_file))
        orbit_files.append(orbit_file)

    #  Fetch the DEM file
    log.info("Getting a DEM file")
    dem_source = 'GLO-30'
    dem_pixel_size = int(alooks) * 40  # typically 160 or 80; IFG pixel size will be half the DEM pixel size (80 or 40)
    get_dem_file_gamma('big.dem', 'big.par', reference_file, pixel_size=dem_pixel_size, water_masking=water_masking)
    log.info("Got dem of type {}".format(dem_source))

    # Figure out which bursts overlap between the two swaths
    burst_tab1, burst_tab2 = get_burst_overlaps(reference_file, secondary_file)

    log.info("Finished calculating overlap - in directory {}".format(os.getcwd()))
    shutil.move(burst_tab1, reference_date_short)
    shutil.move(burst_tab2, secondary_date_short)

    # Mosaic the swaths together and copy SLCs over
    log.info("Starting SLC_copy_S1_fullSW.py")
    reference = reference_date_short
    secondary = secondary_date_short

    os.chdir(reference)
    SLC_copy_S1_fullSW(wrk, reference, "SLC_TAB", burst_tab1, mode=1, dem="big", dempath=wrk, raml=rlooks, azml=alooks)
    os.chdir("..")
    os.chdir(secondary)
    SLC_copy_S1_fullSW(wrk, secondary, "SLC_TAB", burst_tab2, mode=2, raml=rlooks, azml=alooks)
    os.chdir("..")

    # Interferogram creation, matching, refinement
    log.info("Starting interf_pwr_s1_lt_tops_proc.py 0")
    hgt = "DEM/HGT_SAR_{}_{}".format(rlooks, alooks)
    interf_pwr_s1_lt_tops_proc(reference, secondary, hgt, rlooks=rlooks, alooks=alooks, iterations=3, step=0)

    log.info("Starting interf_pwr_s1_lt_tops_proc.py 1")
    interf_pwr_s1_lt_tops_proc(reference, secondary, hgt, rlooks=rlooks, alooks=alooks, step=1)

    log.info("Starting interf_pwr_s1_lt_tops_proc.py 2")
    interf_pwr_s1_lt_tops_proc(reference, secondary, hgt, rlooks=rlooks, alooks=alooks, iterations=3, step=2)

    g = open("offsetfit3.log")
    offset = 1.0
    for line in g:
        if "final azimuth offset poly. coeff.:" in line:
            offset = line.split(":")[1]
    if float(offset) > 0.02:
        log.error("ERROR: Found azimuth offset of {}!".format(offset))
        sys.exit(1)
    else:
        log.info("Found azimuth offset of {}!".format(offset))

    output = reference_date_short + "_" + secondary_date_short

    log.info("Starting s1_coreg_overlap")
    execute(f"S1_coreg_overlap SLC1_tab SLC2R_tab {output} {output}.off.it {output}.off.it.corrected",
            uselogging=True)

    log.info("Starting interf_pwr_s1_lt_tops_proc.py 3")
    interf_pwr_s1_lt_tops_proc(reference, secondary, hgt, rlooks=rlooks, alooks=alooks, step=3)

    # Perform phase unwrapping and geocoding of results
    log.info("Starting phase unwrapping and geocoding")
    unwrapping_geocoding(reference, secondary, step="man", rlooks=rlooks, alooks=alooks)

    #  Generate metadata
    log.info("Collecting metadata and output files")

    os.chdir(wrk)

    # Move the outputs to the product directory
    pixel_spacing = int(alooks) * 20
    product_name = get_product_name(reference_file, secondary_file, orbit_files, pixel_spacing, water_masking)
    os.mkdir(product_name)
    move_output_files(output, reference, product_name, product_name, include_los_displacement, include_look_vectors,
                      include_wrapped_phase, include_inc_map, water_masking, include_dem)

    create_readme_file(reference_file, secondary_file, f'{product_name}/{product_name}.README.md.txt', pixel_spacing)

    execute(f"base_init {reference}.slc.par {secondary}.slc.par - - base > baseline.log", uselogging=True)
    make_parameter_file(igramName, f'{product_name}/{product_name}.txt', alooks, rlooks, dem_source)

    log.info("Done!!!")
    return product_name


def main():
    """Main entrypoint"""
    parser = argparse.ArgumentParser(
        prog='ifm_sentinel.py',
        description=__doc__,
    )
    parser.add_argument("reference", help="Reference input file")
    parser.add_argument("secondary", help="Secondary input file")
    parser.add_argument("-r", "--rlooks", default=20, help="Number of range looks (def=20)")
    parser.add_argument("-a", "--alooks", default=4, help="Number of azimuth looks (def=4)")
    parser.add_argument("-d", action="store_true", help="Add DEM file to product bundle")
    parser.add_argument("-i", action="store_true", help="Create incidence angle map")
    parser.add_argument("-l", action="store_true", help="Create look vector theta and phi files")
    parser.add_argument("-m", action="store_true", help="Create and apply water mask")
    parser.add_argument("-s", action="store_true", help="Create line of sight displacement file")
    parser.add_argument("-w", action="store_true", help="Create wrapped phase file")
    args = parser.parse_args()

    logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p', level=logging.INFO)

    insar_sentinel_gamma(args.reference, args.secondary, rlooks=args.rlooks, alooks=args.alooks,
                         include_look_vectors=args.l, include_los_displacement=args.s,
                         include_wrapped_phase=args.w, include_inc_map=args.i,
                         water_masking=args.m, include_dem=args.d)


if __name__ == "__main__":
    main()
