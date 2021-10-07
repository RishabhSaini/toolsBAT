import numpy as np
import time
import datetime
from astropy.time import Time
from urllib.request import urlopen
from bs4 import BeautifulSoup
import math
import ephem
from astropy.time import Time
import requests
import urllib.parse
import os, sys, json


#Time funcs
def unix2mjd(unixtime):
    """Convert UNIX timestamp into MJD"""
    dt = Time(unixtime,format='unix',scale='utc')
    return dt.mjd

def get_attitude(trigtime_obj):
  unixtime=time.mktime(trigtime_obj.timetuple())
  mjd=unix2mjd(unixtime)
  base = 'https://www.swift.psu.edu/operations/afst_json.php?mjdstart='
  url = base+str(mjd)
  r = requests.get(url = url)
  settle_time = datetime.datetime.strptime(json.loads(r.text)['api_data']['entries'][0]['api_data']['settle'], "%Y-%m-%d %H:%M:%S")
  if trigtime_obj<settle_time:
    print('WARNING: Trigger time in slew. ObsSchedule pointing info not reliable. Must use attitude file.')
    return -1,-1,-1
  else:
    ra=float(json.loads(r.text)['api_data']['entries'][0]['api_data']['ra'])
    dec=float(json.loads(r.text)['api_data']['entries'][0]['api_data']['dec'])
    roll=float(json.loads(r.text)['api_data']['entries'][0]['api_data']['roll'])
    return ra,dec,roll


#TLE stuff
def getDataFromTLE(datetime, tleLatOffset=0, tleLonOffset=0.21):
    # Get TLE and parse
    url = "https://celestrak.com/satcat/tle.php?CATNR=28485"
    data = urlopen(url)
    tle_raw=data.read()
    clean_tle = ''.join(BeautifulSoup(tle_raw, "html.parser").stripped_strings)
    tle_obj = clean_tle.split('\r\n')

    # Print age of TLE
    year = "20"+tle_obj[1][18:20]
    day = tle_obj[1][20:32]
    #print("\nTLE most recently updated "+year+"DOY"+ day)
    
    # Create spacecraft instance
    Swift = ephem.readtle(tle_obj[0],tle_obj[1],tle_obj[2])
    
    # Create observer
    observer_Swift = ephem.Observer()
    observer_Swift.lat = '0'
    observer_Swift.long = '0'
    
    # Initialize lists
    lat = []
    lon = []
    elevation = []
    
    # Iterate through time; compute predicted locations
    observer_Swift.date = ephem.date(datetime)
    
    try:
        Swift.compute(observer_Swift)
    except:
        print('Cant compute observer, maybe TLE too old?')
        return False, False, False

    lat.append(np.degrees(Swift.sublat.znorm))
    lon.append(np.degrees(Swift.sublong.norm))
    elevation.append(Swift.elevation)
    # Correct for inaccuracy
    lon = np.array(lon) - tleLonOffset
    lat = np.array(lat) - tleLatOffset
    
    if len(lon) == 1: # Return single value instead of array if length is 1
        lon = lon[0]
        lat = lat[0]
        elevation = elevation[0]

    return lon, lat, elevation

def deg2dm(deg):
    sign = np.sign(deg)
    deg = np.abs(deg)
    d = np.floor(deg)
    m = (deg - d) * 60
    return int(sign*d), m

def getGeoCenter(datetime, lon, lat):
    # Define the observer to be at the location of the spacecraft
    observer = ephem.Observer()

    # Convert the longitude to +E (-180 to 180)
    if lon > 180:
       lon = (lon % 180) - 180

    lon_deg, lon_min = deg2dm(lon)
    lat_deg, lat_min = deg2dm(lat)

    lon_string = '%s:%s' % (lon_deg, lon_min)
    lat_string = '%s:%s' % (lat_deg, lat_min)

    observer.lon = lon_string
    observer.lat = lat_string
    
    # Set the time of the observations
    observer.date = ephem.date(datetime)
    
    # Get the ra and dec (in radians) of the point in the sky at altitude = 90 (directly overhead)
    ra_zenith_radians, dec_zenith_radians = observer.radec_of('0', '90')
    
    # convert the ra and dec to degrees
    ra_zenith = np.degrees(ra_zenith_radians)
    dec_zenith = np.degrees(dec_zenith_radians)
    
    ra_geocenter =  (ra_zenith+180) % 360
    dec_geocenter = -1 * dec_zenith
    
    return ra_geocenter, dec_geocenter

def getearthsatpos(datetime):
    tleLonOffset=0.21
    tleLatOffset = 0

    try:
        lon, lat, elevation= getDataFromTLE(datetime, tleLatOffset=tleLatOffset, tleLonOffset=tleLonOffset)
    except:
        return False, False, False

    if lon == False and lat == False and elevation == False:
        return False, False, False

    # Get the geo center coordinates in ra and dec
    ra_geocenter, dec_geocenter = getGeoCenter(datetime, lon, lat)

    EARTH_RADIUS = 6378.140 * 1000 #in meters
    dtor = math.pi/180
    elev = elevation + EARTH_RADIUS
    earthsize_rad = np.arcsin(EARTH_RADIUS/elev)/dtor

    return ra_geocenter,dec_geocenter, earthsize_rad

#Plotting
def makeEarthContour(ra,dec,radius):
    thetas = np.linspace(0, -2*np.pi, 200)
    ras = radius * np.cos(thetas)
    decs = radius * np.sin(thetas)
    contour = np.c_[ras,decs]
    Earthcont = project_footprint(contour, ra, dec, 0)
    return Earthcont