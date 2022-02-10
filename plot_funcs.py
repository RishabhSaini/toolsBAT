import numpy as np
import ephem
import datetime
import time
from shapely.geometry import Polygon
from shapely.geometry import Point
from urllib.request import Request, urlopen, urlretrieve
from bs4 import BeautifulSoup

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.offsetbox import OffsetImage, AnnotationBbox


def getPosFromTLE(trigtime):
	# Get TLE and parse
	url = "https://celestrak.com/satcat/tle.php?CATNR=28485"
	data = urlopen(url)
	tle_raw=data.read()
	clean_tle = ''.join(BeautifulSoup(tle_raw, "html.parser").stripped_strings)
	tle_obj = clean_tle.split('\r\n')

	# Print age of TLE
	year = "20"+tle_obj[1][18:20]
	day = tle_obj[1][20:32]
	print("\nTLE most recently updated "+year+"DOY"+ day)
	
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


	start = trigtime-datetime.timedelta(seconds=1000)
	end = trigtime+datetime.timedelta(seconds=1000)
	timestamps = [start + datetime.timedelta(seconds=x) for x in range(0, (end-start).seconds)]
	
	# Iterate through time; compute predicted locations
	for timestamp in timestamps:
	  observer_Swift.date = ephem.date(timestamp)
	  
	  try:
		  Swift.compute(observer_Swift)
	  except:
		  return False, False, False

	  lat.append(np.degrees(Swift.sublat.znorm))
	  lon.append(-(360-np.degrees(Swift.sublong.norm)))
	  elevation.append(Swift.elevation)

	
	if len(lon) == 1: # Return single value instead of array if length is 1
		lon = lon[0]
		lat = lat[0]
		elevation = elevation[0]

	return lon, lat, elevation


class BATSAA:
	def __init__(self,lat,long):
		self.lat = lat
		self.long = long
		self.points = np.array([[-8.50000000e+01, -1.99999935e+01],
	   [-8.49999999e+01, -21.5],
	   [-8.40000000e+01, -21.5],
	   [-8.30000000e+01, -21.5],
	   [-8.20000000e+01, -21.5],
	   [-8.10000000e+01, -21.5],
	   [-8.00000000e+01, -21.5],
	   [-7.90000000e+01, -21.5],
	   [-7.80000000e+01, -21.5],
	   [-7.70000000e+01, -21.5],
	   [-7.60000000e+01, -21.5],
	   [-7.50000000e+01, -21.5],
	   [-7.40000000e+01, -21.5],
	   [-7.30000000e+01, -21.5],
	   [-7.20000000e+01, -21.5],
	   [-7.10000000e+01, -21.5],
	   [-7.00000000e+01, -21.5],
	   [-6.90000000e+01, -21.5],
	   [-6.80000000e+01, -21.5],
	   [-6.70000000e+01, -21.5],
	   [-6.60000000e+01, -21.5],
	   [-6.50000000e+01, -21.5],
	   [-6.40000000e+01, -21.5],
	   [-6.30000000e+01, -21.5],
	   [-6.20000000e+01, -21.5],
	   [-6.10000000e+01, -21.5],
	   [-6.00000000e+01, -21.5],
	   [-5.90000000e+01, -21.5],
	   [-5.80000000e+01, -21.5],
	   [-5.70000000e+01, -21.5],
	   [-5.60000000e+01, -21.5],
	   [-5.50000000e+01, -21.5],
	   [-5.40000000e+01, -21.5],
	   [-5.30000000e+01, -21.5],
	   [-5.20000000e+01, -21.5],
	   [-5.10000000e+01, -21.5],
	   [-5.00000000e+01, -21.5],
	   [-4.90000000e+01, -21.5],
	   [-4.80000000e+01, -21.5],
	   [-4.70000000e+01, -21.5],
	   [-4.60000000e+01, -21.5],
	   [-4.50000000e+01, -21.5],
	   [-4.40000000e+01, -21.5],
	   [-4.30000000e+01, -21.5],
	   [-4.20000000e+01, -21.5],
	   [-4.10000000e+01, -21.5],
	   [-4.00000000e+01, -21.5],
	   [-3.90000000e+01, -21.5],
	   [-3.80000000e+01, -21.5],
	   [-3.70000000e+01, -21.5],
	   [-3.60000000e+01, -21.5],
	   [-3.50000000e+01, -21.5],
	   [-3.40000000e+01, -21.5],
	   [-3.30000000e+01, -21.5],
	   [-3.20000000e+01, -21.5],
	   [-3.10000000e+01, -21.5],
	   [-3.00000000e+01, -21.5],
	   [-2.90000000e+01, -21.5],
	   [-2.80000000e+01, -21.5],
	   [-2.70000000e+01, -21.5],
	   [-2.60000000e+01, -21.5],
	   [-2.50000000e+01, -21.5],
	   [-2.40000000e+01, -21.5],
	   [-2.30000000e+01, -21.5],
	   [-2.20000000e+01, -21.5],
	   [-2.10000000e+01, -21.5],
	   [-21.5, -21.5],
	   [-1.90000000e+01, -21.5],
	   [-1.80000000e+01, -21.5],
	   [-1.70000000e+01, -21.5],
	   [-1.60000000e+01, -21.5],
	   [-1.50000000e+01, -21.5],
	   [-1.40000000e+01, -21.5],
	   [-1.30000000e+01, -21.5],
	   [-1.20000000e+01, -21.5],
	   [-1.10000000e+01, -21.5],
	   [-1.00000000e+01, -21.5],
	   [-9.00000000e+00, -21.5],
	   [-8.00000000e+00, -21.5],
	   [-7.00000000e+00, -21.5],
	   [-6.00000000e+00, -21.5],
	   [-5.00000000e+00, -21.5],
	   [-4.00000000e+00, -21.5],
	   [-3.00000000e+00, -21.5],
	   [-2.00000000e+00, -21.5],
	   [-1.00000000e+00, -21.5],
	   [-4.15588488e-08, -1.90000000e+01],
	   [-3.35000038e-07, -1.90000000e+01],
	   [-1.00000000e+00, -1.80000003e+01],
	   [-2.00000000e+00, -1.80000001e+01],
	   [-3.00000000e+00, -1.80000000e+01],
	   [-4.00000000e+00, -1.80000000e+01],
	   [-5.00000000e+00, -1.80000000e+01],
	   [-5.00000058e+00, -1.80000000e+01],
	   [-6.00000000e+00, -1.70000006e+01],
	   [-7.00000000e+00, -1.70000001e+01],
	   [-7.00000030e+00, -1.70000000e+01],
	   [-8.00000000e+00, -1.60000003e+01],
	   [-9.00000000e+00, -1.60000001e+01],
	   [-1.00000000e+01, -1.60000000e+01],
	   [-1.00000005e+01, -1.60000000e+01],
	   [-1.10000000e+01, -1.50000005e+01],
	   [-1.20000000e+01, -1.50000000e+01],
	   [-1.30000000e+01, -1.50000000e+01],
	   [-1.30000000e+01, -1.50000000e+01],
	   [-1.40000000e+01, -1.40000000e+01],
	   [-1.40000001e+01, -1.40000000e+01],
	   [-1.50000000e+01, -1.30000001e+01],
	   [-1.60000000e+01, -1.30000000e+01],
	   [-1.70000000e+01, -1.30000000e+01],
	   [-1.70000000e+01, -1.30000000e+01],
	   [-1.80000000e+01, -1.20000000e+01],
	   [-1.90000000e+01, -1.20000000e+01],
	   [-1.90000000e+01, -1.20000000e+01],
	   [-1.90000000e+01, -1.19999999e+01],
	   [-1.80000001e+01, -1.10000000e+01],
	   [-1.80000000e+01, -1.10000000e+01],
	   [-1.70000000e+01, -1.00000000e+01],
	   [-1.70000001e+01, -9.00000000e+00],
	   [-1.80000000e+01, -8.00000013e+00],
	   [-1.90000000e+01, -8.00000002e+00],
	   [-1.90000001e+01, -8.00000000e+00],
	   [-21.5, -7.00000014e+00],
	   [-2.10000000e+01, -7.00000003e+00],
	   [-2.20000000e+01, -7.00000001e+00],
	   [-2.20000000e+01, -7.00000000e+00],
	   [-2.30000000e+01, -6.00000003e+00],
	   [-2.40000000e+01, -6.00000001e+00],
	   [-2.50000000e+01, -6.00000000e+00],
	   [-2.50000000e+01, -6.00000000e+00],
	   [-2.60000000e+01, -5.00000002e+00],
	   [-2.70000000e+01, -5.00000001e+00],
	   [-2.70000000e+01, -5.00000000e+00],
	   [-2.80000000e+01, -4.00000004e+00],
	   [-2.80000000e+01, -4.00000000e+00],
	   [-2.90000000e+01, -3.00000003e+00],
	   [-3.00000000e+01, -3.00000002e+00],
	   [-3.00000000e+01, -3.00000000e+00],
	   [-3.00000003e+01, -2.00000000e+00],
	   [-3.10000000e+01, -1.00000032e+00],
	   [-3.10000003e+01, -1.00000000e+00],
	   [-3.20000000e+01, -2.79999995e-07],
	   [-3.30000000e+01, -4.84444485e-08],
	   [-3.40000000e+01, -1.17647119e-08],
	   [-3.50000000e+01, -1.62000049e-08],
	   [-3.50000001e+01,  0.00000000e+00],
	   [-3.60000000e+01,  9.99999907e-01],
	   [-3.70000000e+01,  9.99999984e-01],
	   [-3.80000000e+01,  9.99999988e-01],
	   [-3.90000000e+01,  9.99999888e-01],
	   [-4.00000000e+01,  9.99999982e-01],
	   [-4.10000000e+01,  9.99999993e-01],
	   [-4.20000000e+01,  9.99999995e-01],
	   [-4.30000000e+01,  9.99999995e-01],
	   [-4.40000000e+01,  9.99999993e-01],
	   [-4.50000000e+01,  9.99999987e-01],
	   [-4.60000000e+01,  9.99999990e-01],
	   [-4.70000000e+01,  9.99999993e-01],
	   [-4.80000000e+01,  9.99999995e-01],
	   [-4.90000000e+01,  9.99999994e-01],
	   [-5.00000000e+01,  9.99999994e-01],
	   [-5.10000000e+01,  9.99999977e-01],
	   [-5.20000000e+01,  9.99999900e-01],
	   [-5.30000000e+01,  9.99999943e-01],
	   [-5.40000000e+01,  9.99999969e-01],
	   [-5.50000000e+01,  9.99999854e-01],
	   [-5.60000000e+01,  9.99999867e-01],
	   [-5.69999999e+01,  0.00000000e+00],
	   [-5.70000000e+01, -7.33845695e-09],
	   [-5.80000000e+01, -1.05333271e-08],
	   [-5.90000000e+01, -1.47777826e-08],
	   [-5.90000007e+01,  0.00000000e+00],
	   [-6.00000000e+01,  9.99999315e-01],
	   [-6.10000000e+01,  9.99999901e-01],
	   [-6.19999999e+01,  0.00000000e+00],
	   [-6.20000000e+01, -2.48888909e-08],
	   [-6.30000000e+01, -4.80000040e-08],
	   [-6.40000000e+01, -1.93000005e-07],
	   [-6.49999998e+01, -1.00000000e+00],
	   [-6.50000000e+01, -1.00000002e+00],
	   [-6.60000000e+01, -1.00000004e+00],
	   [-6.70000000e+01, -2.00000000e+00],
	   [-6.70000000e+01, -2.00000001e+00],
	   [-6.80000000e+01, -2.00000002e+00],
	   [-6.90000000e+01, -2.00000011e+00],
	   [-6.99999999e+01, -3.00000000e+00],
	   [-7.00000000e+01, -3.00000002e+00],
	   [-7.10000000e+01, -3.00000007e+00],
	   [-7.19999999e+01, -4.00000000e+00],
	   [-7.20000000e+01, -4.00000001e+00],
	   [-7.30000000e+01, -4.00000003e+00],
	   [-7.40000000e+01, -4.00000021e+00],
	   [-7.49999998e+01, -5.00000000e+00],
	   [-7.50000000e+01, -5.00000002e+00],
	   [-7.60000000e+01, -6.00000000e+00],
	   [-7.60000000e+01, -6.00000002e+00],
	   [-7.70000000e+01, -6.00000025e+00],
	   [-7.79999998e+01, -7.00000000e+00],
	   [-7.80000000e+01, -7.00000012e+00],
	   [-7.89999999e+01, -8.00000000e+00],
	   [-7.90000000e+01, -8.00000001e+00],
	   [-8.00000000e+01, -8.00000005e+00],
	   [-8.10000000e+01, -9.00000000e+00],
	   [-8.10000000e+01, -9.00000002e+00],
	   [-8.20000000e+01, -1.00000000e+01],
	   [-8.20000000e+01, -1.10000000e+01],
	   [-8.20000000e+01, -1.10000000e+01],
	   [-8.30000000e+01, -1.10000002e+01],
	   [-8.39999998e+01, -1.20000000e+01],
	   [-8.40000000e+01, -1.30000000e+01],
	   [-8.40000000e+01, -1.30000001e+01],
	   [-8.49999999e+01, -1.40000000e+01],
	   [-8.50000000e+01, -1.50000000e+01],
	   [-8.50000000e+01, -1.50000004e+01],
	   [-8.59999996e+01, -1.60000000e+01],
	   [-8.59999993e+01, -1.70000000e+01],
	   [-8.59999999e+01, -1.80000000e+01],
	   [-8.59999935e+01, -1.90000000e+01],
	   [-8.50000000e+01, -20],
	   [-8.50000000e+01, -21.5]])
		self.saapoly = Polygon(self.points)

	def insaa(self):
		# '''For a given time, are we inside the BAT SAA polygon'''
		# if not self.eph:
		#     self.eph = STKReadEph(latestephem(self.year,self.day))
		# index = self.eph.ephindex(utime)

		# self.long = self.eph.long[index]
		# self.lat = self.eph.lat[index]

		return self.saapoly.contains(Point(self.long,self.lat))

def calc_mcilwain_l(latitude, longitude):
	"""Estimate the McIlwain L value given the latitude (-30, +30) and 
	East Longitude.  This uses a cubic polynomial approximation to the full 
	calculation and is similar to the approach used by the GBM FSW.
   
	Args:
		latitude (np.array): Latitude in degrees from -180 to 180
		longitude (np.array): East longitude in degrees from 0 to 360
	
	Returns:
		np.array: McIlwain L value
	"""

	latitude = np.asarray([latitude])
	longitude = np.asarray([longitude])
	orig_shape = latitude.shape
	latitude = latitude.flatten()
	longitude = longitude.flatten()
	# numPts = latitude.shape[0]
	#coeffs_file = os.path.join(os.path.dirname(__file__),'McIlwainL_Coeffs.npy')
	#poly_coeffs = np.load(coeffs_file)
	
	poly_coeffs = np.array([[ 1.11534e+00, -9.08140e-03,  4.25143e-04,  9.38231e-07],
	   [ 1.10493e+00, -9.14453e-03,  4.53037e-04,  1.86429e-07],
	   [ 1.09338e+00, -8.81792e-03,  4.75715e-04, -4.69379e-07],
	   [ 1.08172e+00, -8.13947e-03,  4.92098e-04, -9.68546e-07],
	   [ 1.06871e+00, -7.18035e-03,  5.02030e-04, -1.59499e-06],
	   [ 1.05403e+00, -6.41215e-03,  5.07259e-04, -2.36934e-06],
	   [ 1.03835e+00, -6.24044e-03,  5.15952e-04, -3.00787e-06],
	   [ 1.02520e+00, -6.45000e-03,  5.27692e-04, -3.56748e-06],
	   [ 1.01458e+00, -6.63289e-03,  5.41228e-04, -4.17381e-06],
	   [ 1.00585e+00, -6.54267e-03,  5.53382e-04, -4.81297e-06],
	   [ 1.00018e+00, -6.27450e-03,  5.58926e-04, -5.26692e-06],
	   [ 9.98888e-01, -5.97149e-03,  5.56298e-04, -5.47440e-06],
	   [ 1.00365e+00, -5.81604e-03,  5.44148e-04, -5.39764e-06],
	   [ 1.01316e+00, -5.77722e-03,  5.25601e-04, -5.09640e-06],
	   [ 1.02513e+00, -5.68398e-03,  5.02120e-04, -4.59079e-06],
	   [ 1.03589e+00, -5.26021e-03,  4.76949e-04, -4.08567e-06],
	   [ 1.04276e+00, -4.44838e-03,  4.53280e-04, -3.59028e-06],
	   [ 1.04528e+00, -3.31088e-03,  4.32882e-04, -3.11920e-06],
	   [ 1.04655e+00, -2.06559e-03,  4.17320e-04, -2.60663e-06],
	   [ 1.04927e+00, -9.81826e-04,  4.04499e-04, -1.95729e-06],
	   [ 1.05376e+00, -1.97174e-04,  3.94790e-04, -1.05819e-06],
	   [ 1.05897e+00,  5.03682e-04,  3.88200e-04, -8.00078e-08],
	   [ 1.06542e+00,  1.31589e-03,  3.85567e-04,  8.21309e-07],
	   [ 1.07396e+00,  2.21685e-03,  3.90154e-04,  1.74743e-06],
	   [ 1.08446e+00,  3.21422e-03,  4.04193e-04,  2.70925e-06],
	   [ 1.09830e+00,  4.35729e-03,  4.25459e-04,  3.62526e-06],
	   [ 1.11814e+00,  5.67198e-03,  4.47946e-04,  4.33789e-06],
	   [ 1.14611e+00,  7.04419e-03,  4.56505e-04,  4.46156e-06],
	   [ 1.17719e+00,  7.88691e-03,  4.41125e-04,  4.02487e-06],
	   [ 1.19453e+00,  7.50620e-03,  4.08137e-04,  3.33238e-06],
	   [ 1.17747e+00,  5.15525e-03,  3.84831e-04,  3.31511e-06],
	   [ 1.14391e+00,  1.55877e-03,  3.67130e-04,  3.55329e-06],
	   [ 1.12344e+00, -1.84304e-03,  3.52753e-04,  3.32335e-06],
	   [ 1.11912e+00, -4.75332e-03,  3.50950e-04,  2.90764e-06],
	   [ 1.12162e+00, -7.00585e-03,  3.65812e-04,  2.38347e-06],
	   [ 1.12179e+00, -8.44449e-03,  3.92892e-04,  1.71072e-06]])

	longitude[longitude < 0.0] += 360.0
	longitude[longitude == 360.0] = 0.0

	bad_idx = (latitude < -30.0) | (latitude > 30.0) | (longitude < 0.0) | (
				longitude >= 360.0)
	if np.sum(bad_idx) != 0:
		raise ValueError(
			'Out of range coordinates for McIlwain L for {0} locations'.format(
				np.sum(bad_idx)))

	idx = np.asarray((longitude / 10.0).astype(int))
	idx2 = np.asarray(idx + 1)
	idx2[idx2 >= 36] = 0
	idx2 = idx2.astype(int)

	longitude_left = 10.0 * idx
	f = (longitude - longitude_left) / 10.0  # interpolation weight, 0 to 1

	try:
		num_pts = len(latitude)
	except:
		num_pts = 1
	mc_l = np.zeros(num_pts)
	for i in range(num_pts):
		mc_l[i] = (1.0 - f[i]) * (
				poly_coeffs[idx[i], 0] + poly_coeffs[idx[i], 1] * latitude[i] +
				poly_coeffs[idx[i], 2] *
				latitude[i] ** 2 + poly_coeffs[idx[i], 3] * latitude[i] ** 3) + \
				  f[i] * (
						  poly_coeffs[idx2[i], 0] + poly_coeffs[idx2[i], 1] *
						  latitude[i] + poly_coeffs[idx2[i], 2] *
						  latitude[i] ** 2 + poly_coeffs[idx2[i], 3] *
						  latitude[i] ** 3)
	mc_l = mc_l.reshape(orig_shape)
	return np.squeeze(mc_l)

def mcilwain_map(lat_range, lon_range, m, ax, saa_mask=False, color=None,
				 alpha=0.5, **kwargs):
	"""Plot the McIlwain L heatmap on a Basemap
	
	Parameters:
	-----------
	lat_range: (float, float)
		The latitude range
	lon_range: (float, float)
		The longitude range
	m: Basemap
		The basemap references
	ax: matplotlib.axes
		The plot axes references
	saa_mask: bool, optional
		If True, mask out the SAA from the heatmap.  Default is False.
	color: str, optional
		The color of the heatmap
	alpha: float, optional
		The alpha opacity of the heatmap
	kwargs: optional
		Other plotting keywords
	
	Returns:
	-----------
	image: QuadContourSet
		The heatmap plot object
	"""
	# do create an array on the earth
	lat_array = np.linspace(*lat_range, 108)
	lon_array = np.linspace(*lon_range, 720)
	LAT, LON = np.meshgrid(lat_array, lon_array)
	# convert to projection coordinates
	mLON, mLAT = m(LON, LAT)
	# mcilwain l over the grid
	mcl = calc_mcilwain_l(LAT, LON)

	# if we want to mask out the SAA
	if saa_mask:
		saa_path = saa_mask.get_path()
		mask = saa_path.contains_points(
			np.array((mLON.ravel(), mLAT.ravel())).T)
		mcl[mask.reshape(mcl.shape)] = 0.0

	# do the plot  
	levels = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
	image = ax.contourf(mLON, mLAT, mcl, levels=levels, alpha=alpha, **kwargs)
	return image

def orbitfromsao(saofile):
	orbitdat=saofile[1].data
	times=[]
	lats=[]
	lons=[]
	alts=[]
	earthras=[]
	earthdecs=[]
	for i in range(len(orbitdat)):
		met = orbitdat[i]['TIME_ADJ']
		lat = orbitdat[i]['SAT_LAT']
		lon = orbitdat[i]['SAT_LON']
		alt = orbitdat[i]['SAT_ALT']
		earthra= orbitdat[i]['EARTH_RA']
		earthdec = orbitdat[i]['EARTH_DEC']
		times.append(met)
		lats.append(lat)
		lons.append(lon)
		alts.append(alt)
		earthras.append(earthra)
		earthdecs.append(earthdec)

	return times, lons,lats, alts, earthras, earthdecs 

def unixtime2sc(unixtime):
	"""Convert Unix time to Spacecraft time"""
	return int(unixtime) - time.mktime((2001,1,1,0,0,0,0,0,0))+(unixtime-int(unixtime))

def earth_line(lat, lon, m, ax, color='black', alpha=0.4):
    """Plot a line on the Earth (e.g. orbit)
    
    Parameters:
    -----------
    lat: np.array
        Array of latitudes
    lon: np.array
        Array of longitudes
    m: Basemap
        The basemap references
    color: str, optional
        The color of the lines
    alpha: float, optional
        The alpha opacity of line
    kwargs: optional
        Other plotting keywords
    
    Returns:
    -----------
    refs: list
        The list of line plot object references
    """
    lat = np.array(lat)
    lon = np.array(lon)
    lon[(lon > 180.0)] = lon[(lon > 180.0)] - 360.0
    path = np.vstack((lon, lat))
    isplit = np.nonzero(np.abs(np.diff(path[0])) > 5.0)[0]
    segments = np.split(path, isplit + 1, axis=1)

    refs = []
    for segment in segments:
        x, y = m(segment[0], segment[1])
        refs.append(m.plot(x, y, ax=ax, color=color, alpha=alpha))
    return refs

def SwiftEarthPlot(trigtime,trigid,  prompt=True):
	if prompt == True:
		sclons, sclats, scalts = getPosFromTLE(trigtime)
	else:
	#get from sao file
		obsid = getobsid(trigtime.timestamp())
		saofile = getsaofile(obsid)
		times, sclats, sclons, scalts = orbitfromsao(saofile)[0:4]  #iknow the order doesnt look right, but it is

	fig, ax = plt.subplots(figsize=(20,10), dpi=100)
	lat_range = (-30.0, 30.0)
	lon_range = (-180.0, 180.0)
	map=Basemap(projection='merc', llcrnrlat=lat_range[0],
							urcrnrlat=lat_range[1], llcrnrlon=lon_range[0],
							urcrnrlon=lon_range[1], lat_ts=0, resolution='c',
							ax= ax)
	map.drawcoastlines()
	map.drawparallels(np.arange(-90., 91., 30.), labels=[1, 0, 0, 0],
						fontsize=12)
	map.drawmeridians(np.arange(-180., 181., 30.), labels=[0, 0, 0, 1],
						fontsize=12)

	#SC trajectory
	artist = earth_line(sclats,sclons, map, ax)

  #pos at trig time
	if prompt:
		triglon, triglat = sclons[1000],sclats[1000]
	else:
		unixtrig = trigtime.timestamp()
		mettrig = unixtime2sc(unixtrig)
		listvalue =  min(times, key=lambda x:abs(x-mettrig))
		index = times.index(listvalue)
		triglon, triglat = sclons[index], sclats[index]
		if triglon > 180.0:
			triglon = -(360-triglon)
	

	# SAA polygon
	lons,lats=list(zip(*BATSAA(triglon,triglat).points)) 
	x,y = map(lons,lats)
	xy = np.column_stack([x,y])
	poly = Polygon( xy, facecolor='darkred', alpha=0.4 )
	ax.add_patch(poly)

	# Plot SC position at trigtime with Icon
	# x_size, y_size = 0.8, 0.4
	# x0, y0 = map(sclons[1000] - x_size/2., sclats[1000] - y_size/2.)
	# x1, y1 = map(sclons[1000] + x_size/2., sclats[1000] + y_size/2.)
	# im = plt.imshow(plt.imread(f,0), extent=(x0, x1, y0, y1))

	f = urlopen("https://www.n2yo.com/inc/saticon.php?t=0&s=28485")
	im = OffsetImage(plt.imread(f,0), zoom=1)
	ab = AnnotationBbox(im, (map(triglon,triglat)), xycoords='data', frameon=False)
	map._check_ax().add_artist(ab)


	#McIlwain L
	lat_range = (-27.0, 27.0)
	lon_range = (-180.0, 180.0)
	artist = mcilwain_map(lat_range, lon_range, map, ax, alpha=0.5, saa_mask=poly)
	cb = plt.colorbar(artist, label='McIlwain L', ax=ax, shrink=0.6,
					pad=0.05,
					orientation='horizontal')
	cb.draw_all()

	if prompt:
		title = f'E Long, Lat={round(triglon,1), round(triglat,1)} FROM TLE'
		ax.set_title(title)
	else:
		title = f'E Long, Lat={round(triglon,1), round(triglat,1)} FINAL '
		ax.set_title(title)

	filename = f'EarthPlot.png'
	plt.savefig(filename)

	return filename