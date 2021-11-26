#From Adam Goldstein

import os
import astropy.io.fits as fits
from astropy import wcs
import numpy as np
import healpy as hp
from scipy.interpolate import interp1d
from scipy.spatial.transform import Rotation
from collections import OrderedDict
from matplotlib.pyplot import contour as Contour
from copy import deepcopy

from gbm.file import GbmFile
from gbm.data import PosHist
from gbm.data.data import DataFile
from gbm import coords

class SwiftFile(GbmFile):
    @classmethod
    def from_path(cls, path):
        dir, base = os.path.split(path)
        #dtype, _, det, tnum,  = base.split('_')
        base = base.split('.')[0]
        tnum = int(base[2:13])
        dtype = base[13:]
        
        obj = cls.create(directory=dir, data_type=dtype, 
                         detector='BAT', uid=tnum)
        return obj


class SwiftDataFile(DataFile):
    def _file_properties(self, filename, ignore_file_check=False):
        super()._file_properties(filename, ignore_file_check=False)
        self._filename_obj = SwiftFile.from_path(filename)


class SwiftPosHist(SwiftDataFile, PosHist):

    def __init__(self):
        PosHist.__init__(self)

    @property
    def time_range(self):
        return (self._times[0], self._times[-1])

    @classmethod
    def open(cls, filename):
        """Open and read a position history file
        
        Args:
            filename (str): The filename of the FITS file
        
        Returns:        
            :class:`PosHist`: The PosHist object
        """
        obj = cls()
        obj._file_properties(filename)
        # open FITS file
        with fits.open(filename) as hdulist:
            for hdu in hdulist:
                obj._headers.update({hdu.name: hdu.header})
            data = hdulist['PREFILTER'].data

        times = data['TIME']
        obj._times = times
        obj._data = data

        # set the interpolators
        obj._set_interpolators()

        return obj

    def _set_interpolators(self):
        times = self._times

        data = self._data
        # Earth inertial coordinates interpolator
        eic = data['POSITION'].T*1000.0 # need to convert to meters
        self._eic_interp = interp1d(times, eic)

        # quaternions interpolator
        quat = data['QUATERNION'].T
        self._quat_interp = interp1d(times, quat)

        # Orbital position interpolator
        # mark TODO: poshist uses the "simple" version of lat/lon calc
        self._lat_interp = interp1d(times, data['SAT_LAT'])
        self._lon_interp = interp1d(times, data['SAT_LON'])
        alt = self._altitude_from_scpos(eic)
        self._alt_interp = interp1d(times, alt)

        # Earth radius and geocenter interpolators
        self._earth_radius_interp = interp1d(times, self._geo_half_angle(alt))
        self._geocenter_interp = interp1d(times,
                                          coords.geocenter_in_radec(eic))


class SwiftBatPartialCoding():
    _file = 'pcode_default.img'
    nside = 128
    def __init__(self):
        hdulist = fits.open(self._file, memmap=False)
        w = wcs.WCS(hdulist[0].header)
        data = hdulist[0].data
        
        num_y, num_x = w.array_shape
        x = np.arange(num_x)
        y = np.arange(num_y)
        x, y = np.meshgrid(x, y)
        ra, dec = w.wcs_pix2world(x, y, 1)
        ra += 360.0
        
        npix = hp.nside2npix(self.nside)
        pix = hp.ang2pix(self.nside, ra, dec, lonlat=True)
        self._hpx = np.zeros(npix)
        self._hpx[pix] = data.reshape(pix.shape)

    def partial_coding_path(self, frac, numpts_ra=360, numpts_dec=180):
        """Return the bounding path for a given partial coding fraction
        
        Args:
            frac (float): The partial coding fraction (valid range 0-1)
            numpts_ra (int, optional): The number of grid points along the RA 
                                       axis. Default is 360.
            numpts_dec (int, optional): The number of grid points along the Dec 
                                        axis. Default is 180.
        
        Returns:
            [(np.array, np.array), ...]: A list of RA, Dec points, where each \
                item in the list is a continuous closed path.
        """
        # create the grid and integrated probability array
        grid_pix, phi, theta = self._mesh_grid(numpts_ra, numpts_dec)
        frac_arr = self._hpx[grid_pix]
        ra = self._phi_to_ra(phi)
        dec = self._theta_to_dec(theta)

        # use matplotlib contour to produce a path object
        contour = Contour(ra, dec, frac_arr, [frac])

        # get the contour path, which is made up of segments
        paths = contour.collections[0].get_paths()

        # extract all the vertices
        pts = [path.vertices for path in paths]

        # unfortunately matplotlib will plot this, so we need to remove
        for c in contour.collections:
            c.remove()

        return pts
    
    def rotate(self, quaternion):
        # use scipy to convert between quaternion and euler angle, which is
        # what healpy uses
        eulers = Rotation.from_quat(quaternion).as_euler('ZYX')
        
        # rotate partial coding map according to euler angle
        rotator = hp.Rotator(rot=np.rad2deg(eulers))
        rot_hpx = rotator.rotate_map_pixel(self._hpx)
        # rotate it again because the reference frame for the map is in 
        # equatorial coordinates instead of spacecraft coordinates
        rotator = hp.Rotator(rot=[0.0, 0.0, 90.0])
        rot_hpx = rotator.rotate_map_pixel(rot_hpx)
        
        obj = deepcopy(self)
        obj._hpx = rot_hpx
        return obj
    
    def _mesh_grid(self, num_phi, num_theta):
        # create the mesh grid in phi and theta
        theta = np.linspace(np.pi, 0.0, num_theta)
        phi = np.linspace(0.0, 2 * np.pi, num_phi)
        phi_grid, theta_grid = np.meshgrid(phi, theta)
        grid_pix = hp.ang2pix(self.nside, theta_grid, phi_grid)
        return (grid_pix, phi, theta)

    @staticmethod
    def _phi_to_ra(phi):
        return np.rad2deg(phi)

    @staticmethod
    def _theta_to_dec(theta):
        return np.rad2deg(np.pi / 2.0 - theta)

        
        