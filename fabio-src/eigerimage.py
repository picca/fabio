#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Reader for Eiger based HDF5 files

License: GPLv2+

Authors:
........
* Jérôme Kieffer:
  European Synchrotron Radiation Facility;
  Grenoble (France)
"""
from __future__ import with_statement, print_function, absolute_import, division

__author__ = "Jerome Kieffer"
__contact__ = "Jerome.Kieffer@ESRF.eu"
__license__ = "GPLv2+"
__copyright__ = "European Synchrotron Radiation Facility, Grenoble, France"
__date__ = "03/02/2015"
__status__ = "development"

# get ready for python3
import os, logging
logger = logging.getLogger("eigerimage")
import numpy
from .fabioimage import FabioImage
from .third_party import six
from .nexus import h5py, Nexus

class EigerImage(FabioImage):

    def _readheader(self, infile):
        """ Read a Eiger image header """

        infile.seek(0)

        self.header = {}
        for name, nbytes, format in GE_HEADER_INFO:
            if format is None:
                self.header[ name ] = infile.read(nbytes)
            else:
                self.header[ name ] = struct.unpack(format,
                                                     infile.read(nbytes))[0]

    def read(self, fname, frame=None):
        """
        Read in header into self.header and
        the data   into self.data
        """
        if frame is None:
            frame = 0
        self.header = {}
        self.resetvals()

        infile = self._open(fname, "rb")
        self.sequencefilename = fname
        self._readheader(infile)
        self.nframes = self.header['NumberOfFrames']
        self._readframe(infile, frame)
        infile.close()
        return self

    def _makeframename(self):
        """ The thing to be printed for the user to represent a frame inside
        a file """
        self.filename = "%s$%04d" % (self.sequencefilename,
                                   self.currentframe)

    def _readframe(self, filepointer, img_num):
        """
        # Load only one image from the sequence
        #    Note: the first image in the sequence 0
        # raises an exception if you give an invalid image
        # otherwise fills in self.data
        """
        if(img_num > self.nframes or img_num < 0):
            raise Exception("Bad image number")
        imgstart = self.header['StandardHeaderSizeInBytes'] + \
                   self.header['UserHeaderSizeInBytes'] + \
                   img_num * self.header['NumberOfRowsInFrame'] * \
                   self.header['NumberOfColsInFrame'] * \
                   self.header['ImageDepthInBits'] // 8
        # whence = 0 means seek from start of file
        filepointer.seek(imgstart, 0)

        self.bpp = self.header['ImageDepthInBits'] // 8  # hopefully 2
        imglength = self.header['NumberOfRowsInFrame'] * \
                    self.header['NumberOfColsInFrame'] * self.bpp
        if self.bpp != 2:
            logging.warning("Using uint16 for GE but seems to be wrong")

        # Guessing it is always unsigned int?
        self.data = numpy.fromstring(filepointer.read(imglength), numpy.uint16)
        self.data.shape = (self.header['NumberOfRowsInFrame'],
                            self.header['NumberOfColsInFrame'])
        self.dim2 , self.dim1 = self.data.shape
        self.currentframe = int(img_num)
        self._makeframename()


    def write(self, fname, force_type=numpy.uint16):
        """ Not yet implemented"""
        raise Exception("Write is not implemented")

    def getframe(self, num):
        """
        Returns a frame as a new fabioimage object
        """
        if num < 0 or num > self.nframes:
            raise Exception("Requested frame number is out of range")
        # Do a deep copy of the header to make a new one
        newheader = {}
        for k in self.header.keys():
            newheader[k] = self.header[k]
        frame = GEimage(header=newheader)
        frame.nframes = self.nframes
        frame.sequencefilename = self.sequencefilename
        infile = frame._open(self.sequencefilename, "rb")
        frame._readframe(infile, num)
        infile.close()
        return frame

    def next(self):
        """
        Get the next image in a series as a fabio image
        """
        if self.currentframe < (self.nframes - 1) and self.nframes > 1:
            return self.getframe(self.currentframe + 1)
        else:
            newobj = GEimage()
            newobj.read(next_filename(
                self.sequencefilename))
            return newobj

    def previous(self):
        """
        Get the previous image in a series as a fabio image
        """
        if self.currentframe > 0:
            return self.getframe(self.currentframe - 1)
        else:
            newobj = GEimage()
            newobj.read(previous_filename(
                self.sequencefilename))
            return newobj

eigerimage = EigerImage