# coding: utf-8
#
#    Project: X-ray image reader
#             https://github.com/silx-kit/fabio
#
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
#
#    Principal author:       Jérôme Kieffer (Jerome.Kieffer@ESRF.eu)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE

"""

Authors: Henning O. Sorensen & Erik Knudsen
         Center for Fundamental Research: Metal Structures in Four Dimensions
         Risoe National Laboratory
         Frederiksborgvej 399
         DK-4000 Roskilde
         email:henning.sorensen@risoe.dk

mods for fabio by JPW
modification for HDF5 by Jérôme Kieffer

"""

from __future__ import with_statement, print_function, absolute_import

import os.path
import logging
logger = logging.getLogger(__name__)
from . import fabioutils
from .fabioutils import FilenameObject, BytesIO, NotGoodReader
from .fabioimage import FabioImage

# Make sure to load all formats
from . import fabioformats  # noqa


MAGIC_NUMBERS = [
    # "\42\5a" : 'bzipped'
    # "\1f\8b" : 'gzipped'
    (b"FORMAT :100", 'bruker100'),
    (b"FORMAT :        86", 'bruker'),
    (b"\x4d\x4d\x00\x2a", 'tif'),
    # The marCCD and Pilatus formats are both standard tif with a header
    # hopefully these byte patterns are unique for the formats
    # If not the image will be read, but the is missing
    (b"\x49\x49\x2a\x00\x08\x00", 'marccd/tif'),
    (b"\x49\x49\x2a\x00\x82\x00", 'pilatus'),
    (b"\x49\x49\x2a\x00", 'tif'),
    # d*TREK must come before edf
    (b"{\nHEA", 'dtrek'),
    (b"{", 'edf'),
    (b"\r{", 'edf'),
    (b"\n{", 'edf'),
    (b"ADEPT", 'GE'),
    (b"OD", 'OXD'),
    (b"IM", 'HiPiC'),
    (b'\x2d\x04', 'mar345'),
    (b'\xd2\x04', 'mar345'),
    (b'\x04\x2d', 'mar345'),  # some machines may need byteswapping
    (b'\x04\xd2', 'mar345'),
    # hint : MASK in 32 bit
    (b'M\x00\x00\x00A\x00\x00\x00S\x00\x00\x00K\x00\x00\x00', 'fit2dmask'),
    (b'\x00\x00\x00\x03', 'dm3'),
    (b"No", "kcd"),
    (b"<", "xsd"),
    (b"\n\xb8\x03\x00", 'pixi'),
    (b"\x89\x48\x44\x46\x0d\x0a\x1a\x0a", "eiger/hdf5"),
    (b"R-AXIS", 'raxis'),
    (b"\x93NUMPY", 'numpy'),
    (b"\\$FFF_START", 'fit2d'),
    # Raw JPEG
    (b"\xFF\xD8\xFF\xDB", "jpeg"),
    # JFIF format
    (b"\xFF\xD8\xFF\xE0", "jpeg"),
    # Exif format
    (b"\xFF\xD8\xFF\xE1", "jpeg"),
    # JPEG 2000 (from RFC 3745)
    (b"\x00\x00\x00\x0C\x6A\x50\x20\x20\x0D\x0A\x87\x0A", "jpeg2k"),
]


def do_magic(byts, filename):
    """ Try to interpret the bytes starting the file as a magic number """
    for magic, format_type in MAGIC_NUMBERS:
        if byts.startswith(magic):
            if "/" in format_type:
                if format_type == "eiger/hdf5":
                    if "::" in filename:
                        return ["hdf5"]
                    else:
                        return ["eiger", "soleil"]
                elif format_type == "marccd/tif":
                    if "mccd" in filename.split("."):
                        return ["marccd"]
                    else:
                        return ["tif"]
            return [format_type]
    return None


def openimage_str(filename, frame):
    logger.debug("Attempting to open %s" % (filename))
    objs = _openimage(filename)
    for obj in objs:
        try:
            return obj.read(filename, frame)
        except NotGoodReader as ex:
            pass
    raise NotGoodReader("Did not found a good reader for {}".format(filename))

def openimage_PathTypes(filename, frame):
    if not isinstance(filename, fabioutils.StringTypes):
        filename = str(filename)
    openimage_str(filename, frame)

def openimage_FilenameObject(filename, frame):
    try:
        return openimage_str(filename.tostring(), frame)
    except Exception as ex:
        logger.debug("Exception %s, trying name %s" % (ex, filename.stem))
        return openimage_str(filename.stem, filename.num)

def openimage(filename, frame=None):
    """Open an image.

    It returns a FabioImage-class instance which can be used as a context manager to close the file
    at the termination.

    .. code-block:: python

        with fabio.open("image.edf") as i:
            print(i.nframes)
            print(i.data)

    :param Union[str,FilenameObject] filename: A filename or a filename
        iterator.
    :param Union[int,None] frame: A specific frame inside this file.
    :rtype: FabioImage
    """
    if isinstance(filename, str):
        return openimage_str(filename, frame)
    elif isinstance(filename, unicode):
        return openimage_str(str(filename), frame)
    elif isinstance(filename, fabioutils.PathTypes):
        return openimage_PathTypes(filename, frame)
    elif isinstance(filename, FilenameObject):
        return openimage_FilenameObject(filename, frame)

def openheader(filename):
    """ return only the header"""
    if isinstance(filename, fabioutils.PathTypes):
        if not isinstance(filename, fabioutils.StringTypes):
            filename = str(filename)

    objs = _openimage(filename)
    for obj in objs:
        try:
            obj.readheader(obj.filename)
            return obj
        except NotGoodReader:
            pass

def _openimage(filename):
    """
    determine which format for a filename
    and return appropriate class which can be used for opening the image

    :param filename: can be an url like:

    hdf5:///example.h5?entry/instrument/detector/data/data#slice=[:,:,5]

    """
    if hasattr(filename, "seek") and hasattr(filename, "read"):
        # Looks to be a file containing filenames
        if not isinstance(filename, BytesIO):
            filename.seek(0)
            actual_filename = BytesIO(filename.read())
    else:
        if os.path.exists(filename):
            # Already a valid filename
            actual_filename = filename
        elif "::" in filename:
            actual_filename = filename.split("::")[0]
        else:
            actual_filename = filename

    try:
        imo = FabioImage()
        with imo._open(actual_filename) as f:
            magic_bytes = f.read(18)
    except IOError:
        logger.debug("Backtrace", exc_info=True)
        raise
    else:
        imo = None

    filetypes = do_magic(magic_bytes, filename)
    if filetypes is None:
        logger.debug("Backtrace", exc_info=True)
        try:
            file_obj = FilenameObject(filename=filename)
            if file_obj is None:
                raise Exception("Unable to deconstruct filename")
            if (file_obj.format is not None) and\
               len(file_obj.format) != 1 and \
               isinstance(file_obj.format, list):
                # one of OXD/ADSC - should have got in previous
                raise Exception("openimage failed on magic bytes & name guess")
            if file_obj.format is not None:
                filetypes = [file_obj.format]
        except Exception:
            logger.debug("Backtrace", exc_info=True)
            raise IOError("Fabio could not identify " + filename)

    if filetypes is None:
        raise IOError("Fabio could not identify " + filename)

    objs = []
    for filetype in filetypes:
        klass_name = "".join(filetype) + 'image'
        try:
            obj = fabioformats.factory(klass_name)
            objs.append(obj)
        except (RuntimeError, Exception):
            pass

    if len(objs) == 0:
        logger.debug("Backtrace", exc_info=True)
        raise IOError("Filename %s can't be read as format %s" % (filename, filetypes))

    for obj in objs:
        obj.filename = filename

    # skip the read for read header
    return objs

def open_series(filenames=None, first_filename=None,
                single_frame=None, fixed_frames=None, fixed_frame_number=None):
    """
    Create an object to iterate frames through a file series.

    This function is a wrapper over :class:`~file_series.FileSeries` to facilitate
    simple uses of file series iterations.

    :param Union[Generator,Iterator,List] filenames: Ordered list of filenames
        to process as a file series. It also can be a generator, and
        iterator, or :class:`~fabio.file_series.filename_series` or
        :class:`~fabio.file_series.file_series` objects.
    :param str first_filename: If provided iterate filenames from this filename
        and try to consecutivelly open next files. If this argument is specified
        the `filenames` have to unspecified. Internally it uses
        :class:`~fabio.file_series.filename_series` to iterate the filenames.
    :param Union[Bool,None] single_frame: If True, all files are supposed to
        contain only one frame.
    :param Union[Bool,None] fixed_frames: If True, all files are supposed to
        contain the same amount of frames (this fixed amount will be reached
        from the first file of the serie).
    :param Union[Integer,None] fixed_frame_number: If set, all files are
        supposed to contain the same amount of frames (sepecified by this
        argument)
    :rtype: :class:`~file_series.FileSeries`
    """
    # Here to avoid recursive import
    from . import file_series

    if filenames is not None and first_filename is not None:
        raise ValueError("'filenames' and 'first_filename' are mutual exclusive")

    if first_filename is not None:
        filenames = file_series.filename_series(filename=first_filename)

    return file_series.FileSeries(filenames=filenames,
                                  single_frame=single_frame,
                                  fixed_frames=fixed_frames,
                                  fixed_frame_number=fixed_frame_number)
