#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Project: Fable Input Output
#             https://github.com/silx-kit/fabio
#
#    Copyright (C) European Synchrotron Radiation Facility, Grenoble, France
#
#    Principal author:       Jérôme Kieffer (Jerome.Kieffer@ESRF.eu)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

"""Test failing files
"""

from __future__ import print_function, with_statement, division, absolute_import

import unittest
import os
import io
import fabio
import tempfile
import shutil


class TestFailingFiles(unittest.TestCase):
    """Test failing files"""

    @classmethod
    def setUpClass(cls):
        cls.tmp_directory = tempfile.mkdtemp()
        cls.createResources(cls.tmp_directory)

    @classmethod
    def createResources(cls, directory):

        cls.txt_filename = os.path.join(directory, "test.txt")
        f = io.open(cls.txt_filename, "w+t")
        f.write(u"Kikoo")
        f.close()

        cls.bad_edf_filename = os.path.join(directory, "bad_edf.edf")
        f = io.open(cls.bad_edf_filename, "w+b")
        f.write(b"\r{")
        f.write(b"\x00\xFF\x99" * 10)
        f.close()

        cls.bad_edf2_filename = os.path.join(directory, "bad_edf2.edf")
        f = io.open(cls.bad_edf2_filename, "w+b")
        f.write(b"\n{\n\n}\n")
        f.write(b"\xFF\x00\x99" * 10)
        f.close()

        cls.bad_msk_filename = os.path.join(directory, "bad_msk.msk")
        f = io.open(cls.bad_msk_filename, "w+b")
        f.write(b'M\x00\x00\x00A\x00\x00\x00S\x00\x00\x00K\x00\x00\x00')
        f.write(b"\x00\xFF\x99" * 10)
        f.close()

        cls.bad_dm3_filename = os.path.join(directory, "bad_dm3.dm3")
        f = io.open(cls.bad_dm3_filename, "w+b")
        f.write(b'\x00\x00\x00\x03')
        f.write(b"\x00\xFF\x99" * 10)
        f.close()

        cls.bad_npy_filename = os.path.join(directory, "bad_numpy.npy")
        f = io.open(cls.bad_npy_filename, "w+b")
        f.write(b"\x93NUMPY")
        f.write(b"\x00\xFF\x99" * 10)
        f.close()

        cls.missing_filename = os.path.join(directory, "test.missing")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp_directory)

    def test_missing_file(self):
        self.assertRaises(IOError, fabio.open, self.missing_filename)

    def test_wrong_format(self):
        self.assertRaises(IOError, fabio.open, self.txt_filename)

    def test_wrong_edf(self):
        self.assertRaises(IOError, fabio.open, self.bad_edf_filename)

    def test_wrong_edf2(self):
        self.assertRaises(IOError, fabio.open, self.bad_edf_filename)

    def test_wrong_msk(self):
        self.assertRaises(ValueError, fabio.open, self.bad_msk_filename)

    def test_wrong_dm3(self):
        self.assertRaises(ValueError, fabio.open, self.bad_dm3_filename)

    def test_wrong_numpy(self):
        self.assertRaises(ValueError, fabio.open, self.bad_npy_filename)


def suite():
    loadTests = unittest.defaultTestLoader.loadTestsFromTestCase
    testsuite = unittest.TestSuite()
    testsuite.addTest(loadTests(TestFailingFiles))
    return testsuite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
