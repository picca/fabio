

"""
Test cases for the fabioimage clas
"""

from fabio import fabioimage
import unittest, os, sys, Numeric, RandomArray

class test50000(unittest.TestCase):
    """ test with 50000 everywhere"""
    def setUp(self):
        """make the image"""
        dat = Numeric.ones((1024, 1024), Numeric.UInt16)
        dat = (dat * 50000).astype(Numeric.UInt16)
        assert dat.typecode() == Numeric.ones((1), Numeric.UInt16).typecode()
        hed = {"Title":"50000 everywhere"}
        self.obj = fabioimage(dat, hed)
      
    def testgetmax(self):
        """check max"""
        self.assertEqual( self.obj.getmax(), 50000)

    def testgetmin(self):
        """check min"""
        self.assertEqual( self.obj.getmin(), 50000)
    
    def testgetmean(self):
        """check mean"""
        self.assertEqual( self.obj.getmean(), 50000)
        
    def getstddev(self):
        """check stddev"""
        self.assertEqual( self.obj.getstddev(), 0)
        
class testslices(unittest.TestCase):
    """check slicing"""
    def setUp(self):
        """make test data"""
        dat2 = Numeric.zeros((1024, 1024), Numeric.UInt16, savespace = 1 )
        hed = {"Title":"zeros and 100"}
        self.cord = [ 256, 256, 790, 768 ]
        self.obj = fabioimage(dat2, hed)
        self.slic = slic = self.obj.make_slice(self.cord)
        # Note - d2 is modified *after* fabioimage is made
        dat2[slic] = dat2[slic] + 100
        assert self.obj.maxval is None
        assert self.obj.minval is None
        self.npix = (slic[0].stop - slic[0].start) * \
            (slic[1].stop - slic[1].start)
        
    def testgetmax(self):
        """check max"""
        self.assertEqual( self.obj.getmax(), 100)

    def testgetmin(self):
        """check min"""
        self.assertEqual( self.obj.getmin(), 0)
        
    def testintegratearea(self):
        """ check integrations"""
        self.obj.resetvals()
        area1 = self.obj.integrate_area(self.cord) 
        self.obj.resetvals()
        area2 = self.obj.integrate_area(self.slic)
        self.assertEqual(area1, area2)
        self.assertEqual(area1, self.npix*100)
        
    
class testopen(unittest.TestCase):
    """check opening compressed files"""
            
    def setUp(self):
        """ create test files"""
        open("testfile","wb").write("{ hello }")
        os.system("gzip testfile")
        open("testfile","wb").write("{ hello }")
        os.system("bzip2 testfile")
        open("testfile","wb").write("{ hello }")
        self.obj = fabioimage()

    def tearDown(self):
        """clean up"""
        for name in ["testfile", "testfile.gz", "testfile.bz2"]:
            if os.path.exists(name):
                os.remove(name)

    def testFlat(self):
        """ no compression"""
        res = self.obj._open("testfile").read()
        self.assertEqual( res , "{ hello }" ) 

    def testgz(self):
        """ gzipped """
        res = self.obj._open("testfile.gz").read()
        self.assertEqual( res , "{ hello }" ) 
    
    def testbz2(self):    
        """ bzipped"""
        res = self.obj._open("testfile.bz2").read()
        self.assertEqual( res , "{ hello }" ) 


NAMES = { Numeric.UInt8 : "Numeric.UInt8",
          Numeric.Int8  : "Numeric.Int8" ,  
          Numeric.UInt16: "Numeric.UInt16",  
          Numeric.Int16 :  "Numeric.Int16" ,  
          Numeric.UInt32:  "Numeric.UInt32" , 
          Numeric.Int32 :  "Numeric.Int32"   ,
          Numeric.Float32: "Numeric.Float32" ,
          Numeric.Float64:"Numeric.Float64"}


class testPILimage(unittest.TestCase):
    """ check PIL creation"""
    def setUp(self):
        """ list of working numeric types"""
        self.okformats = [Numeric.UInt8,
                          Numeric.Int8,
                          Numeric.UInt16,
                          Numeric.Int16,
                          Numeric.UInt32,
                          Numeric.Int32,
                          Numeric.Float32]

    def mkdata(self, shape, typ):
        """ generate [01] testdata """
        return (RandomArray.random(shape)).astype(typ)

        
    def testpil(self):
        """ check all formats with random data"""
        for typ in self.okformats:
            name = NAMES[typ]
            for shape in [(10, 20), (431, 1325)]:
                testdata = self.mkdata(shape, typ)
                img = fabioimage(testdata, {"title":"Random data"})
                pim = img.toPIL16()
                for i in [ 0, 5, 6, shape[1]-1 ]:
                    for j in [0, 5, 7, shape[0]-1 ]:
                        errstr = name + " %d %d %f %f t=%s" % (
                            i, j, testdata[j, i], pim.getpixel((i, j)), typ)
                        
                        er1 = img.data[j, i] - pim.getpixel((i, j))
                        er2 = img.data[j, i] + pim.getpixel((i, j))
                        
                        # difference as % error in case of rounding
                        if er2 != 0.:
                            err = er1/er2
                        else:
                            err = er1

                        self.assertAlmostEquals( err,
                                                 0,
                                                 6, 
                                                 errstr)


class testPILimage2(testPILimage):
    """ check with different numbers"""
    def mkdata(self, shape, typ):
        """ positive and big"""
        return (RandomArray.random(shape) * sys.maxint / 10).astype(typ)

class testPILimage3(testPILimage):
    """ check with different numbers"""
    def mkdata(self, shape, typ):
        """ positive, negative and big"""
        return ((RandomArray.random(shape) - 0.5) * sys.maxint / 10).astype(typ)
        
    
if __name__ == "__main__":
    unittest.main()
