import unittest
from libradtranpy.libsimulateVisible import Dict_Of_sitesAltitudes,Dict_Of_sitesPressures
from libradtranpy.libsimulateVisible import UVspec3,ApplyAerosols,ProcessSimulation,ProcessSimulationaer
from libradtranpy.libsimulateVisible import CleanSimDir
import numpy as np 
import os

# python -m unittest tests/libradtranpy/test_libsimulateVisible.py

class libsimulateVisible(unittest.TestCase):
    """A test case for the libsimulateVisible package."""


    def test_libradtranenvironment(self):
        """
        test if LIBRADTRANDIR environnement variable is defined
        """
        var = 'LIBRADTRANDIR'
        self.assertTrue(var in os.environ)

    def test_libradtranuvspecexecpath(self):
        """
        test if libRadtran uvspec executable exists
        """
        libradtranpath = os.getenv('LIBRADTRANDIR')+ '/'
        libradtranbinpath = libradtranpath + "/bin"
        libradtranexecutablepath = os.path.join(libradtranbinpath,'uvspec')
        self.assertTrue(os.path.exists(libradtranexecutablepath))

    def test_libradtrandatapath(self):
        """
        test if libradtran data-path exists
        """
        libradtranpath = os.getenv('LIBRADTRANDIR')+ '/'
        libradtrandatapath = libradtranpath + "/share/libRadtran/data"
        self.assertTrue(os.path.exists(libradtrandatapath))

    def test_ProcessSimulation_outputfilename(self):
        """
        test if libradtran output simulation exists 
        """
        CleanSimDir()
        outputdir,outputfilename = ProcessSimulation(1.2,4.0,300.,0.)
        fullnameoutputfilename=os.path.join(outputdir,outputfilename)
        self.assertTrue(os.path.isfile(fullnameoutputfilename))

    def test_ProcessSimulation_outputdataexists(self):
        """
        test if libradtran simulation output data is readable
        """
        CleanSimDir()
        outputdir,outputfilename = ProcessSimulation(1.2,4.0,300.,0.)
        fullnameoutputfilename=os.path.join(outputdir,outputfilename)
        data = np.loadtxt(fullnameoutputfilename)
        self.assertTrue(np.any(data))


    def test_ProcessSimulation_outputdatahascorrectshape(self):
        """
        test if libradtran simulation output data has correct shape
        """
        CleanSimDir()
        outputdir,outputfilename = ProcessSimulation(1.2,4.0,300.,0.)
        fullnameoutputfilename=os.path.join(outputdir,outputfilename)
        data = np.loadtxt(fullnameoutputfilename)
        self.assertEqual(data.shape[1], 2)
        self.assertEqual(data.shape[0], 951)
    




if __name__ == "__main__":
    unittest.main()