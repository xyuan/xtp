import sys
sys.path.append("/home/pranav/dev/votca/src/scripts/")

from XtpPy.OrbReader import OrbReader2

of = OrbReader2("../build/src/tests/xtp_testing.hdf5")
of.Print()
