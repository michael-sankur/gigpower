import sys
import powerflowpy
from powerflowpy.utils import init_from_dss

if __name__ == '__main__':
    # Runs from command line with path to dss file as an argument
    dss_file = sys.argv[1]
    network = init_from_dss(dss_file)
    sys.stdout.write(str(network))
