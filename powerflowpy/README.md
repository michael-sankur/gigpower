# powerflowpy
__About__: powerflowpy is a python module that parses dss files into Python data structures, implements a forward/backward sweep solver in python, and compares powerflow solution results with opendss.  

__Authors__: [@elaguerta](https://github.com/elaguerta) is the primary author of this module, with code review by [@msankur](https://github.com/msankur/) and [@toanngosy](https://github.com/toanngosy).

__History__: The main modules `utils.py`, `network.py`, `solution.py` and `fbs.py` were written in August/September 2020.  

## Testing Suite
Tests were written using `pytest`. These tests live in `powerflowpy/tests` and can be run manually by invoking `pytest` on the powerflowpy directory:  
`$ pytest ./powerflowpy`

Below is a brief description of the test files and their coverage:  
1. `test_compare_matlab.py`  
Previous work on an fbs solver was written in MATLAB. The MATLAB workflow maps a text file representation of the network to a MATLAB struct. In contrast, powerflowpy's network mapper, defined in `powerflowpy.utils.init_from_dss` maps a dss file representing the network to a Python Network object, defined in `powerflowpy.network`.   This test compares the results of the MATLAB mapper to the powerflowpy's network mapper on a 5 node test case. Please note that 

2. `test_pre_fbs.py`  
This test compares 





