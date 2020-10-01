# powerflowpy
__About__: powerflowpy is a python module that parses dss files into Python data structures, implements a forward/backward sweep solver in python, and compares powerflow solution results with opendss.  

__Authors__: [@elaguerta](https://github.com/elaguerta) is the primary author of this module, with code review by [@msankur](https://github.com/msankur/) and [@toanngosy](https://github.com/toanngosy).

__History__: The main submodules `utils.py`, `network.py`, `solution.py` and `fbs.py` were written in August/September 2020.  

## Usage
Please see below for use cases and the module functions available to the user.  
Use Case| Module function | Notes
---|---|---

1. Read data from a radial network represented in a `.dss` file into Python data structures  
    a. see `powerflowpy.utils.init_from_dss()`
2. Run a forward/backward sweep (FBS) solver on the network
3. Output the results
4. Compare  

### Testing Suite
Tests were written using `pytest`. These tests live in `powerflowpy/tests` and can be run manually by invoking `pytest` on the powerflowpy directory:  
`$ pytest ./powerflowpy`

Below is a brief description of the test files and their coverage:  
1. `test_compare_matlab.py`  
Previous work on an fbs solver was written in MATLAB. The MATLAB solver relies on `network_mapper.m`. This 'network mapper' takes as input a text file representation of the network (in a different format from .dss) and outputs the relevant network parameters to a MATLAB struct.  
The analogous network mapper for `powerflowpy` is defined in `powerflowpy.utils.init_from_dss`, which maps a `.dss` file representing the network to a Python Network object. The `Network` class definition is in `powerflowpy.network`.   This test compares the results of the MATLAB mapper to the powerflowpy's network mapper on a 5 node test case. Please note that 

2. `test_pre_fbs.py`  
This test compares 





