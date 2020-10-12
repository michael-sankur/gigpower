# powerflowpy
__About__: powerflowpy is a python module that parses dss files into Python data structures, implements a forward/backward sweep solver in python, and compares powerflow solution results with opendss.  

__Authors__: [@elaguerta](https://github.com/elaguerta) is the primary author of this module, with code review by [@msankur](https://github.com/msankur/) and [@toanngosy](https://github.com/toanngosy).

__History__: The main submodules `utils.py`, `network.py`, `solution.py` and `fbs.py` were written in August/September 2020.  

## Usage
Please see below for helpful module functions available to the user.  
Goal    | Module function   
---     |---                
Initialize a `Network` object from a `.dss` file. | `utils.init_from_dss()`
Run forward/backward sweep (FBS) to solve powerflow for a `Network.` | `fbs.fbs()`  
Get results from FBS solution | `fbs.get_solution()`
Change the kW and kvar of a Load within a Network. | `network.Load.set_kW()`, and `network.Load.set_kvar()`. See also [this gist](https://gist.github.com/elaguerta/acab7589b770fe32359e336a6759a2b9).
Alter a `Network` object, or other objects that compose the `Network`. | See class definitions in `network.py`  
Compare opendss' solution to the FBS solution | `tests/test_compare_matlab.py`
Conduct a timing experiment on opendss vs. FBS | `main.py`

## Testing Suite
Tests were written using `pytest`. These tests live in `fbs/tests` and can be run manually by invoking `pytest` on the powerflowpy directory:  
`$ pytest ./fbs`

Below is a brief description of the test files and their coverage:
<hr>

### `test_compare_matlab.py`  
Previous work on an fbs solver was written in MATLAB. The MATLAB solver relies on `network_mapper.m`. This 'network mapper' takes as input a text file representation of the network (in a different format from .dss) and outputs the relevant network parameters to a MATLAB struct. 

The analogous network mapper for `powerflowpy` is defined in `powerflowpy.utils.init_from_dss`, which maps a `.dss` file representing the network to a Python Network object. The `Network` class definition is in `powerflowpy.network`.

    
This test compares the results of the MATLAB mapper to the powerflowpy's network mapper on a 5 node test case. You can find this test case in `fbs/tests/05n3ph_unbal`. The test was constructed by the following process: 

1. Run `network_mapper.m` on the txt file representation of the network.
    - See `compare_lindist3flow_05node_threephase_unbalanced_oscillation_03.txt` 
2. Save the resulting MATLAB structs as .mat files.
    - See all the `*.mat` files in `fbs/tests/05n3ph_unbal`
3. Run `utils.init_from_dss()` on the corresponding .dss file of the network. 
     - See fbs/tests/05n3ph_unbal/`compare_opendss_05node_threephase_unbalanced_oscillation_03.dss`
4. Compare the MATLAB structs with the resulting Network object obtained in 3.

<hr>

### `test_pre_fbs.py`  

This test file has minimal coverage and was written as a smoke-test. It has two main goals:
1. Compare the Python `Network` object with the opendss `Circuit` object. We check that the Number of nodes and number of Lines match. 
2. Checks that `fbs.topo_sort()` has correctly sorted the network for the FBS algorithm.

<hr>

### `test_compare_opendss.py`

This test file compares the solution values obtained by opendss vs. FBS. It checks that the values are comparable within a specified tolerance.





