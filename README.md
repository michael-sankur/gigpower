# gigpower

## Version 0.1.0

This Version is the initial version of gigpower, containing older code that was used as a a model. 
The Version establishes a baseline for test results described in `/tests/v1_test_report.csv`.


### Contents
- `/src/fbs` contains code for a forward-backward-sweep library completed at LBL GIG from 2020 - 2021, by Daniel Arnold, Michael Sankur, Sy-Toan Ngo, Elaine Laguerta, 
and others. Primary contact: [@elaguerta](https://github.com/elaguerta)
  - `/tests/test_solution_fbs.py` tests `/src/gigpower` against this library. 
- `/src/nr3_python` contains code for a Newton-Raphson library completed at LBL GIG from 2020-2021, by Daniel Arnold, Michael Sankur, Sy-Toan Ngo, Kathleen Chang,
and others. Primary contact: [@kathleenchang](https://github.com/kathleenchang)
  - `/tests/test_mapping_nr3.py` tests `/src/gigpower` against this library. 
  - This library also contains code for volt-var controllers. Future implementers of volt-var controllers can refer to this library. 
  Initial tests for volt-var-control are located in `/tests/test_vvc_nr3.py`
  
### Development
To work with this version:
1. Create a python virtual environment, and activate it. 
```
python3 -m venv venv
source bin/venv/activate
```
2. Navigate to `./gigpower/src`. Install requirements.
```
cd ./gigpower/src
pip install -r requirements.txt
pip install -r requirements-dev.txt
```
3. Navigate to `./gigpower/src`. Install modules locally in editable mode
```
pip install -e .
```
4. To run tests, run pytest as a module from `./gigpower`L
```
python3 -m pytest tests
```

Contact [Daniel Arnold](dbarnold@lbl.gov) for any questions about this version. 

