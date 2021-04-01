# LinDist3Flow

LinDist3Flow is a python package used to solve power flow for radial and meshed networks. The package relies on OpenDSSPy to parse text files in *.dss format into custom package objects. The package uses Newton-Raphson and Forward/Backward Sweep algorithms to solve power flow. Calculations rely on the CPython standard library and on the Numpy package. The package allows users to solve iteratively over time-varying models, with volt-var control and regulator controls. The package provides methods for users to access and calculate common solution parameters through Pandas Data

## Usage

TODO:
`$ pip install [NAME]`  

## Getting started
From the project root (the parent of `src/`): 

1. Get thee a virtual env, and activate it:  
`$ python -m venv .`  
`$ source bin/activate`  

2. Update pip and install dependencies with pip:  
`$ pip install --upgrade pip`  
`$ pip install -r requirements.txt`

3. Install python packages in the repo, with editable flag:  
`$ pip install -e .`  
Note: this step is temporary, for developers working on the package.


## Running tests
1. Install pytest:  
`$ pip install pytest`

2. Run pytest on the `tests` folder:  
`$ pytest tests`

3. To capturing all prints to stdout, run with `-s`:  
`$ pytest tests -s`


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
