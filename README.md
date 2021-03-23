# LinDist3Flow

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
