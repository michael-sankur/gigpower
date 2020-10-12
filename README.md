# LinDist3Flow

## Getting started with the powerflow python package
The python package that will house all the stable python module for power flow work is currently in `./fbs`.
(`powerflowpy` is a working title, and the module is under active development)
We are using [conda](https://docs.conda.io/projects/conda/en/latest/index.html) to manage packages and environments.  
If you would like to run or alter any code in `powerflowpy` here is how to get set up: 
1. If you don't have Conda, [Install Conda or Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Create the `pfp-env` virtual environment  
Open a terminal and navigate to your `LinDist3Flow` root.  
Enter the following at the command line:  
`conda env create -f environment.yml`
3. Activate the pfp-env environment.  
Run the following at the command line:  
`conda activate pfp-env`  
You should see your command line prompt prepended with `(pfp-env)`, like so:
`(pfp-env) YourComputer:LinDist3Flow username$`
4. To deactivate the environment, enter  
`conda deactivate`  
