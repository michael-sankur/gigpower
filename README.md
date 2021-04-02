# LinDist3Flow

LinDist3Flow is a python package used to solve power flow for radial and meshed networks. The package relies on OpenDSSPy to parse text files in *.dss format into custom package objects. The package uses Newton-Raphson and Forward/Backward Sweep algorithms to solve power flow. The package allows users to solve iteratively over time-varying models, with volt-var control and regulator controls. The package provides methods for users to access and calculate common solution parameters through Pandas DataFrames.

This project was motivated by looking for alternatives to OpenDSS to conduct certain high speed simulations with a custom Python solver, in hopes that it might be faster than OpenDSS. We did not accomplish the goal of developing a faster solver. However, this package offers the option for the programmer to solve powerflow with either NR3 or FBS. It also offers a basic, lean representation of feeders with a lean Python-based object model that may be preferred to the OpenDSS representation for some applications. Please leave a comment on here [insert issue number] about how you are using this package!

## Installing

TODO:
LinDist3Flow is available on PyPI:

`
$ pip install [NAME]
`  

## Simple Example
```python
def new_nr3_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionNR3(fp)
    # solution.maxiter = 1
    solution.solve()
    return solution

nr3_vals = new_nr3_solution.get_data_frame(solution_param)
print(nr3_vals)

def new_fbs_solution(dss_file):
    fp = str(Path(DSS_FILE_DIR, dss_file))
    solution = SolutionFBS(fp)
    # solution.maxiter = 1
    solution.solve()
    return solution

fbs_vals = new_fbs_solution.get_data_frame(solution_param)
print(fbs_vals)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## Disclosure

This project is related to [Topology Reconfiguration and distributed Volt Var Control (TR-dVVC) v1](https://www.osti.gov/biblio/1567734).
