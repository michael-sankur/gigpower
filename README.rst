gigpower
========

``gigpower`` solves power flow for radial and meshed networks. It supports
solving power flow with both Newton-Raphson and Forward Backward Sweep algorithms,
written completely in Python.

Features
--------
- Parse `*.dss` files into `Circuit` objects, which you can interact with directly in Python. 
- Solve power flow for radial or meshed networks using Newton-Raphson method.
- Solve power flow radial networks using Forward Backward Sweep.
- Compare solutions to OpenDSSDirect.py_.

.. _OpenDSSDirect.py: https://github.com/dss-extensions/OpenDSSDirect.py

Installation
------------

Install gigpower by running::

    > pip install gigpower


Solve Power Flow
________________

Solve with Newton-Raphson::

    from gigpower.SolutionNR3 import SolutionNR3

    my_solution_nr3 = SolutionNR3('my_feeder.dss')
    my_solution_nr3.solve()

Solve with Forward Backward Sweep::

    from gigpower.SolutionFBS import SolutionFBS

    my_solution_fbs = SolutionFBS('my_feeder.dss')
    my_solution_fbs.solve()

Solve with _OpenDSSDirect::

    from gigpower.SolutionDSS import SolutionDSS

    my_solution_dss = SolutionDSS('my_feeder.dss')
    my_solution_dss.solve()

Compare solutions::

    # compares solved voltages
    from gigpower.pretty_print import compare_data_frames

    fbs_V = my_solution_fbs.get_data_frame('V')
    dss_V = my_solution_dss.get_data_frame('V')
    compare_data_frames(fbs_V, dss_V)

See https://github.com/LBNL-ETA/gigpower/docs/examples for Jupyter notebooks with examples.

Contribute
----------

- Issue Tracker: https://github.com/LBNL-ETA/gigpower/issues
- Source Code: https://github.com/LBNL-ETA/gigpower

Development
-----------

Follow the steps below to develop this project.

1. Create a python virtual environment, and activate it.::

    python3 -m venv venv
    source bin/venv/activate

2. Navigate to `./gigpower/src`. Install requirements.::

    cd ./gigpower/src
    pip install -r requirements.txt
    pip install -r requirements-dev.txt

3. Navigate to `./gigpower/src`. Install modules locally in editable mode::

    pip install -e .

4. To run tests, run pytest as a module from `./gigpower`::

    python3 -m pytest tests


Support
-------
This project is maintained by the `Lawrence Berkeley National Lab Grid Integration Group <https://gridintegration.lbl.gov/>`_. 
For support, contact Daniel Arnold at: dbarnold@lbl.gov

License
-------

The project is licensed under the BSD license.
