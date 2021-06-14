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

gigpower Copyright (c) 2021, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

(1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

(2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

(3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Dept. of Energy, nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

You are under no obligation whatsoever to provide any bug fixes, patches, or upgrades to the features, functionality or performance of the source code ("Enhancements") to anyone; however, if you choose to make your Enhancements available either publicly, or directly to Lawrence Berkeley National Laboratory, without imposing a separate written license agreement for such Enhancements, then you hereby grant the following license: a non-exclusive, royalty-free perpetual license to install, use, modify, prepare derivative works, incorporate into other computer software, distribute, and sublicense such enhancements or derivative works thereof, in binary and source code form.
