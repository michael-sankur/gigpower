# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: 23 March 2020
# Compare solutions between opendss and SolutiionFBS
# based on tests/test_dss_compare_fbs.py

from circuit_mapper.solution_dss import SolutionDSS
from circuit_mapper.solution_fbs import SolutionFBS
from circuit_mapper.pretty_print import compare_data_frames

import pandas as pd
import numpy as np

import sys
import os
from pathlib import Path
import requests
import pytest

DSS_FILE_REPO = 'https://raw.githubusercontent.com/lbnl-cybersecurity/python-powerflow/master/'
OUT_DIR = Path('./tests/test_compare_results')

# dictionaries with FILE_STEM: API TOKEN
# API token is used to download source files from raw.githubuser.com
DSS_CKT_FILES = {
                'IEEE_13_Node/IEEE_13_Bus_allwye.dss': 'ADCUMKIJTZZUEN5R3OXNX3DAMQCHY',
                # 'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss': 'ADCUMKM27O7EDMUSGZZXHEDAIF2EQ',
                # 'IEEE_34_feeder_UB/IEEE_34_Bus_allwye.dss': 'ADCUMKPSDEI4H6FJ64E3WP3AIF2V6',
                # 'IEEE_34_feeder_UB/IEEE_34_Bus_allwye_noxfm_noreg.dss': 'ADCUMKKKHOWAGZAFIYF74K3AIF2ZG',
                # 'IEEE_37_feeder_UB/IEEE_37_Bus_allwye.dss': 'ADCUMKMRXJ6XFKXF2MPS4OTAIF3BE',
                # 'IEEE_37_feeder_UB/IEEE_37_Bus_allwye_noxfm_noreg.dss': 'ADCUMKP36GFYRQDGQKDLZK3AIF3FA'
                }
DSS_CONFIG_FILES = {
                    'IEEE_13_Node/IEEELineCodes.dss': 'ADCUMKNWINSZLV323HRWDSTAMQCPQ',
                    'IEEE_13_Node/IEEE13Node_BusXY.csv': 'ADCUMKLAYMADVC4CEN5PFP3AMQCVK',
                    # 'IEEE_34_feeder_UB/IEEE34_BusXY.csv': 'ADCUMKPOJ2UU4GZUNGDS6B3AIF2NC',
                    # 'IEEE_34_feeder_UB/IEEELineCodes.dss': 'ADCUMKMONSF4CEDEPYC2RUTAIF2RG',
                    # 'IEEE_37_feeder_UB/IEEELineCodes.dss': 'ADCUMKNPPPEQ3RYI6MA7CK3AIF26M'
                    }


def test_all():
    """
    Compare the python FBS solution to the opendss solution.
    Writes output to OUTPUT FOLDER.
    """
    tolerance = 10**-1

    # download all files
    for file, token in DSS_CKT_FILES.items():
        if not get_local_source_fp(file).exists():
            download_file(file, token)
    for file, token in DSS_CONFIG_FILES.items():
        if not get_local_source_fp(file).exists():
            download_file(file, token)

    # test and compare all files
    if not OUT_DIR.exists():
        os.makedirs(OUT_DIR)

    for file in DSS_CKT_FILES:
        infile = get_local_source_fp(file)
        out_file = Path(OUT_DIR, infile.stem).with_suffix('.out.txt')
        sys.stdout = open(out_file, 'w')
        compare_solution(str(infile), tolerance)


def compare_solution(dss_file, tolerance):
    fbs_solution = SolutionFBS(dss_file)
    fbs_solution.solve()

    dss_solution = SolutionDSS(dss_file)
    dss_solution.solve()

    for param in SolutionFBS.SOLUTION_PARAMS:
        fbs_vals = fbs_solution.get_data_frame(param)
        dss_vals = dss_solution.get_data_frame(param)
        compare_data_frames(fbs_vals, dss_vals, 'fbs', 'dss', param)
        assert ((fbs_vals - dss_vals).abs().max() <= tolerance).all() 


def get_local_source_fp(file: str):
    return Path(OUT_DIR, 'source_files', file)

def download_file(file: str, token: str):
    url = f'{DSS_FILE_REPO}{file}?token={token}'
    download = requests.get(url).text
    local_fp = get_local_source_fp(file)
    if not local_fp.parent.exists():
        os.makedirs(local_fp.parent)
    f = open(get_local_source_fp(file), 'w')
    f.write(download)
