# Elaine Laguerta (github: @elaguerta)
# LBNL GIG
# File created: September 2020
# Compare solutions from opendss and powerflowpy.fbs

from fbs.utils import init_from_dss
from fbs.fbs import fbs, get_solution as get_fbs_solution
from fbs.dss_solve import solve_with_dss, getVMag, getNominalNodePower, getTotalNodePowers

import pandas as pd
import numpy as np

import sys
import os.path
import requests

DSS_FILE_REPO = 'https://raw.githubusercontent.com/lbnl-cybersecurity/python-powerflow/master/'
OUT_DIR = './test_compare_results'

# dictionaries with FILE_STEM: API TOKEN
# API token is used to download source files from raw.githubuser.com
DSS_CKT_FILES = {
                'IEEE_13_Node/IEEE_13_Bus_allwye.dss': 'ADCUMKN5VISDG4XWF2CHE6DAIFWYI',
                'IEEE_13_Node/IEEE_13_Bus_allwye_noxfm_noreg.dss': 'ADCUMKM27O7EDMUSGZZXHEDAIF2EQ',
                'IEEE_34_feeder_UB/IEEE_34_Bus_allwye.dss': 'ADCUMKPSDEI4H6FJ64E3WP3AIF2V6',
                'IEEE_34_feeder_UB/IEEE_34_Bus_allwye_noxfm_noreg.dss': 'ADCUMKKKHOWAGZAFIYF74K3AIF2ZG',
                'IEEE_37_feeder_UB/IEEE_37_Bus_allwye.dss': 'ADCUMKMRXJ6XFKXF2MPS4OTAIF3BE',
                'IEEE_37_feeder_UB/IEEE_37_Bus_allwye_noxfm_noreg.dss': 'ADCUMKP36GFYRQDGQKDLZK3AIF3FA'
                }
DSS_CONFIG_FILES = {
                    'IEEE_13_Node/IEEELineCodes.dss': 'ADCUMKKVFD42HGKO3CSUIXTAIFXXY',
                    'IEEE_13_Node/IEEE13Node_BusXY.csv': 'ADCUMKI7QNQV76IMUUMPLEDAIFZXI',
                    'IEEE_34_feeder_UB/IEEE34_BusXY.csv': 'ADCUMKPOJ2UU4GZUNGDS6B3AIF2NC',
                    'IEEE_34_feeder_UB/IEEELineCodes.dss': 'ADCUMKMONSF4CEDEPYC2RUTAIF2RG',
                    'IEEE_37_feeder_UB/IEEELineCodes.dss': 'ADCUMKNPPPEQ3RYI6MA7CK3AIF26M'
                    }


def get_local_fp(file: str):
    return os.path.join(OUT_DIR, file.split('/')[-1])


def download_file(file: str, token: str):
    url = f'{DSS_FILE_REPO}{file}?token={token}'
    download = requests.get(url).text
    f = open(get_local_fp(file), 'w')
    f.write(download)


def test_all():
    """
    Compare the python FBS solution to the opendss solution.
    Writes output to OUTPUT FOLDER.
    """
    tolerance = 10**-1

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)

    # download all files
    for file, token in DSS_CKT_FILES.items():
        download_file(file, token)
    for file, token in DSS_CONFIG_FILES.items():
        download_file(file, token)

    # test and compare all files
    for file in DSS_CKT_FILES:
        infile = get_local_fp(file)
        out_file = file.split('/')[-1].split('.')[0] + '.out.txt'

        sys.stdout = open(os.path.join(OUT_DIR, out_file), 'w')
        compare_fbs_sol(infile, tolerance)


# construct the python FBS solution
def compare_fbs_sol(dss_file, tolerance):
    # TODO: Figure out how to get Inode (iNR) and sV (sNR) at each node from dss and perform a compare
    print(f'\n\nTesting File: {dss_file}')
    # solve with fbs
    network = init_from_dss(dss_file)
    fbs_sol = get_fbs_solution(fbs(network))

    # solve with dss
    dss_sol = solve_with_dss(dss_file)

    fbsV, fbsI, fbsStx, fbsSrx, fbs_sV = fbs_sol.V_df(), fbs_sol.I_df(), fbs_sol.Stx_df(), \
        fbs_sol.Srx_df(), fbs_sol.sV_df()
    dssV, dssI, dssStx, dssSrx, dssNodePowers, dss = dss_sol

    print(f"FBS iterations: {fbs_sol.iterations}\t FBS convergence:\
        {fbs_sol.diff}\t FBS tolerance: {fbs_sol.tolerance}")
    V_maxDiff = compare_dfs(fbsV, dssV, "COMPARE V")

    dssVMag = getVMag(dss)
    fbsVMag = fbs_sol.VMag_df()
    VMag_maxDiff = compare_dfs(fbsVMag, dssVMag, "COMPARE VMag")

    I_maxDiff = compare_dfs(fbsI, dssI, "COMPARE I")

    Stx_maxDiff = compare_dfs(fbsStx, dssStx, "COMPARE Stx")

    Srx_maxDiff = compare_dfs(fbsSrx, dssSrx, "COMPARE Srx")

    fbsNomPwrs = fbs_sol.nomNodePwrs_df() * 1000  # convert to kW
    dssNomPwrs = getNominalNodePower(dss)
    compare_dfs(fbsNomPwrs, dssNomPwrs, "TOTAL NOMINAL NODE POWERS (kW/kVAR)")

    print("Total node power comparison, from dss.CktElement.Powers()")
    fbsNodePowers = fbs_sV * 1000  # convert to kW
    nodePowers_maxDiff1 = compare_dfs(fbsNodePowers, dssNodePowers, "TOTAL NODE POWERS (kW/kVAR)")
    print("Total node power comparison, from dss.CktElement.Powers()")
    compare_dfs(fbs_sV, dssNodePowers * 1000/network.Sbase, "TOTAL NODE POWERS(pu)")

    print("Total node power comparison, from opendss voltage solution")
    dssNodePowers_2 = getTotalNodePowers(dss) * 1000
    nodePowers_maxDiff2 = compare_dfs(fbsNodePowers, dssNodePowers_2, "TOTAL NODE POWERS (kW/kVAR)")
    print("Total node power comparison, from opendss voltage solution")
    compare_dfs(fbs_sV, dssNodePowers_2 * 1000/network.Sbase, "TOTAL NODE POWERS(pu)")

    fbsNodePowers_sumPhase = fbsNodePowers.sum(axis=0)
    dssNodePowers_sumPhase = dssNodePowers.sum(axis=0)
    compare_dfs(fbsNodePowers_sumPhase, dssNodePowers_sumPhase, "TOTAL NODE POWERS - SUM OVER PHASE (kW/kVAR)")
    compare_dfs(fbs_sV.sum(axis=0), dssNodePowers_sumPhase * 1000 /
                network.Sbase, "TOTAL NODE POWERS = SUM OVER PHASE (pu)")

    # assert (V_maxDiff <= tolerance).all()
    # assert (VMag_maxDiff <= tolerance).all()
    # assert (I_maxDiff <= tolerance).all()
    # assert (Stx_maxDiff <= tolerance).all()
    # assert (Srx_maxDiff <= tolerance).all()
    # assert (nodePowers_maxDiff1/1000 <= tolerance).all()
    # assert (nodePowers_maxDiff2/1000 <= tolerance).all()


def compare_dfs(fbs_df: pd.DataFrame, dss_df: pd.DataFrame, title: str) -> None:
    """ helper method to compare fbs vs. dss and print comparisons """
    compare = fbs_df - dss_df
    pd.options.display.float_format = '{:.4f}'.format
    if len(compare.axes) == 2:
        compare_cols = ['A.(fbs - dss)', 'B.(fbs - dss)', 'C.(fbs - dss)']
        compare.columns = compare_cols
        print(f"{title} - SUMMARY STATS")
        summary_stats = np.asarray([compare.abs().max(), compare.abs().sum(), compare.abs().mean()])
        summary = pd.DataFrame(summary_stats,
            ['MAX |fbs - dss|', 'SUM |fbs - dss|', 'MEAN |fbs - dss|'], ['A', 'B', 'C'])
        print(summary, '\n')
    print(f"{title} - ERROR, |fbs - dss|:")
    print(compare, "\n")
    print(f"{title} - FBS RESULTS")
    print(fbs_df, "\n")
    print(f"{title} - DSS RESULTS")
    print(dss_df, "\n")
    print("~"*100, "\n")
    return compare.abs().max()
