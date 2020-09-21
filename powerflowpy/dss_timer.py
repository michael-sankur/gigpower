from powerflowpy.dss_solve_detailed_timing import *
import numpy as np
import csv

dss_file = 'powerflowpy/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss'
out_csv = './powerflowpy/tests/timings/dss_detailed_timing.csv'
TRIALS = 20
if __name__ == '__main__':
    with open(out_csv, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(
            [f'TIMING DATA FOR DSS. INPUT FILE: {dss_file}. All timing data in ms.'])  # title
        # column headers
        writer.writerow(['num_iterations', 'total_time', 'init_circuit_time', 'save_loads_time',
                         'reset_loads_total_time', 'run_solvesnap_total_time', 'final_calcs_time'])

        for i in range(TRIALS):
            result = np.zeros(7) # initialize result row for this trial set
            for j in range(i+1):
                trial_result = solve_with_dss(dss_file)
                total_time = sum(trial_result)
                # add this trial to the result array
                result = result + np.asarray([1, total_time] + trial_result)
            writer.writerow(result)



