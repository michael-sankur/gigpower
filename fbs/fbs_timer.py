from powerflowpy.fbs_detailed_timing import fbs
import numpy as np # type: ignore
import csv

dss_file = 'fbs/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss'
out_csv = './fbs/tests/timings/fbs_detailed_timing.csv'
TRIALS = 20

if __name__ == '__main__':
    with open(out_csv, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(
            [f'TIMING DATA FOR FBS. INPUT FILE: {dss_file}. All timing data in ms.'])  # title
        # column headers
        writer.writerow(['num_trials', 'total_time', 'init_from_dss_time', 'fbs_time', 'topo_sort_time', 'solution_init_time', 'fwd_sweep_total_time',
                         'pre_bwd_sweep_load_update_total_time', 'bwd_sweep_total_time', 'convergence_check_total_time', 'final_calcs_time'])

        for i in range(TRIALS):
            result = np.zeros(11) # initialize result row for this trial set
            for j in range(i+1):
                trial_result = fbs(dss_file)
                total_time = sum(trial_result)
                fbs_time = sum(trial_result[1:]) # get total time of all subroutines except init_from_dss_time
                # add this trial to the result array
                result = result + np.asarray([1, total_time, fbs_time] + trial_result)
            writer.writerow(result)



