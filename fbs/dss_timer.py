from fbs.dss_solve import setup, solve
import numpy as np  # type: ignore
import csv
import time

dss_file = 'fbs/tests/06n3ph_rad_unbal/06node_threephase_radial_unbalanced.dss'
out_csv = './fbs/tests/timings/dss_detailed_timing.csv'
TRIALS = 20
if __name__ == '__main__':
    with open(out_csv, 'w') as out:
        writer = csv.writer(out)
        writer.writerow(
            [f'TIMING DATA FOR DSS. INPUT FILE: {dss_file}. All timing data in ms.'])  # title
        # column headers
        writer.writerow(['num_trials', 'setup_time', 'solve_time', 'total_time'])

        for i in range(TRIALS):
            result = np.zeros(4)  # initialize result row for this trial set
            for j in range(i+1):

                t1 = time.perf_counter()
                setup(dss_file)
                t2 = time.perf_counter()

                solve()
                t3 = time.perf_counter()

                trial_result = [(t2-t1) * 1000, (t3 - t2) * 1000, (t3 - t1) * 1000]
                # add this trial to the result array
                result = result + np.asarray([1] + trial_result)  # type: ignore
            writer.writerow(result)
