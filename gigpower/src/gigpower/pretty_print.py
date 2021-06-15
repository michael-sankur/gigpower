import pandas as pd
import numpy as np
from . solution import Solution


def compare_solutions(sol_left: Solution, sol_right: Solution, left_name:str,
                      right_name: str):
    """ Prints a comparison of two Solutions."""
    for param in Solution.SOLUTION_PARAMS:
        left_df = sol_left.get_data_frame(param, 'rows')
        right_df = sol_left.get_data_frame(param, 'rows')
        compare_data_frames(left_df, right_df, left_name, right_name, f'{param}')


def compare_data_frames(df1: pd.DataFrame, df2: pd.DataFrame, df1_label:str, 
                        df2_label: str, title: str, cols = ['A', 'B', 'C']):
    """ 
    Prints a comparison of 2 dataframes, df1 - df2
    """
    summary_labels = ['MAX', 'SUM', 'MEAN']
    compare = df1 - df2
    col_suffix = f'{df1_label} - {df2_label}'
    pd.options.display.float_format = '{:.2f}'.format
    if len(compare.axes) == 2:
        compare.columns = [f'{col}.{col_suffix}' for col in cols]
        print(f"{title} - SUMMARY STATS")
        summary_stats = np.asarray([compare.abs().max(), compare.abs().sum(), compare.abs().mean()])
        summary = pd.DataFrame(summary_stats,
            [f'{stat} |{col_suffix}|' for stat in summary_labels], cols)
        print(summary, '\n')
    print(f"{title} - ERROR, |{col_suffix}|:")
    print(compare, "\n")
    print(f"{title} - {df1_label} Results")
    print(df1, "\n")
    print(f"{title} - {df2_label} RESULTS")
    print(df2, "\n")
    print("~"*100, "\n")
