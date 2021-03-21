import pandas as pd


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
