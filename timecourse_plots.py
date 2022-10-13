"""_summary_

_extended_summary_
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def calculate_metrics(results, answers):
    diffdf = results.subtract(answers)
    percentDiffdf = diffdf.divide(answers).multiply(100)
    percentDiffdf.replace([np.inf, -np.inf], np.nan, inplace=True)

    exclude_list = ['Cr', 'PCr', 'PCho', 'GPC']
    diffdf.drop(exclude_list, axis=1, inplace=True)
    percentDiffdf.drop(exclude_list, axis=1, inplace=True)

    metrics = {
        'mm_all': np.mean(diffdf.abs().to_numpy()),
        'per_all': np.nanmean(percentDiffdf.abs().to_numpy())}

    combined_list = ['NAA+NAAG', 'Cr+PCr', 'Glu+Gln', 'Ins+Gly', 'GPC+PCho']
    diffdf = diffdf[combined_list]
    percentDiffdf = percentDiffdf[combined_list]

    metrics.update({
        'mm_high': np.mean(diffdf.abs().to_numpy()),
        'per_high': np.nanmean(percentDiffdf.abs().to_numpy())})

    return pd.Series(metrics)


def generate_tc_plots(metrics_list, versions_list, name, output):
    full_version_df = pd.concat(metrics_list, axis=1, keys=versions_list, names='versions').T

    def plot_lineplot_ver(df, ylabel):
        df_plot = df.stack().reset_index().rename(columns={'level_1': 'metric', 0: 'value'})

        sns.lineplot(
            data=df_plot,
            x='versions',
            y='value',
            hue='metric',
            style='metric',
            markers=True)
        plt.ylabel(ylabel)

    fig1 = plt.figure(figsize=(8, 6))
    plot_lineplot_ver(full_version_df.filter(regex='per'), 'Absolute error (%)')
    per_error_file = output / f'percent_error_{name}.svg'
    fig1.savefig(per_error_file, format='svg', bbox_inches='tight', pad_inches=0.0)
    plt.close()

    fig2 = plt.figure(figsize=(8, 6))
    plot_lineplot_ver(full_version_df.filter(regex='mm'), 'Absolute error (mM)')
    mm_error_file = output / f'mm_error_{name}.svg'
    fig2.savefig(mm_error_file, format='svg', bbox_inches='tight', pad_inches=0.0)
    plt.close()

    root = output.parent
    return {'per_error': str(per_error_file.relative_to(root)), 'mm_error': str(mm_error_file.relative_to(root))}
