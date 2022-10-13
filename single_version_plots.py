"""_summary_

_extended_summary_

"""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# Plot settings
SMALL_SIZE = 10
MEDIUM_SIZE = 16
BIGGER_SIZE = 20

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


plot_order = ['True', 'Estimated']
colors = sns.color_palette(palette='bright', n_colors=1)
colors.insert(0, (0.0, 0.0, 0.0))

default_plot_list = ['Cr+PCr', 'Glu+Gln', 'Ins+Gly', 'GPC+PCho', 'GABA', 'GSH']


def generate_version_plots(results, answers, name, output):
    answers.columns = pd.MultiIndex.from_product(
        [['True'], answers.columns])
    results.columns = pd.MultiIndex.from_product(
        [['Estimated'], results.columns])
    plotDF = results.join(answers)
    plotDF.columns.names = ['valtype', 'dataset']
    plotDF.index.name = 'metabolite'
    plotDF = plotDF.T.reset_index()

    fig_lines = _plot_lineplots(plotDF)
    metabolite_plots_file = output / ('lineplt_' + name + '.svg')
    fig_lines.savefig(
        metabolite_plots_file,
        format='svg',
        bbox_inches='tight',
        pad_inches=0.0)
    plt.close()

    fig_mat = _plot_matrix(results, answers)
    matrix_plot_file = output / ('matplt_' + name + '.svg')
    fig_mat.savefig(
        matrix_plot_file,
        format='svg',
        bbox_inches='tight',
        pad_inches=0.0)
    plt.close()

    root = output.parent
    return {
        'matrix_plot': str(matrix_plot_file.relative_to(root)),
        'metabolite_plots': str(metabolite_plots_file.relative_to(root))}


def _lineplot(ax, data_df, metab, annotate=False):
    vpos = [2.5, 5.5, 9.5, 10.5, 13.5, 16.5, 19.5]
    for pos in vpos:
        ax.axvline(x=pos, color=(0.4, 0.4, 0.4), ls=':')

    sns.lineplot(
        x='dataset',
        y=metab,
        data=data_df,
        hue='valtype',
        style='valtype',
        hue_order=plot_order,
        palette=colors,
        markers=True,
        dashes=True,
        legend=annotate,
        ax=ax)

    ax.set_title(metab)
    ax.set_ylabel('Conc (mM)')
    ax.set_xticks(np.arange(0, 21))
    if annotate:
        ax.set_ylim([8, 25])
        ax.legend(loc='upper right')
        annoYpos = \
            ax.get_ylim()[0] + 0.05 * (ax.get_ylim()[1] - ax.get_ylim()[0])

        annotationlist = [
            ('Lorentzian\n(SNR = 80)', (-0.25, annoYpos)),
            ('Gaussian\n(SNR = 80)', (2.75, annoYpos)),
            ('GABA/GSH\n(SNR =160/80)', (5.75, annoYpos)),
            ('No\nMM', (9.55, annoYpos)),
            ('SNR=20', (11.0, annoYpos)),
            ('SNR=30', (14.0, annoYpos)),
            ('SNR=40', (17.0, annoYpos)),
            ('SNR\n=160', (19.6, annoYpos))]
        for anno in annotationlist:
            plt.annotate(anno[0], anno[1], fontsize=12)
        ax.set_xticks(np.arange(0, 21))
    else:
        ax.set_xticks([])
        ax.set_xlabel('')


def _plot_lineplots(plotDF, plot_list=default_plot_list):
    # Figure setup
    fig = plt.figure(figsize=(10, 11))
    gs = fig.add_gridspec(4, 2)
    f1_ax1 = fig.add_subplot(gs[0, 0:2])

    # Subplot 1
    _lineplot(f1_ax1, plotDF, 'NAA+NAAG', annotate=True)

    # Subsidary plots
    plt.subplots_adjust(hspace=0.3)
    axesPos = [gs[1, 0], gs[1, 1], gs[2, 0], gs[2, 1], gs[3, 0], gs[3, 1]]
    axes = [fig.add_subplot(ap) for ap in axesPos]

    for currax, metab in zip(axes, plot_list):
        _lineplot(currax, plotDF, metab, annotate=False)
    return fig


def _showMat(
        ax,
        res,
        ans,
        absolute=True,
        percent=True,
        sort=True,
        title=None,
        disable_xlabel=False):
    if percent:
        df = res.divide(ans).multiply(100)
        df.replace([np.inf, -np.inf], np.nan, inplace=True)
        vmin = -50.0
        vmax = 50.0
    else:
        df = res
        vmin = None
        vmax = None
    if sort:
        tosort = df.abs().mean().to_numpy()
        SI = np.argsort(tosort)
        df = df[df.columns[SI]]
        mat = df.to_numpy().T
    if absolute:
        mat = np.abs(mat)
        vmin = 0.0
        cmap = plt.cm.Greens
    else:
        cmap = plt.cm.PiYG

    cax = ax.matshow(mat, cmap=cmap, vmin=vmin, vmax=vmax)
    # ax.set_xticklabels(res.index.to_numpy())
    if disable_xlabel:
        ax.set_xticks([])
    else:
        ax.set_xticks(np.arange(res.shape[0]))
        ax.set_xticklabels(res.index.to_numpy())
        ax.set_xlabel('Dataset')
        ax.xaxis.set_label_position('top')
    ax.set_yticks(np.arange(0, len(res.columns)))
    ax.set_yticklabels(df.columns)
    if disable_xlabel:
        ax.set_title(title, y=1.0)
    else:
        ax.set_title(title, y=1.08)

    return ax, cax


def _plot_matrix(results, answers):
    fig = plt.figure(figsize=(6, 5))

    MHforSub = results.filter(like='Estimated').droplevel(0, axis=1).T
    MHforSub.drop(['Mac'], axis=1, inplace=True)
    ansforSub = answers.droplevel(0, axis=1).T
    ansforSub.drop(['lipids', 'Ace'], axis=1, inplace=True)
    diffdf = MHforSub.subtract(ansforSub)

    _, cax = _showMat(
        fig.gca(),
        diffdf,
        ansforSub,
        absolute=False,
        percent=True,
        title='Estimated')

    axColor = plt.axes([0.8, 0.21, 0.02, 0.58])
    plt.colorbar(cax, cax=axColor, orientation="vertical")
    axColor.set_title('%')

    return fig
