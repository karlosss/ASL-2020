#!/usr/bin/env python3

import os
import sys
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import python_lib.csv_cols as csv_cols
import python_lib.constants as constants
import python_lib.cpu_info as cpu_info
import random
import argparse


plt.rcParams.update(plt.rcParamsDefault)
plt.style.use('ggplot')
colormap = 'tab10'

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]


def is_exponential(vals):
    return (
        len(vals) > 3 
        and vals[0] * 2 == vals[1]  
        and vals[1] * 2 == vals[2]  
        and vals[2] * 2 == vals[3]
    )


def extract_NP_data(data, M, T, section):
    f =  data[
          (data[csv_cols.PARAM_M] == M)
        & (data[csv_cols.PARAM_T] == T)
        & (data[csv_cols.VARIABLE] == 0)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_N], f[csv_cols.PERFORMANCE]

def extract_NP_data_cache(data, M, T, section):
    f =  data[
          (data[csv_cols.PARAM_M] == M)
        & (data[csv_cols.PARAM_T] == T)
        & (data[csv_cols.VARIABLE] == 0)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_N], f[csv_cols.MISS_RATE]


def extract_MP_data(data, N, T, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_T] == T)
        & (data[csv_cols.VARIABLE] == 1)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_M], f[csv_cols.PERFORMANCE]

def extract_MP_data_cache(data, N, T, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_T] == T)
        & (data[csv_cols.VARIABLE] == 1)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_M], f[csv_cols.MISS_RATE]


def extract_TP_data(data, N, M, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_M] == M)
        & (data[csv_cols.VARIABLE] == 2)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_T], f[csv_cols.PERFORMANCE]

def extract_TP_data_cache(data, N, M, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_M] == M)
        & (data[csv_cols.VARIABLE] == 2)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_T], f[csv_cols.MISS_RATE]


def adjust_param(data, param, value):
    vals = data[param].unique()
    if value is None:
        max_val = vals.max()
        print(f"Parameter {param} not specified. Setting {param}:={max_val} (max value)")
        return max_val
    elif value not in vals:
        nearest_val = find_nearest(vals, value)
        print(f"No data found for {param}={value}. Setting {param}:={nearest_val} (nearest value)")
        return nearest_val
    else:
        return value


def get_experiment_info(data):
    return data[csv_cols.BINARY].iloc[0], data[csv_cols.COMPILER].iloc[0], data[csv_cols.FLAGS].iloc[0]


def get_sections(data):
    return data[csv_cols.SECTION].unique()

def get_fixed_variables(data):
    N_fix = data[data[csv_cols.VARIABLE] == 1][csv_cols.PARAM_N].iloc[0]
    M_fix = data[data[csv_cols.VARIABLE] == 2][csv_cols.PARAM_M].iloc[0]
    T_fix = data[data[csv_cols.VARIABLE] == 0][csv_cols.PARAM_T].iloc[0]

    return N_fix, M_fix, T_fix


def comparison_data_generator_NP(data_generator, M, T):
    return (extract_NP_data(data, M, T, "baum_welch") for data in data_generator)


def format_plot(ax, xlabel, ylabel, title, is_exp, min_exp=None, max_exp=None):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, rotation=0)
    ax.yaxis.set_label_coords(-0.05, 1.0)
    ax.grid(axis='x')
    #ax.legend(loc='lower right')
    if is_exp:
        ax.set_xticks([2 ** i for i in range(min_exp, max_exp)])
        ax.set_xscale('log', basex=2)
    else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(ymin=-0.2)


def plot_series(ax, x, y, label):
    ax.plot(x, y, '-o', label=label)

def plot_series_color(ax, x, y, label, color):
    ax.plot(x, y, '-o', label=label, color=color)

def plot_NP_sections(ax, data, M, T, colors, sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_NP_data(data, M, T, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_N,
        ylabel=f"Perf [F/C]",
        title=f"M = {M}, T = {T}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )

def plot_NP_sections_cache(ax, data, M, T, colors, sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_NP_data_cache(data, M, T, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_N,
        ylabel=f"Miss Rate",
        title=f"M = {M}, T = {T}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )


def plot_MP_sections(ax, data, N, T, colors, sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_MP_data(data, N, T, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_M,
        ylabel=f"Perf [F/C]",
        title=f"N = {N}, T = {T}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )

def plot_MP_sections_cache(ax, data, N, T, colors, sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_MP_data_cache(data, N, T, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_M,
        ylabel=f"Miss Rate",
        title=f"N = {N}, T = {T}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )

def plot_TP_sections(ax, data, N, M, colors,sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_TP_data(data, N, M, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_T,
        ylabel=f"Perf [F/C]",
        title=f"N = {N}, M = {M}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )

def plot_TP_sections_cache(ax, data, N, M, colors,sections=constants.SECTIONS):
    exp = None
    for i,section in enumerate(sections):
        x, y = extract_TP_data_cache(data, N, M, section)
        plot_series_color(ax, x, y, label=section, color=colors[i])
        if exp is None:
            exp = is_exponential(x.tolist())
            min_exp, max_exp = int(math.log2(x.min())), int(math.log2(x.max()))

    format_plot(ax, 
        xlabel=csv_cols.PARAM_T,
        ylabel=f"Miss Rate",
        title=f"N = {N}, M = {M}",
        is_exp=exp, 
        min_exp=min_exp, 
        max_exp=max_exp
    )


def multiplot_NP_M_comparison(csv_files, N=None, M=None, T=None):
    plt.figure(figsize=(15, 12), facecolor='w')

    fig = plt.gcf()
    fig.suptitle(f"Comparison", fontsize=16)

    ax_NP = plt.subplot(2, 3, 1)
    ax_MP = plt.subplot(2, 3, 2)
    ax_TP = plt.subplot(2, 3, 3)
    ax_table = plt.subplot(2,3,4)

    exp_NP = exp_MP = exp_TP = None

    for csv_file in csv_files:
        data = pd.read_csv(csv_file)

        # get the fixed parameters in the experiment
        N_fix, M_fix, T_fix = get_fixed_variables(data)

        binary_name, compiler, flags = get_experiment_info(data)
        label = f"{binary_name}, {compiler}, {flags}"

        x_NP, y_NP = extract_NP_data(data, M_fix, T_fix, section="baum_welch")
        x_MP, y_MP = extract_MP_data(data, N_fix, T_fix, section="baum_welch")
        x_TP, y_TP = extract_TP_data(data, N_fix, M_fix, section="baum_welch")

        plot_series(ax_NP, x_NP, y_NP, label=label)
        plot_series(ax_MP, x_MP, y_MP, label=label)
        plot_series(ax_TP, x_TP, y_TP, label=label)

        if exp_NP is None:
            exp_NP = is_exponential(x_NP.tolist())
            min_exp_NP, max_exp_NP = int(math.log2(x_NP.min())), int(math.log2(x_NP.max()))
        if exp_MP is None:
            exp_MP = is_exponential(x_MP.tolist())
            min_exp_MP, max_exp_MP = int(math.log2(x_MP.min())), int(math.log2(x_MP.max()))
        if exp_TP is None:
            exp_TP = is_exponential(x_TP.tolist())
            min_exp_TP, max_exp_TP = int(math.log2(x_TP.min())), int(math.log2(x_TP.max()))


    format_plot(ax_NP, 
        xlabel=csv_cols.PARAM_N,
        ylabel=f"Perf [F/C]",
        title=f"M = {M_fix}, T = {T_fix}",
        is_exp=exp_NP, 
        min_exp=min_exp_NP, 
        max_exp=max_exp_NP
    )
    format_plot(ax_MP, 
        xlabel=csv_cols.PARAM_M,
        ylabel=f"Perf [F/C]",
        title=f"N = {N_fix}, T = {T_fix}",
        is_exp=exp_MP, 
        min_exp=min_exp_MP, 
        max_exp=max_exp_MP
    )
    format_plot(ax_TP, 
        xlabel=csv_cols.PARAM_T,
        ylabel=f"Perf [F/C]",
        title=f"N = {N_fix}, M = {M_fix}",
        is_exp=exp_TP, 
        min_exp=min_exp_TP, 
        max_exp=max_exp_TP
    )
    plot_cpu_info_table(ax_table, "CPU Info")
    fig.tight_layout(pad=4.0, rect=[0, 0.0, 1, 0.95])
    handles, labels = ax_NP.get_legend_handles_labels()
    fig.legend(
        handles, 
        labels, 
        loc='center',
        bbox_to_anchor=(0.6, 0., 0.5, 0.5),
        fontsize=16
    )
    plt.show()


def multiplot_NP_M_comparison_cache(csv_files, N=None, M=None, T=None):
    plt.figure(figsize=(15, 12), facecolor='w')

    fig = plt.gcf()
    fig.suptitle(f"Comparison", fontsize=16)

    ax_NP = plt.subplot(2, 3, 1)
    ax_MP = plt.subplot(2, 3, 2)
    ax_TP = plt.subplot(2, 3, 3)
    ax_table = plt.subplot(2,3,4)

    exp_NP = exp_MP = exp_TP = None

    for csv_file in csv_files:
        data = pd.read_csv(csv_file)

        # get the fixed parameters in the experiment
        N_fix, M_fix, T_fix = get_fixed_variables(data)

        binary_name, compiler, flags = get_experiment_info(data)
        label = f"{binary_name}, {compiler}, {flags}"

        x_NP, y_NP = extract_NP_data_cache(data, M_fix, T_fix, section="baum_welch")
        x_MP, y_MP = extract_MP_data_cache(data, N_fix, T_fix, section="baum_welch")
        x_TP, y_TP = extract_TP_data_cache(data, N_fix, M_fix, section="baum_welch")

        plot_series(ax_NP, x_NP, y_NP, label=label)
        plot_series(ax_MP, x_MP, y_MP, label=label)
        plot_series(ax_TP, x_TP, y_TP, label=label)

        if exp_NP is None:
            exp_NP = is_exponential(x_NP.tolist())
            min_exp_NP, max_exp_NP = int(math.log2(x_NP.min())), int(math.log2(x_NP.max()))
        if exp_MP is None:
            exp_MP = is_exponential(x_MP.tolist())
            min_exp_MP, max_exp_MP = int(math.log2(x_MP.min())), int(math.log2(x_MP.max()))
        if exp_TP is None:
            exp_TP = is_exponential(x_TP.tolist())
            min_exp_TP, max_exp_TP = int(math.log2(x_TP.min())), int(math.log2(x_TP.max()))


    format_plot(ax_NP, 
        xlabel=csv_cols.PARAM_N,
        ylabel=f"Miss Rate",
        title=f"M = {M_fix}, T = {T_fix}",
        is_exp=exp_NP, 
        min_exp=min_exp_NP, 
        max_exp=max_exp_NP
    )
    format_plot(ax_MP, 
        xlabel=csv_cols.PARAM_M,
        ylabel=f"Miss Rate",
        title=f"N = {N_fix}, T = {T_fix}",
        is_exp=exp_MP, 
        min_exp=min_exp_MP, 
        max_exp=max_exp_MP
    )
    format_plot(ax_TP, 
        xlabel=csv_cols.PARAM_T,
        ylabel=f"Miss Rate",
        title=f"N = {N_fix}, M = {M_fix}",
        is_exp=exp_TP, 
        min_exp=min_exp_TP, 
        max_exp=max_exp_TP
    )
    plot_cpu_info_table(ax_table, "CPU Info")
    fig.tight_layout(pad=4.0, rect=[0, 0.0, 1, 0.95])
    handles, labels = ax_NP.get_legend_handles_labels()
    fig.legend(
        handles, 
        labels, 
        loc='center',
        bbox_to_anchor=(0.6, 0., 0.5, 0.5),
        fontsize=16
    )
    plt.show()

def multiplot_NP_MP_TP_S(csv_file, save_dir, N=None, M=None, T=None, show_plot=False):
    data = pd.read_csv(csv_file)
    N = adjust_param(data, csv_cols.PARAM_N, N)
    M = adjust_param(data, csv_cols.PARAM_M, M)
    T = adjust_param(data, csv_cols.PARAM_T, T)
    plt.figure(figsize=(15, 12), facecolor='w')
    binary_name,compiler, flags = get_experiment_info(data)

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Compiler: {compiler}, Flags: {flags}", fontsize=16)

    ax_NP = plt.subplot(2, 3, 1)
    ax_MP = plt.subplot(2, 3, 2)
    ax_TP = plt.subplot(2, 3, 3)
    ax_S  = plt.subplot(2, 3, 5)
    ax_table = plt.subplot(2, 3, 4)
    
    sections = get_sections(data)
    cmap = plt.get_cmap(colormap)
    colors = cmap(np.linspace(0, 1, len(sections)))

    plot_NP_sections(ax_NP, data, M, T, colors, sections)
    plot_MP_sections(ax_MP, data, N, T, colors, sections)
    plot_TP_sections(ax_TP, data, N, M, colors, sections)
    plot_cpu_info_table(ax_table, "CPU Info")

    plot_regions_pie(ax_S, data, "Sections", M, T, colors)

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    fig_dir = save_dir
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_name = "performance"
    fig_save_path = os.path.join(fig_dir, fig_name)
    fig = plt.gcf()
    handles, labels = ax_NP.get_legend_handles_labels()
    fig.legend(
        handles, 
        labels, 
        loc='center',
        bbox_to_anchor=(0.6, 0., 0.5, 0.5),
        fontsize=16
    )

    plt.savefig(fig_save_path)
    if(show_plot):
        plt.show()
    print(f"Figure saved to: {fig_save_path}.png")


def multiplot_NP_MP_TP_Cache(csv_file, save_dir, N=None, M=None, T=None, show_plot=False):
    data = pd.read_csv(csv_file)
    N,M,T = get_fixed_variables(data)
    plt.figure(figsize=(15, 12), facecolor='w')
    binary_name,compiler, flags = get_experiment_info(data)

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Compiler: {compiler}, Flags: {flags}", fontsize=16)

    ax_NP = plt.subplot(2, 3, 1)
    ax_MP = plt.subplot(2, 3, 2)
    ax_TP = plt.subplot(2, 3, 3)
    ax_table = plt.subplot(2, 3, 4)
    
    sections = get_sections(data)
    cmap = plt.get_cmap(colormap)
    colors = cmap(np.linspace(0, 1, len(sections)))

    plot_NP_sections_cache(ax_NP, data, M, T, colors, sections)
    plot_MP_sections_cache(ax_MP, data, N, T, colors, sections)
    plot_TP_sections_cache(ax_TP, data, N, M, colors, sections)
    plot_cpu_info_table(ax_table, "CPU Info")

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    fig_dir = save_dir
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_name = "miss_rate"
    fig_save_path = os.path.join(fig_dir, fig_name)
    fig = plt.gcf()
    handles, labels = ax_NP.get_legend_handles_labels()
    fig.legend(
        handles, 
        labels, 
        loc='center',
        bbox_to_anchor=(0.6, 0., 0.5, 0.5),
        fontsize=16
    )

    plt.savefig(fig_save_path)
    if(show_plot):
        plt.show()
    print(f"Figure saved to: {fig_save_path}.png")


def plot_regions_pie(ax, data, title, M_fix, T_fix, colors):
    N_max = data[csv_cols.PARAM_N].max()

    region_perf = data[
          (data[csv_cols.PARAM_N] == N_max)
        & (data[csv_cols.PARAM_M] == M_fix)
        & (data[csv_cols.PARAM_T] == T_fix)
        & (data[csv_cols.VARIABLE] == 0)
        & (data[csv_cols.SECTION] != 'baum_welch')
    ]
    ax.set_title(title)
    _, _, autotexts = ax.pie(
        region_perf[csv_cols.NUM_CYCLES], 
        labels=region_perf[csv_cols.SECTION],
        colors = colors[1:,:],
        autopct='%1.1f%%'
    )
    ax.axis('equal')
    for autotext in autotexts:
        autotext.set_color('white')


def plot_cpu_info_table(ax, title):
    row_labels = [
        "Name:",
        "Number:",
        "Base Freq:",
        "Max Freq:",
        "Turbo Boost:",
        "L1 cache:",
        "L2 cache:",
        "L3 cache:"
    ]
    data = [[item] for item in cpu_info.info_list]
    table = ax.table(
        cellText=data ,
        rowLabels=row_labels,
        rowLoc='left', 
        colLoc='left',
        edges='horizontal',
        loc='center',
        bbox=[0.4, 0.0, 0.7, 1.0])
    table.set_fontsize(18)
    # table.scale(0.5, 3)
    ax.axis('off')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare the performance of the specified experiments')

    parser.add_argument('--experiment_dir', '-e', nargs='+', help='A list of experiment directory paths with data to compare.')
    parser.add_argument('--directory', '-d', help='Directory containing experiment directories with data to compare.')
    parser.add_argument('--recreate', '-r', help='Path to .csv experiment report file to recreate the plot from.')
    
    args = parser.parse_args()
    csv_files = []
    directory = None
    if args.experiment_dir is not None:
        for experiment_dir in args.experiment_dir:
            for filename in os.listdir(experiment_dir):
                if filename.endswith(".csv"):
                    csv_files.append(os.path.join(experiment_dir, filename))

    elif args.recreate is not None:
        csv_path = args.recreate
        data = pd.read_csv(csv_path)
        N_fix, M_fix, T_fix = get_fixed_variables(data)
        dir_path = os.path.dirname(csv_path)
        multiplot_NP_MP_TP_S(csv_path, dir_path, N_fix, M_fix, T_fix, False)
        multiplot_NP_MP_TP_Cache(csv_path, dir_path, N_fix, M_fix, T_fix, False)
        exit(0)
    else:
        if args.directory is not None:
            directory = args.directory
        else:
            print(f"No arguments specified. Comparing all csv files in directory '{constants.OUTPUT_DIR}'.")
            directory = constants.OUTPUT_DIR
        for experiment_dir in os.listdir(directory):
            experiment = os.path.join(directory, experiment_dir)
            for filename in os.listdir(experiment):
                if filename.endswith(".csv"):
                    csv_files.append(os.path.join(experiment, filename))

    multiplot_NP_M_comparison(csv_files)
    multiplot_NP_M_comparison_cache(csv_files)
