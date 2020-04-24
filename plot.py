#!/usr/bin/python3

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
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_N], f[csv_cols.PERFORMANCE]


def extract_MP_data(data, N, T, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_T] == T)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_M], f[csv_cols.PERFORMANCE]


def extract_TP_data(data, N, M, section):
    f =  data[
          (data[csv_cols.PARAM_N] == N)
        & (data[csv_cols.PARAM_M] == M)
        & (data[csv_cols.SECTION] == section)
    ]
    return f[csv_cols.PARAM_T], f[csv_cols.PERFORMANCE]


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
    return data[csv_cols.BINARY].iloc[0], data[csv_cols.FLAGS].iloc[0]


def get_sections(data):
    return data[csv_cols.SECTION].unique()


def comparison_data_generator_NP(data_generator, M, T):
    return (extract_NP_data(data, M, T, "baum_welch") for data in data_generator)


def format_plot(ax, xlabel, ylabel, title, is_exp, min_exp=None, max_exp=None):
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, rotation=0)
    ax.yaxis.set_label_coords(-0.05, 1.0)
    ax.grid(axis='x')
    ax.legend(loc='lower right')
    if is_exp:
        ax.set_xticks([2 ** i for i in range(min_exp, max_exp)])
        ax.set_xscale('log', basex=2)
    else:
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(ymin=-0.2)


def plot_series(ax, x, y, label):
    ax.plot(x, y, '-o', label=label)

def plot_NP_sections(ax, data, M, T, sections=constants.SECTIONS):
    exp = None
    for section in sections:
        x, y = extract_NP_data(data, M, T, section)
        plot_series(ax, x, y, label=section)
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


def plot_MP_sections(ax, data, N, T, sections=constants.SECTIONS):
    exp = None
    for section in sections:
        x, y = extract_MP_data(data, N, T, section)
        plot_series(ax, x, y, label=section)
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


def plot_TP_sections(ax, data, N, M, sections=constants.SECTIONS):
    exp = None
    for section in sections:
        x, y = extract_TP_data(data, N, M, section)
        plot_series(ax, x, y, label=section)
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


def multiplot_NP_M_comparison(csv_files, N=None, M=None, T=None):
    plt.figure(figsize=(15, 6), facecolor='w')

    fig = plt.gcf()
    fig.suptitle(f"Comparison", fontsize=16)

    ax_NP = plt.subplot(1, 3, 1)
    ax_MP = plt.subplot(1, 3, 2)
    ax_TP = plt.subplot(1, 3, 3)
    
    exp_NP = exp_MP = exp_TP = None

    for csv_file in csv_files:
        data = pd.read_csv(csv_file)

        N = adjust_param(data, csv_cols.PARAM_N, N)
        M = adjust_param(data, csv_cols.PARAM_M, M)
        T = adjust_param(data, csv_cols.PARAM_T, T)

        binary_name, flags = get_experiment_info(data)
        label = f"{binary_name}, {flags}"

        x_NP, y_NP = extract_NP_data(data, M, T, section="baum_welch")
        x_MP, y_MP = extract_MP_data(data, N, T, section="baum_welch")
        x_TP, y_TP = extract_TP_data(data, N, M, section="baum_welch")

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
        title=f"M = {M}, T = {T}",
        is_exp=exp_NP, 
        min_exp=min_exp_NP, 
        max_exp=max_exp_NP
    )
    format_plot(ax_MP, 
        xlabel=csv_cols.PARAM_M,
        ylabel=f"Perf [F/C]",
        title=f"N = {N}, T = {T}",
        is_exp=exp_MP, 
        min_exp=min_exp_MP, 
        max_exp=max_exp_MP
    )
    format_plot(ax_TP, 
        xlabel=csv_cols.PARAM_T,
        ylabel=f"Perf [F/C]",
        title=f"N = {N}, M = {M}",
        is_exp=exp_TP, 
        min_exp=min_exp_TP, 
        max_exp=max_exp_TP
    )
    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    plt.show()


def multiplot_NP_MP_TP_S(csv_file, N=None, M=None, T=None):
    data = pd.read_csv(csv_file)
    N = adjust_param(data, csv_cols.PARAM_N, N)
    M = adjust_param(data, csv_cols.PARAM_M, M)
    T = adjust_param(data, csv_cols.PARAM_T, T)
    plt.figure(figsize=(15, 12), facecolor='w')
    binary_name, flags = get_experiment_info(data)

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Flags: {flags}", fontsize=16)

    ax_NP = plt.subplot(2, 3, 1)
    ax_MP = plt.subplot(2, 3, 2)
    ax_TP = plt.subplot(2, 3, 3)
    ax_S  = plt.subplot(2, 3, 5)
    ax_table = plt.subplot(2, 3, 4)
    
    sections = get_sections(data)
    plot_NP_sections(ax_NP, data, M, T, sections)
    plot_MP_sections(ax_MP, data, N, T, sections)
    plot_TP_sections(ax_TP, data, N, M, sections)
    plot_cpu_info_table(ax_table, "CPU Info")

    plot_regions_pie(ax_S, data, "Sections")

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    fig_dir = "figures"
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)
    fig_name, suffix = os.path.splitext(os.path.basename(csv_file))
    fig_save_path = os.path.join(fig_dir, fig_name)
    plt.savefig(fig_save_path)
    plt.show()
    print(f"Figure saved to: {fig_save_path}{suffix}")


def plot_regions_pie(ax, data, title):
    N_max = data[csv_cols.PARAM_N].max()
    M_max = data[csv_cols.PARAM_M].max()
    T_max = data[csv_cols.PARAM_T].max()

    region_perf = data[
          (data[csv_cols.PARAM_N] == N_max)
        & (data[csv_cols.PARAM_M] == M_max)
        & (data[csv_cols.PARAM_T] == T_max)
        & (data[csv_cols.SECTION] != 'baum_welch')
    ]
    
    ax.set_title(title)
    _, _, autotexts = ax.pie(
        region_perf[csv_cols.PERFORMANCE]/region_perf[csv_cols.PERFORMANCE].sum(), 
        labels=region_perf[csv_cols.SECTION], 
        autopct='%1.1f%%'
    )
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
        bbox=[0.4, 0.0, 0.6, 1.0])
    table.set_fontsize(16)
    # table.scale(0.5, 3)
    ax.axis('off')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare the performance of the specified experiments')

    parser.add_argument('--csv_files', '-f', nargs='+', help='A list of csv files with data to compare.')
    parser.add_argument('--directory', '-d', help='Directory containing csv files with data to compare.')
    
    args = parser.parse_args()
    csv_files = []
    directory = None
    if (args.csv_files is not None):
        csv_files = args.csv_files
    else:
        if (args.directory is not None):
            directory = args.directory
        else:
            print(f"No arguments specified. Comparing all csv files in directory '{constants.OUTPUT_DIR}'.")
            directory = constants.OUTPUT_DIR

        for filename in os.listdir(directory):
            if filename.endswith(".csv"):
                csv_files.append(os.path.join(directory, filename))

    multiplot_NP_M_comparison(csv_files)