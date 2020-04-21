#!/usr/bin/python3

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv_cols
import random

GENERATE = False

SECTIONS = [
    'baum_welch',
    'forward_vars',
    'backward_vars',
    'update_initial',
    'update_transition',
    'update_emission'
] 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

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


def comparison_data_generator_NP(data_generator, M, T):
    return (extract_NP_data(data, M, T, "baum_welch") for data in data_generator)


def plot_NP_sections(ax, data, M, T, sections=SECTIONS):
    for section in sections:
        x, y = extract_NP_data(data, M, T, section)
        ax.set_xticks(x)
        ax.plot(x, y, label=section)
    # ax formatting
    ax.set_title(f"M = {M}, T = {T}")
    ax.set_xlabel(csv_cols.PARAM_N)
    ax.set_ylabel(f"Perf [F/C]", rotation=0)
    ax.yaxis.set_label_coords(-0.05, 1.0)

    ax.grid(axis='x')
    ax.legend(loc='lower right')


def plot_MP_sections(ax, data, N, T, sections=SECTIONS):
    for section in sections:
        x, y = extract_MP_data(data, N, T, section)
        ax.set_xticks(x)
        ax.plot(x, y, label=section)
    # ax formatting
    ax.set_title(f"N = {N}, T = {T}")
    ax.set_xlabel(csv_cols.PARAM_M)
    ax.set_ylabel(f"Perf [F/C]", rotation=0)
    ax.yaxis.set_label_coords(-0.05, 1.0)

    ax.grid(axis='x')
    ax.legend(loc='lower right')


def get_experiment_info(data):
    return data[csv_cols.BINARY].iloc[0], data[csv_cols.FLAGS].iloc[0]


# def plot_NP_comparison(ax, data_it):
#     # ax formatting
#     for data in data_it:
#         binary, _ = get_experiment_info(data)
#         ax.plot(*extract_NP_data(data, M, T, section="baum_welch"), label=binary)



def plot_NP_MP_S(data, title, N, M, T=None):
    N = adjust_param(data, csv_cols.PARAM_N, N)
    M = adjust_param(data, csv_cols.PARAM_M, M)
    T = adjust_param(data, csv_cols.PARAM_T, T)
    plt.figure(figsize=(15, 6), facecolor='w')
    binary_name, flags = get_experiment_info(data)

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Flags: {flags}", fontsize=16)

    ax_NP = plt.subplot(1, 3, 1)
    ax_MP = plt.subplot(1, 3, 2)
    ax_S  = plt.subplot(1, 3, 3)
 
    plot_NP_sections(ax_NP, data, M, T)
    plot_MP_sections(ax_MP, data, N, T)


    plot_regions_pie(ax_S, data, "Sections")

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    plt.show()
    return


def plot_regions_pie(ax, data, title):
    N_max = data[csv_cols.PARAM_N].max()
    M_max = data[csv_cols.PARAM_M].max()
    T_max = data[csv_cols.PARAM_T].max()

    region_perf = data[
          (data[csv_cols.PARAM_N] == N_max)
        & (data[csv_cols.PARAM_M] == M_max)
        & (data[csv_cols.PARAM_T] == T_max)
    ]
    
    ax.set_title(title)
    _, _, autotexts = ax.pie(region_perf[csv_cols.PERFORMANCE], labels=region_perf[csv_cols.SECTION], autopct='%1.1f%%')
    for autotext in autotexts:
        autotext.set_color('white')


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



def generate_dummie_values(binary, flags, N_range, M_range, T_range, iterations, sections):
    cols = [
        csv_cols.BINARY, 
        csv_cols.FLAGS, 
        csv_cols.PARAM_N, 
        csv_cols.PARAM_M,
        csv_cols.PARAM_T,
        csv_cols.NUM_ITERATIONS,
        csv_cols.SECTION,
        csv_cols.NUM_CYCLES,
        csv_cols.PERFORMANCE
    ]

    df = pd.DataFrame(columns=cols)
    for n in range(*N_range):
        for m in range(*M_range):
            for t in range(*T_range):
                for section in sections:
                    df = df.append({
                        csv_cols.BINARY : binary,
                        csv_cols.FLAGS : flags,
                        csv_cols.PARAM_N : n,
                        csv_cols.PARAM_M : m,
                        csv_cols.PARAM_T : t,
                        csv_cols.NUM_ITERATIONS : iterations,
                        csv_cols.SECTION : section,
                        csv_cols.NUM_CYCLES : random.randint(9000, 10000),
                        csv_cols.PERFORMANCE : random.uniform(0.0, 1.0)
                    }, ignore_index=True)
    return df



if __name__ == "__main__":

    if GENERATE:
        data = generate_dummie_values(
            binary="no_opt",
            flags="-ftree-no-vectorize",
            N_range=(0,5),
            M_range=(0,10,2),
            T_range=(0,2),
            iterations=100,
            sections=[
                "baum_welch", 
                "forward_vars", 
                "backward_vars", 
                "update_initial", 
                "update_transition", 
                "update_emission"
            ]
        )
        data.to_csv('data/dummie_df.csv', index=False)


    data = pd.read_csv('../output_data/no_opt%O0.csv')
    data2 = pd.read_csv('../output_data/no_opt%-O3.csv')

    plt.rcParams.update(plt.rcParamsDefault)
    plt.style.use('ggplot')

    # plot_perf_NM_and_section_pie_chart(data, "Performance Plot Variants", N=27, M=37, all_sections=True)

    plot_NP_MP_S(data, "Blubb",  N=27, M=37)


    

    
