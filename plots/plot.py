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

def extract_N_plot_data(data, M, T, section="baum_welch"):
    M_constraint, T_constraint = data[csv_cols.PARAM_M] == M, data[csv_cols.PARAM_T] == T
    return data[M_constraint & T_constraint]

def extract_M_plot_data(data, N, T, section="baum_welch"):
    N_constraint, T_constraint = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_T] == T
    return data[N_constraint & T_constraint]
    

def extract_data(data, x_axis, y_axis=csv_cols.PERFORMANCE, section="baum_welch", N=None, M=None, T=None):
    if x_axis == csv_cols.PARAM_N:
        fixed_0, fixed_1 = data[csv_cols.PARAM_M] == M, data[csv_cols.PARAM_T] == T

    elif x_axis == csv_cols.PARAM_M: 
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_T] == T

    elif x_axis == csv_cols.PARAM_T:
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_M] == M
    
    section_constraint = data[csv_cols.SECTION] == section
    return  data[fixed_0 & fixed_1 & section_constraint]


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



def plot_perf_NM_and_section_pie_chart(data, title, N, M, T=None, all_sections=False):
    N = adjust_param(data, csv_cols.PARAM_N, N)
    M = adjust_param(data, csv_cols.PARAM_M, M)
    T = adjust_param(data, csv_cols.PARAM_T, T)
    plt.figure(figsize=(15, 6), facecolor='w')
    binary_name = data[csv_cols.BINARY].iloc[0]
    flags = data[csv_cols.FLAGS].iloc[0]

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Flags: {flags}", fontsize=16)

    ax_N = plt.subplot(1, 3, 1)
    ax_M = plt.subplot(1, 3, 2)
    ax_T = plt.subplot(1, 3, 3)


    sections = SECTIONS if all_sections else ["baum_welch"]
    plot_N_P(ax_N, data, "Plot N", M, T, sections=sections)
    plot_M_P(ax_M, data, "Plot M", N, T, sections=sections)
    plot_regions_pie(ax_T, data, "Sections")

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    plt.show()
    return


def plot_perf_NM(data, title, N, M, T=None):
    N = adjust_param(data, csv_cols.PARAM_N, N)
    M = adjust_param(data, csv_cols.PARAM_M, M)
    T = adjust_param(data, csv_cols.PARAM_T, T)
    plt.figure(figsize=(15, 6), facecolor='w')
    binary_name = data[csv_cols.BINARY].iloc[0]
    flags = data[csv_cols.FLAGS].iloc[0]

    fig = plt.gcf()
    fig.suptitle(f"Binary: {binary_name}, Flags: {flags}", fontsize=16)

    ax_N = plt.subplot(1, 2, 1)
    ax_M = plt.subplot(1, 2, 2)

    plot_N_P(ax_N, data, "Plot N", M, T)
    plot_M_P(ax_M, data, "Plot M", N, T)

    fig.tight_layout(pad=3.0, rect=[0, 0.0, 1, 0.95])
    plt.show()
    return


def plot_N_P(ax, data, title, M, T, sections=["baum_welch"]):
    ax.set_title(f"M = {M}, T = {T}")
    plot_data(ax, data, csv_cols.PARAM_N, csv_cols.PERFORMANCE, sections=sections, M=M, T=T)

    
def plot_M_P(ax, data, title, N, T, sections=["baum_welch"]):
    ax.set_title(f"M = {N}, T = {T}")
    plot_data(ax, data, csv_cols.PARAM_M, csv_cols.PERFORMANCE, sections=sections, N=N, T=T)


def plot_T_P(ax, data, title, N, M, sections=["baum_welch"]):
    ax.set_title(title)
    plot_data(ax, data, csv_cols.PARAM_M, csv_cols.PERFORMANCE, sections=sections, N=N, M=M)



def plot_data(ax, data, x_axis, y_axis=csv_cols.PERFORMANCE, sections=["baum_welch"], N=None, M=None, T=None):
    for section in sections:
        extracted_data = extract_data(data, x_axis, y_axis, section, N, M, T)
        # ax.set_title(extracted_data[csv_cols.BINARY].iloc[0])
        ax.set_xlabel(x_axis)
        ax.set_ylabel(f"Perf [F/C]", rotation=0)
        ax.yaxis.set_label_coords(-0.05, 1.0)
        ax.set_xticks(extracted_data[x_axis])
        ax.plot(extracted_data[x_axis], extracted_data[y_axis], '-o', label=section)
    ax.grid(axis='x')
    ax.legend(loc='lower right')


# def compare_data(csv_files):
#     ax_N = plt.subplot(1, 2, 1)
#     ax_M = plt.subplot(1, 2, 2)

#     for csv_file in csv_files:
#         data = pd.read_csv(csv_file)
#         plot_N_P(ax_N, data, )


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

    plot_perf_NM_and_section_pie_chart(data, "Performance Plot Variants", N=27, M=37, all_sections=True)

    


    

    
