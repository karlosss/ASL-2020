#!/usr/bin/python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import csv_cols
import random

GENERATE = False


def extract_data(data, x_axis, y_axis=csv_cols.PERFORMANCE, section="baum_welch", N=None, M=None, T=None):
    if x_axis == csv_cols.PARAM_N:
        if M is None: 
            M = data[csv_cols.PARAM_M].max()
            print(f"Warning: No M specified. Setting M := {M}")
        if T is None: 
            T = data[csv_cols.PARAM_T].max()
            print(f"Warning: No T specified. Setting T := {T}")
        fixed_0, fixed_1 = data[csv_cols.PARAM_M] == M, data[csv_cols.PARAM_T] == T

    elif x_axis == csv_cols.PARAM_M: 
        if N is None: 
            N = data[csv_cols.PARAM_N].max()
            print(f"Warning: No N specified. Setting N := {N}")
        if T is None: 
            T = data[csv_cols.PARAM_T].max()
            print(f"Warning: No T specified. Setting T := {T}")
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_T] == T

    elif x_axis == csv_cols.PARAM_T:
        if M is None: 
            M = data[csv_cols.PARAM_M].max()
            print(f"Warning: No M specified. Setting M := {M}")
        if N is None: 
            N = data[csv_cols.PARAM_N].max()
            print(f"Warning: No N specified. Setting N := {N}")
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_M] == M
    
    section_constraint = data[csv_cols.SECTION] == section
    return  data[fixed_0 & fixed_1 & section_constraint]

    # plt.plot(extracted_data[x_axis], extracted_data[csv_cols.PERFORMANCE])
    # plt.show()    

def plot_perf_NMT(data, title, N, M, T):
    plt.figure(figsize=(20, 10), facecolor='w')
    ax_N = plt.subplot(1, 3, 1)
    ax_M = plt.subplot(1, 3, 2)
    ax_T = plt.subplot(1, 3, 3)
    plot_N_P(ax_N, data, "Plot N", M, T)
    plot_M_P(ax_M, data, "Plot M", N, T)
    plot_T_P(ax_T, data, "Plot T", N, M)
    plt.show()
    return

def plot_N_P(ax, data, title, M, T):
    ax.set_title(title)
    plot_data(ax, data, csv_cols.PARAM_N, csv_cols.PERFORMANCE, section="baum_welch", M=M, T=T)

    
def plot_M_P(ax, data, title, N, T):
    ax.set_title(title)
    plot_data(ax, data, csv_cols.PARAM_M, csv_cols.PERFORMANCE, section="baum_welch", N=N, T=T)


def plot_T_P(ax, data, title, N, M):
    ax.set_title(title)
    plot_data(ax, data, csv_cols.PARAM_M, csv_cols.PERFORMANCE, section="baum_welch", N=N, M=M)


def plot_data(ax, data, x_axis, y_axis=csv_cols.PERFORMANCE, section="baum_welch", N=None, M=None, T=None):
    extracted_data = extract_data(data, x_axis, y_axis, section, N, M, T)
    # ax.set_title(extracted_data[csv_cols.BINARY].iloc[0])
    ax.set_xlabel(x_axis)
    ax.set_ylabel(f"Perf [F/C]", rotation=0)
    ax.yaxis.set_label_coords(-0.05, 1.0)
    ax.grid(axis='x')
    ax.set_xticks(extracted_data[x_axis])
    ax.plot(extracted_data[x_axis], extracted_data[y_axis])
    



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


    data = pd.read_csv('data/dummie_df.csv')

    plt.rcParams.update(plt.rcParamsDefault)
    plt.style.use('ggplot')

    plot_perf_NMT(data, "Performance Plot Variants", N=4, M=6, T=1)

    


    

    
