#!/usr/bin/python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import csv_cols
import random

GENERATE = False

def extract_data(data, x_axis, section="baum_welch", N=None, M=None, T=None):
    if  x_axis == csv_cols.PARAM_N: 
        fixed_0, fixed_1 = data[csv_cols.PARAM_M] == M, data[csv_cols.PARAM_T] == T
    elif x_axis == csv_cols.PARAM_M: 
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_T] == T
    elif x_axis == csv_cols.PARAM_T:
        fixed_0, fixed_1 = data[csv_cols.PARAM_N] == N, data[csv_cols.PARAM_M] == M
    
    section_constraint = data[csv_cols.SECTION] == section
    return  data[fixed_0 & fixed_1 & section_constraint]

    # plt.plot(extracted_data[x_axis], extracted_data[csv_cols.PERFORMANCE])
    # plt.show()    


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

    fixed_M = data[csv_cols.PARAM_M] == 6
    fixed_T = data[csv_cols.PARAM_T] == 1
    forward_vars = data[csv_cols.SECTION] == "forward_vars"

    

    
