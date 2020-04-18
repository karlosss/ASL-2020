#!/usr/bin/python3

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import csv_cols


if len(sys.argv) == 1 :
    csv_file = "data/dummie_data.csv"


data = pd.read_csv(csv_file)
print(data.head())


plt.rcParams.update(plt.rcParamsDefault)
plt.style.use('ggplot')

plt.plot(data[csv_cols.PARAM_N], data[csv_cols.PERFORMANCE])
plt.show()