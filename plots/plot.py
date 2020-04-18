import os
import pandas as pd
import matplotlib.pyplot as plt

def extract_columns(data, x, y):
    return data


data = pd.read_csv("data/dummie_data.csv")
print(data.head())


plt.rcParams.update(plt.rcParamsDefault)
plt.style.use('ggplot')

plt.plot(data["N"], data["PERFORMANCE"])
plt.show()