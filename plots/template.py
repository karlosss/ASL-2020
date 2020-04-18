
#%%
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import os
import csv

def get_data(path, use_rdtsc=True):
    # with open(os.path.dirname(__file__) + '/../data/2b_1.txt', mode='r') as csv_file:
    with open(path, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        print(csv_reader.fieldnames)
        fields = csv_reader.fieldnames
        triplets = [[row[f] for f in fields] for row in csv_reader]
        lists = list(map(list, zip(*triplets)))
        
        ns = list(map(int, lists[0]))
        rdtsc = list(map(float, lists[1]))
        c_clocks = list(map(float, lists[2]))
    return ns, (rdtsc if use_rdtsc else c_clocks)


def plot_data(ns, data, title, y_start, y_end, id):
    

    plt.subplot(1, 3, id)
    plt.plot(ns, data, '-o')
    ax = plt.gca()
    ax.set_title(title, pad=20)
    ax.yaxis.set_label_coords(-0.05, 1.06)
    # ax.yaxis.set_major_formatter(StrMethodFormatter('{x:,.1f}'))
    ax.set_ylabel('[F/C]', rotation=0)
    ax.set_xlabel('Input size')
    ax.tick_params(axis ='x', width=1)
    plt.xticks(range(0, 4001, 500))
    plt.ylim((y_start, y_end))     
    plt.grid(axis='x')

plt.rcParams.update(plt.rcParamsDefault)
plt.style.use('ggplot')


plot_data(*get_data('data/2c_1.csv', use_rdtsc=True), '-O0, no optimizations', 0.0, 0.3, 1)
plot_data(*get_data('data/2c_2.csv', use_rdtsc=True), '-O3, without SIMD', 0.0, 0.6, 2)
plot_data(*get_data('data/2c_3.csv', use_rdtsc=True), '-O3, with SIMD', 0.0, 2.0, 3)

# plot_data(*get_data('../data/2c_3.csv', use_rdtsc=True), 'Test', 0.0, 2.0, 3)

plt.show()  
    

# %%
