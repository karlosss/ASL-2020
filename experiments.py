#! /usr/bin/env python3

#IMPORTS
import argparse
import datetime
import os
import glob
import subprocess
import json
import plot
import python_lib.csv_cols as csv_cols
from python_lib.constants import OUTPUT_DIR, SECTIONS 
#########################

#GLOBALS
CSV_HEADER=(
    f"{csv_cols.BINARY},"
    f"{csv_cols.FLAGS},"
    f"{csv_cols.PARAM_N},"
    f"{csv_cols.PARAM_M},"
    f"{csv_cols.PARAM_T},"
    f"{csv_cols.NUM_ITERATIONS},"
    f"{csv_cols.SECTION},"
    f"{csv_cols.NUM_CYCLES},"
    f"{csv_cols.PERFORMANCE}\n"
)
T_FACTOR=100

#########################

def read_input():
    parser = argparse.ArgumentParser(description='Run experiment on a specified binary and write to a csv, returns the absolute path to the csv')

    parser.add_argument('-b','-binary',help='The binary you want the experiments to be run on',dest='bin',default='no_opt')
    parser.add_argument('-n_min',help='The minimal value for N, 20 if none specified',dest='n_min',default=20)
    parser.add_argument('-n_max',help='The maximum value for N, 21 if none specified',dest='n_max',default=21)
    parser.add_argument('-n_step',help='The stepping value for N, 10 if none specified',dest='n_step',default=10)
    parser.add_argument('-m_min',help='The minimum value for M, 20 if none specified',dest='m_min',default=20)
    parser.add_argument('-m_max',help='The maximum value for M, 21 if none specified',dest='m_max',default=21)
    parser.add_argument('-m_step',help='The stepping value for M, 10 if none specified',dest='m_step',default=10)
    parser.add_argument('-i',help='The number of iterations, 50 if none specified',dest= 'iters', default=50)
    parser.add_argument('-t_min',help='The minimum value of T (min is 10000*min(N_max,M_max))',dest='t_max',default=0)
    parser.add_argument('-t_max',help='The maximal value of T (min is 10000*min(N_max,M_max))',dest='t_max',default=0)
    parser.add_argument('-t_step',help='The stepping value of T, default is 10',default=10)
    parser.add_argument('-f',help='The flags to use for compilation (as a single string, and without leading "-", i.e. O0 fno-tree-vectorize), default is -O0',dest='flags',default='O0')
    return parser.parse_args()


#takes the name of a binary and a set of flags, creates a csv to write experiment results to
#returns the path to the csv
def create_csv(binary, flags, csv_header):

    #make string of flags
    flag_string= ('_').join(flags.split(' '))

    #create name for csv
    csv_name=str("{0}%{1}.csv".format(binary, flag_string))

    #get absolute path for csv
    path = os.path.dirname(__file__)
    path = os.path.join(path, OUTPUT_DIR)
    if not os.path.exists(path):
        os.makedirs(path)
    csv_path = os.path.join(path,csv_name)

    #create file
    with open(csv_path,'w+') as csv:
        csv.write(csv_header)

    return csv_path


#takes a list of lines with comma separated values and writes it to the output csv
def append_csv(csv_path, lines, binary, flags, n, m, t, iter):
    csv_lines = ['{0},{1},{2},{3},{4},{5},{6}\n'.format(binary, flags, n, m, t, iter, dat) for dat in lines]
    with open(csv_path,'a') as csv:
        for line in csv_lines:
            csv.write(line)


def get_flops_from_binary(binary, n, m, t, iters):
    #call binary directly to get flop count
    stream = os.popen('./bin/{0} flops <<< \"{1} {2} {3} {4}\"'.format(binary, n, m, t, iters))
    exp_out = stream.read()

    return int(exp_out)

#Takes the path to PAPI logs, parses jsons and return a list of csv lines; one line for every region
def get_data(json_path, binary, n, m, t, iters):
    out=[]
    for file in os.scandir(json_path):
        with open(file) as f:
            
            perf = json.load(f)
            regions = perf['threads'][0]['regions']

            try:
                i=0
                for r in regions:

                    name = SECTIONS[i]
                    reg = r[name]
                    cycles = reg['cycles']  
                    scalar = reg['FP_ARITH:SCALAR_DOUBLE']
                    vectorized = reg['FP_ARITH:256B_PACKED_DOUBLE']
                    performance = str(int(scalar + 4*vectorized)/int(cycles))
                    out.append('{0},{1},{2}'.format(name,cycles,performance))

                    i=i+1

            except KeyError:    
                print("Your machine does not support PAPI FLOP counters, proceeding with less acurate count")
                
                name = SECTIONS[0]
                reg = regions[0][name]
                cycles = reg['cycles']
                flops = get_flops_from_binary(binary, n, m, t, iters)
                performance=flops/int(cycles)
                out.append('{0},{1},{2}'.format(name,cycles,performance))

    files = glob.glob(f"{json_path}/*")
    for f in files:
        os.remove(f)
    return out



#compiles all binaries present with provided flags
def compile_all(flags):
    print("Compiling with {0}".format(flags))
    #append "-" to flags
    flags = ['-{0}'.format(f) for f in flags.split(' ')]
    #call compile.sh
    comp = subprocess.Popen(['./compile.sh'], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)



#calls run.sh binary with n, m, iters and t and writes a line with experiment results to file at csv_path
def run_experiment(binary, flags, n, m, iters, t, csv_path):

        #run experiment
        exp = subprocess.Popen(['./run.sh',binary], stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        exp_out = exp.communicate(input=('{0} {1} {2} {3}'.format(n,m,iters,t)).encode('utf-8'))

        #retrieve data from PAPI json and form a comma separated value line
        data = get_data('logs/papi_hl_output', binary, n, m, t, iters)

        #append the line to csv
        append_csv(csv_path, data, binary, flags, n, m, t, iters)



def main(binary, n_min, n_max, n_step, m_min, m_max, m_step, iters, t_min, t_max, t_step, flags):

    #create a csv
    csv_path = create_csv(binary, flags, CSV_HEADER)
    print("CSV Created")

    #compilation
    compile_all(flags)
    print("Files compiled")

    #run all experiments
    for n in range(n_min, n_max, n_step):
        for m in range(m_min, m_max, m_step):
            for t in range(t_min, t_max, t_step):

                print("running for N: {0}, M: {1}, T: {2}".format(n,m,t))

                run_experiment(binary, flags, n, m, iters, t, csv_path)
    
    print("All experiments done")
    print("Find your output in: {0}".format(csv_path))

    print(f"Plotting summary...")
    plot.multiplot_NP_MP_S(csv_path)
    


if __name__=='__main__':

    args = read_input()
    t_min = T_FACTOR*int(args.m_max) if int(args.m_max)<int(args.n_max) else T_FACTOR*int(args.n_max)
    t_max = int(args.t_max) if int(args.t_max) > int(t_min) else int(t_min)+1

    main(
        args.bin, 
        int(args.n_min), 
        int(args.n_max), 
        int(args.n_step), 
        int(args.m_min), 
        int(args.m_max), 
        int(args.m_step), 
        int(args.iters),
        int(t_min), 
        int(t_max), 
        int(args.t_step), 
        args.flags)
