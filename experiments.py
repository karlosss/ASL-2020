#!/usr/bin/env python3

#IMPORTS
import argparse
import datetime
import os
import math
import glob
import subprocess
import json
import plot
import sys
from ast import literal_eval as make_tuple
import python_lib.csv_cols as csv_cols
from python_lib.constants import OUTPUT_DIR, SECTIONS 
#########################

#GLOBALS
CSV_HEADER=(
    f"{csv_cols.BINARY},"
    f"{csv_cols.COMPILER},"
    f"{csv_cols.FLAGS},"
    f"{csv_cols.PARAM_N},"
    f"{csv_cols.PARAM_M},"
    f"{csv_cols.PARAM_T},"
    f"{csv_cols.NUM_ITERATIONS},"
    f"{csv_cols.SECTION},"
    f"{csv_cols.NUM_CYCLES},"
    f"{csv_cols.PERFORMANCE},"
    f"{csv_cols.CACHE_MISS},"
    f"{csv_cols.CACHE_ACCESS},"
    f"{csv_cols.MISS_RATE},"
    f"{csv_cols.NUM_FLOPS},"
    f"{csv_cols.OP_INTENSITY},"
    f"{csv_cols.VARIABLE}\n"
)
T_FACTOR=10

#########################

def read_input():
    parser = argparse.ArgumentParser(
        description=(
            "Run experiment on a specified source file and write to a csv."
        ),
        formatter_class=argparse.RawTextHelpFormatter   
    )
    parser.add_argument(
        "-s", "-sourcefile",
        help='The source code you want to be compiled and the experiments to be .',
        dest='s',default='no_opt'
    )
    parser.add_argument(
        "-f", "--flags",
        help=(
            "The flags to use for compilation (as a single string, and without\n"
            "leading '-', i.e. O0 fno-tree-vectorize), default is -O0."
        ),
        dest="flags", default="O0"
    )
    parser.add_argument(
        "-c","--compiler",
        help=(
           "The compiler to use, i.e. clang or g++; default is g++.\n"
           "Supports absolute path to the compiler binary as well."           
        ),
        dest="compiler",default="g++"
    )
    parser.add_argument(
        "-N", "--hidden-states",
        help=(
            "-N='min, max, step, fix' for range or simply '-N=10' for fixed N.\n"
            "If only three args are provided, the third one is assumed to be step.\n"
            "except for when -e is set, in which case it will be interpreted as fix."
        ),
        dest="N", default="1, 11, 1, 11"
    )
    parser.add_argument(
        "-M", "--observation-alphabet-size", 
        help=(
            "-M='min, max, step, fix' for range or simply '-M=10' for fixed M.\n"
            "If only three args are provided, the third one is assumed to be step.\n"
            "except for when -e is set, in which case it will be interpreted as fix."
            ),
        dest="M", default="1, 11, 1, 11"
    )
    parser.add_argument(
        "-T", "--sequence-length", 
        help=(
            "-T='min, max, step, fix' for range or simply '-T=10' for fixed T.\n"
            "If only three args are provided, the third one is assumed to be step,\n"
            "except for when -e is set, in which case it will be interpreted as fix."
            ),
        dest="T", default="10, 20, 1, 20"
    )
    parser.add_argument(
        '-e', "--exponential", 
        help=(
            "If set, the 'min' and 'max' of the N, M and T ranges are interpreted\n"
            "as base 2 exponents. The third entry dentotes the fixed value used\n"
            "when the parameter does not vary. \n"
            "\n"
            "Example: -e -N='2,11,5' --> N ranges over [2^2, 2^3, ..., 2^10].\n"
            "In the M plot and T plot, a fixed value of N=5 is used."
        ),
        dest="e", action='store_true'
    )
    parser.add_argument(
        "-i", "--iterations",
        help='The number of iterations, 50 if none specified',
        dest= "iters", default=50)
    parser.add_argument(
        "-p", "--show-plot", 
        help=(
            "If set, the plot will pop up at the end (the .png file will be created \n"
            "either way)."
        ),
        dest="p", action='store_true'
    )
    return parser.parse_args()


#takes the name of a binary and a set of flags, creates a csv to write experiment results to
#returns the path to the csv
def create_csv(binary, compiler, flags, csv_header):

    #make string of flags
    flag_string= ('_').join(flags.split(' '))

    timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')

    #create name for csv
    dir_name = str("{0}@{1}_{2}@{3}".format(binary, compiler.split("/")[-1], flag_string, timestamp))
    csv_name= "report.csv"

    #get absolute path for csv
    path = os.path.dirname(__file__)

    path = os.path.join(path, OUTPUT_DIR, dir_name)

    if not os.path.exists(path):
        os.makedirs(path)
    csv_path = os.path.join(path,csv_name)

    #create file
    with open(csv_path,'w+') as csv:
        csv.write(csv_header)

    return csv_path, path


#takes a list of lines with comma separated values and writes it to the output csv
def append_csv(csv_path, lines, binary, flags, n, m, t, iter, compiler, variable):
    csv_lines = ['{0},{1},{2},{3},{4},{5},{6},{7},{8}\n'.format(binary, compiler.split("/")[-1], flags, n, m, t, iter, dat, variable) for dat in lines]
    with open(csv_path,'a') as csv:
        for line in csv_lines:
            csv.write(line)


def get_flops_from_binary(binary, compiler, flags, n, m, t, iters):
    #call binary directly to get flop count

    #derive path for binary
    flags_arr = flags.split(' ')
    flags_arr = ['-{0}'.format(f) for f in flags_arr]
    flag_string = '{1}_{0}'.format('_'.join(flags_arr), compiler)
    exec_path = './bin/{0}/{1} flops <<< \"{2} {3} {4} {5}\"'.format(flag_string,binary, n, m, t, iters)

    stream = os.popen(exec_path)
    exp_out = stream.read()

    return int(exp_out)

#Takes the path to PAPI logs, parses jsons and return a list of csv lines; one line for every region
def get_data(json_path, binary, compiler, flags, n, m, t, iters):
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
                    flops = str(int(scalar) + 4*int(vectorized))
                    performance = str(int(flops)/int(cycles))
                    cache_miss = reg['PAPI_L3_TCM']
                    op_intensity = int(flops)/(64*float(cache_miss))
                    cache_access = reg['PAPI_L3_TCA']
                    cache_miss_rate = str(float(cache_miss)/float(cache_access))
                    out.append('{0},{1},{2},{3},{4},{5},{6},{7}'.format(name,cycles,performance, 
                        cache_miss, cache_access, cache_miss_rate,flops,op_intensity))

                    i=i+1

            except KeyError:    
                print("Your machine does not support PAPI FLOP counters, proceeding with less acurate count.")
                
                name = SECTIONS[0]
                reg = regions[0][name]
                cycles = reg['cycles']
                flops = get_flops_from_binary(binary, compiler, flags, n, m, t, iters)
                performance=flops/int(cycles)
                cache_miss = reg['PAPI_L3_TCM']
                cache_access = reg['PAPI_L3_TCA']
                cache_miss_rate = str(float(cache_miss)/float(cache_access))
                out.append('{0},{1},{2},{3},{4},{5}'.format(name,cycles,performance, 
                    cache_miss, cache_access, cache_miss_rate))

    files = glob.glob(f"{json_path}/*")
    for f in files:
        os.remove(f)
    return out

# display realtime output of a process and return its exit code when it is finished.
def get_realtime_output(process):
    while True:
        if process.poll() is not None:
            return process.poll()
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip().decode("utf8"))


#compiles all binaries present with provided flags
def compile_all(compiler, flags):
    print("Compiling with {0}".format(flags))
    #append "-" to flags
    flags = ['-{0}'.format(f) for f in flags.split(' ')]

    args = flags
    args.insert(0,compiler)
    args.insert(0,'./compile.sh')
    #call compile.sh
    comp = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
    retcode = get_realtime_output(comp)
    if retcode != 0:
        print("Some files could not be compiled or validated!", file=sys.stderr)
        exit(1)


#calls run.sh binary with n, m, iters and t and writes a line with experiment results to file at csv_path
def run_experiment(binary, compiler, flags, n, m, iters, t, csv_path, variable):

        #derive path for binary 
        flag_arr = flags.split(' ')
        flag_arr = ['-{0}'.format(f) for f in flag_arr]
        args = flag_arr
        args.insert(0,compiler)
        args.insert(0,binary)
        args.insert(0,'./run.sh')

        #run experiment
        exp = subprocess.Popen(args, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT)
        exp_out = exp.communicate(input=('{0} {1} {2} {3}'.format(n,m,iters,t)).encode('utf-8'))

        #retrieve data from PAPI json and form a comma separated value line
        data = get_data('logs/papi_hl_output', binary, compiler, flags, n, m, t, iters)

        #append the line to csv
        append_csv(csv_path, data, binary, flags, n, m, t, iters, compiler.split("/")[-1], variable)


def main(binary, N_iter, M_iter, T_iter, iters, N_fix, M_fix, T_Fix, flags, compiler):
    #create a csv
    csv_path, dir_path = create_csv(binary, compiler, flags, CSV_HEADER)
    print("CSV Created")


    #compilation
    compile_all(compiler, flags)
    print("Files compiled")
                
    #run all experiments

    print("=== Computing N plot data ===")
    for n in N_iter:
        print("running for N: {0}, M: {1}, T: {2}".format(n,M_fix,T_Fix))
        run_experiment(binary, compiler, flags, n, M_fix, iters, T_Fix, csv_path, 0)
    
    print("\n=== Computing M plot data ===")
    for m in M_iter:
        print("running for N: {0}, M: {1}, T: {2}".format(N_fix,m,T_Fix))
        run_experiment(binary, compiler, flags, N_fix, m, iters, T_Fix, csv_path, 1)
        
    print("\n=== Computing T plot data ===")
    for t in T_iter:
        print("running for N: {0}, M: {1}, T: {2}".format(N_fix,M_fix,t))
        run_experiment(binary, compiler, flags, N_fix, M_fix, iters, t, csv_path, 2)
    
    print("All experiments done")
    print("Find your output in: {0}".format(csv_path))

    print(f"Plotting summary...")
    plot.multiplot_NP_MP_TP_S(csv_path, dir_path, N_fix, M_fix, T_Fix, args.p)
    plot.multiplot_NP_MP_TP_Cache(csv_path, dir_path, N_fix, M_fix, T_Fix, args.p)
    plot.plot_roofline(csv_path, dir_path, N_fix, M_fix, T_Fix, args.p)
    plot.multiplot_runtime(csv_path, dir_path, N_fix, M_fix, T_Fix, args.p)
    

def parse_tuple(arg_name, arg_string, exp):
    t = make_tuple(arg_string)
    # In case of scalar argument return a tuple representing a 1-number range
    if type(t) == int:
        if not exp: return (t, t+1, 1, t)
        else: return (t, t+1, t)

    elif type(t) == tuple:
        if not exp:
            if len(t) == 4:
                    _min, _max, _step, _fix = t
            elif len(t) == 3:
                    _min, _max, _step = t
                    _fix = _max
                    print(
                        f"Warning for argumet '{arg_name}': No fixed value provided. "
                        f"Using max value."
                    )
            elif len(t) == 2: 
                (_min, _max), _step = t, 1
                _fix = _max
                print(
                    f"Warning for argument '{arg_name}': No step size p and no fixed value provided. "
                    f"Using stepsize 1 and fixed value max."
                )
            else: 
                print(f"Error in argument '{arg_name}': Tuple {t} too long.")
                return None

            if _min > _max:
                print(f"Error in argument '{arg_name}': Min={_min} is larger than max={_max}.")
                return None

            if _step < 0:
                print(f"Error in argument '{arg_name}': Negative stepsize.")
                return None
            return (_min, _max+1, _step, _fix)
        else:
            if len(t) == 4:
                _min, _max, _, _fix = t
                print(
                    f"Warning for argument '{arg_name}': Stepsize will be ignored because"
                    f"'-e' flag is set."
                )
            elif len(t) == 3:
                _min, _max, _fix = t

            elif len(t) == 2:
                _min, _max = t
                _fix = _max
            else: 
                print(f"Error in argument '{arg_name}': Tuple {t} too long.")
                return None


            if _min > _max:
                print(f"Error in argument '{arg_name}': Min={_min} is larger than max={_max}.")
                return None
            return (_min, _max+1, _fix)



if __name__=='__main__':

    args = read_input()

    N_tuple = parse_tuple("-N", args.N, args.e)    
    M_tuple = parse_tuple("-M", args.M, args.e)    
    T_tuple = parse_tuple("-T", args.T, args.e)

    if None in [N_tuple, M_tuple, T_tuple]:
        print(f"Abort")
        exit(1)   

    if len(N_tuple) == 4:
        n_min, n_max, _, n_fix = N_tuple
    else:
        n_min, n_max, n_fix = N_tuple

    if len(M_tuple) == 4:        
        m_min, m_max, _, m_fix = M_tuple
    else:
        m_min, m_max, m_fix = M_tuple

    if len(T_tuple) == 4:
        t_min, t_max, t_step, t_fix = T_tuple
    else:
        t_min, t_max, t_fix = T_tuple

    # t_min = T_FACTOR * (max(m_max, n_max))
    if args.e:
        T_FACTOR = 2
        lower_bound = (2 ** max(n_fix, m_fix)) * T_FACTOR
        t_min = max(t_min, math.ceil(math.log2(lower_bound)))
        t_max = max(t_min+1, t_max)
        lower_bound = (2 ** (max(n_max, m_max) - 1)) * T_FACTOR
        t_fix = max(t_fix, math.ceil(math.log2(lower_bound)))
        n_fix = 2**n_fix
        m_fix = 2**m_fix
        t_fix = 2**t_fix
        
    else:
        t_min = max(t_min, T_FACTOR * max(n_fix, m_fix))
        t_max = max(t_min+1, t_max)
        t_fix = max(t_fix, T_FACTOR * (max(n_max, m_max) - 1))
    
    # Update T_tuple with adusted parameters.
    if args.e:
        t_step = 1
    T_tuple = tuple([t_min, t_max, t_step, t_fix])
    def create_iterator(t, exp):
        if not exp:
            if len(t) == 4:
                t = (t[0],t[1],t[2])
            return list(range(*t))
        else: 
            t = (t[0],t[1])
            return [2**i for i in range(*t)]
    

    N_iter = create_iterator(N_tuple, args.e)    
    M_iter = create_iterator(M_tuple, args.e)    
    T_iter = create_iterator(T_tuple, args.e)
    
    main(
        args.s, 
        N_iter, M_iter, T_iter,
        int(args.iters), 
        n_fix, m_fix, t_fix,
        args.flags,
        args.compiler
    )
