# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 11:16:56 2025

@author: Administrator
"""

import numpy as np 
import threading
import os
import pandas as pd
import itertools
import argparse

def avg(v):
    return sum(v)/float(len(v))

def gini(array):
    array = np.array(array)
    if np.amin(array) < 0:
        array -= np.amin(array)
    array += 0.000001
    array = np.sort(array)
    index = np.arange(1, array.shape[0]+1)
    n = array.shape[0]
    uncorrectedGini = ((np.sum((2 * index - n - 1) * array)) / (n * np.sum(array)))
    correctedGini = ((np.sum((2 * index - n - 1) * array)) / ((n-1) * np.sum(array)))
    return uncorrectedGini, correctedGini

def G4ratio(motif_name,G4list):
    ratio = []
    goodtract = 0
    #num = len(G4list)
    for G4tract in G4list:
        index_list = [i for i in range(len(G4tract)-1) if G4tract[i] < G4tract[i+1]]
        if not index_list:
            continue
        else:
            goodtract += 1
            beforeG = G4tract[index_list[-1]]
            maxG = G4tract[index_list[-1]+1]
            tempratio = float(maxG/(beforeG+maxG))
            ratio.append(tempratio)
    if goodtract >= 4:
        return np.mean(ratio)
    else:
        return 0
               
def GposOnly(motseq, signal):
    gcounts = []
    gcount_temp = []
    for s, y in zip(motseq, signal):
        if s.upper() == 'G':
            gcount_temp.append(y)
        else:
            if len(gcount_temp) > 1:
                gcounts.append(gcount_temp)
            gcount_temp = []
    if gcount_temp and len(gcount_temp) > 1:
        gcounts.append(gcount_temp)
    return [np.array(g) for g in gcounts],gcounts

def load_count_file(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
        data = {}
        for i in range(0, len(lines), 4):
            gene = lines[i].strip()
            seq = lines[i+1].strip().split()
            count = list(map(float, lines[i+2].strip().split()))
            data5 = lines[i+3].strip()
            data[gene] = {'seq': seq, 'count': count, 'data5': data5}
        return data

def process_data(data, rg4_list, output_file, file_name):
    result = []
    with open(rg4_list, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            items = line.strip().split('\t')
            gene_id = items[0]
            name = items[2]
            beg_num = int(items[4]) - 1
            end_num = int(items[5])
            if gene_id in data:
                seq = data[gene_id]['seq'][beg_num:end_num]
                count = data[gene_id]['count'][beg_num:end_num]
                avgCount = avg([c for s, c in zip(seq, count) if s.upper() == 'G'])
                Garray, Glist = GposOnly(seq, count)
                print(Glist)
                Gcluster = list(itertools.chain.from_iterable(Garray))
                giniIndex = gini(Gcluster)[0]
                ratio = G4ratio(name,Glist)
                print(ratio)
                result.append([gene_id] + items[1:] + [avgCount, giniIndex, ratio])
    return result

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-RG4list', type=str, required=True, help='RG4 list file')
    parser.add_argument('-C', '--count', type=str, nargs='+', required=True, help='Count files')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output file')
    parser.add_argument('-t', '--threads', type=int, default=10, help='Number of threads')
    args = parser.parse_args()

    with open(args.RG4list, 'r') as f:
        header = f.readline().strip()
        with open(args.output, 'w') as f_out:
            f_out.write(f'{header}')
            for file_name in args.count:
                f_out.write(f'\tavgCount_{os.path.basename(file_name)}\tgini_{os.path.basename(file_name)}')
            f_out.write('\n')

    results = []
    threads = []
    for file_name in args.count:
        data = load_count_file(file_name)
        t = threading.Thread(target=lambda q, arg1, arg2, arg3, arg4: q.append(process_data(arg1, arg2, arg3, arg4)), args=(results, data, args.RG4list, args.output, os.path.basename(file_name)))
        threads.append(t)
        t.start()
        if len(threads) >= args.threads:
            for t in threads:
                t.join()
            threads = []

    for t in threads:
        t.join()

    df = pd.DataFrame(results[0], columns=['gene_id', 'motif', 'motif_name', 'site', 'beg_num', 'end_num', 'symbol1', 'symbol2', 'symbol3', 'matched_seq', 'avgCount_' + os.path.basename(args.count[0]), 'gini_' + os.path.basename(args.count[0]), 'G4ratio_'+ os.path.basename(args.count[0])])
    for i in range(1, len(results)):
        df_temp = pd.DataFrame(results[i], columns=['gene_id', 'motif', 'motif_name', 'site', 'beg_num', 'end_num', 'symbol1', 'symbol2', 'symbol3', 'matched_seq', 'avgCount_' + os.path.basename(args.count[i]), 'gini_' + os.path.basename(args.count[i]), 'G4ratio_'+ os.path.basename(args.count[i])])
        df = pd.merge(df, df_temp, on=['gene_id', 'motif', 'motif_name', 'site', 'beg_num', 'end_num', 'symbol1', 'symbol2', 'symbol3', 'matched_seq'])

    df.to_csv(args.output, index=False, sep='\t')

if __name__ == "__main__":
    main()
