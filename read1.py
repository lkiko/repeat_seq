# python
# -*- encoding: utf-8 -*-
'''
@File        :read.py
@Time        :2022/02/10 03:30:02
@Author        :charles kiko
@Version        :1.0
@Contact        :charles_kiko@163.com
@Desc        :扫描基因组内的重复序列 python read.py lens fasta out_file
'''
import os
import gc
import sys
import csv
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
from Bio import SeqIO
import multiprocessing # Step I : 导入模块
from multiprocessing import cpu_count#读取CPU核心数用于匹配线程数

def count_str(seq):
    str_dic = {'A':0,'T':0,'C':0,'G':0}
    for i in seq:
        str_dic[i] = str_dic[i] + 1
    if str_dic['G'] + str_dic['C'] == 0:
        return 0
    else:
        return (str_dic['G'] + str_dic['C']) / len(seq)



def find(seq,name,start,end):
    global DNA_dict# 全局变量
    global contig
    global contigs_0
    global dup
    # print(seq)
    file = open('out.out', 'a+')
    file.write(name+'_'+str(start)+'_'+str(end)+'\n')
    file.close()
    name1 = name+str(start)+str(end)
    GC = count_str(seq)
    for contig_ in contig:
        if contig_ not in DNA_dict.keys():
            print(contig_,'not in keys!')
            gc.collect()
            # continue
            return 0
        top = str(DNA_dict[contig_])
        left = str(seq)
        for i in range(len(top) - len(left) + 1):
            if i + len(left) >= len(top):
                continue
            if top[i:i + len(left)] == left:
                if contig_ == name and start == i:
                    continue

                name2 = contig_+str(i)+str(i+len(left))
                if name1 not in dup:
                    dup.append(name1)
                    lt = [contig_,name1,start,end,len(left),GC,seq]
                    file = open(sys.argv[3], 'a+')
                    file.write('\t'.join([str(m) for m in lt])+'\n')
                if name2 not in dup:
                    dup.append(name2)
                    lt = [contig_,name2,i,i+len(left),len(left),GC,seq]
                    file = open(sys.argv[3], 'a+')
                    file.write('\t'.join([str(m) for m in lt])+'\n')
                # print(seq)
                file = open('out.out', 'a+')
                file.write(contig_+'*'+str(i)+'*'+str(i+len(left))+'\n')
                file.close()
                del lt
                gc.collect()
                file.close()
        del top,left
        gc.collect()
    return 0

def find1(seq,name,length):
    global DNA_dict# 全局变量
    global contig
    global contigs_0
    global dup
    # print(seq)
    GC = count_str(seq)
    for contig_ in contig:
        top = list(DNA_dict[contig_])
        left = list(seq)
        df = pd.DataFrame(np.arange((len(top)+1)*(len(left)+1)).reshape(len(left)+1, len(top)+1), index = ['GAP']+left, columns=['GAP']+top)
        df.loc[:,:] = 0
        print(df)
        for i in range(len(top)+1):
            for j in range(len(left)+1):
                if i == 0 or j == 0:
                    continue
                else:

                    if top[i-1] == left[j-1]:
                        df.iloc[j, i] = df.iloc[j-1, i-1] + 1
                    else:
                        continue
        if name == contig_:
            for x in range(len(top)):
                df.iloc[x,x] = 0
        print(df)
        # df.to_csv('forward.csv', index=True, header=True)
        for i in range(length,len(left)+1):
            for j in range(length,len(top)+1):
                if i == len(left) or j == len(top):
                    # print(i,j,df.iloc[i,j])
                    if df.iloc[i,j] >= length:
                        length_ = df.iloc[i,j]
                        ti,tj = j+1-length_,j+1
                        li,lj = i+1-length_,i+1
                        seq = top[ti:tj]
                        GC = count_str(seq)
                        name1 = contig_+str(ti)+str(tj)
                        if name1 not in dup:
                            dup.append(name1)
                            lt = [contig_,name1,ti,tj,length_,GC,"".join(seq)]
                            file = open(sys.argv[3], 'a+')
                            file.write('\t'.join([str(m) for m in lt])+'\n')
                            print(seq)
                            del lt
                            gc.collect()
                            file.close()
                        name2 = name+str(li)+str(lj)
                        if name2 not in dup:
                            dup.append(name2)
                            lt = [name,name2,li,lj,length_,GC,"".join(seq)]
                            file = open(sys.argv[3], 'a+')
                            file.write('\t'.join([str(m) for m in lt])+'\n')
                            print(seq)
                            del lt
                            gc.collect()
                            file.close()
                else:
                    # print('ij',i,j,df.iloc[i,j])
                    if df.iloc[i+1,j+1] >= length:
                        continue
                    else:
                        length_ = df.iloc[i,j]
                        if length_ < length:
                            del length_
                            gc.collect()
                            continue
                        ti,tj = j+1-length_,j+1
                        li,lj = i+1-length_,i+1
                        seq = top[ti:tj]
                        GC = count_str(seq)
                        name1 = contig_+str(ti)+str(tj)
                        if name1 not in dup:
                            dup.append(name1)
                            lt = [contig_,name1,ti,tj,length_,GC,"".join(seq)]
                            file = open(sys.argv[3], 'a+')
                            file.write('\t'.join([str(m) for m in lt])+'\n')
                            print(seq)
                            del lt
                            gc.collect()
                            file.close()
                        name2 = name+str(li)+str(lj)
                        if name2 not in dup:
                            dup.append(name2)
                            lt = [name,name2,li,lj,length_,GC,"".join(seq)]
                            file = open(sys.argv[3], 'a+')
                            file.write('\t'.join([str(m) for m in lt])+'\n')
                            print(seq)
                            del lt
                            gc.collect()
                            file.close()

        gc.collect()
    gc.collect()

if __name__ == '__main__' :#多进程
    length = 5
    dup=[]

    contigs = {}
    contigs_0 = {}

    for i in open(sys.argv[1],'r'):
        if i != '\n':
            lt = i.strip('\n').split()
            contigs[str(lt[0])] = int(lt[1])

    DNA_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[2], "fasta"))# 提取之后直接返回字典
    for key in DNA_dict.keys():
        DNA_dict[key] = str(DNA_dict[key].seq).upper()

    contigs_ = {k: v for k, v in sorted(contigs.items(), key=lambda item: item[1], reverse=True)}

    contig = list(contigs_.keys())
    # print(contig)

    file = open('out.out', 'w')
    file.close()


    file = open(sys.argv[3], 'w')
    file.write(
        '\t'.join(['contig_name','name','start', 'end','seq_length','GC','SEQ'])+'\n')
    file.close()
    pool = multiprocessing.Pool(processes = 8) # Step II : 进程池
    for i in contig:
        if i not in DNA_dict.keys():
            continue
        # pool.apply_async(find1, (DNA_dict[i],i,length,), )  # Step III : 异步（并行）计算
        for j in range(5,len(DNA_dict[i])):
            for m in range(len(DNA_dict[i]) - j + 1):
                contigs_0[str(DNA_dict[i][m:m + j])+str(m)+str(m+j)] = 1
                pool.apply_async(find, (DNA_dict[i][m:m + j],i,m,m+j,), )  # Step III : 异步（并行）计算
    pool.close() # Step IV : 准备结束
    pool.join() # Step IV : 完全结束
