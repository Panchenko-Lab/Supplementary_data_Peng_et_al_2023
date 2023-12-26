#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""      
Created on Fri Jul 12 13:10:51 2019

@author: pengy10
"""
import pandas as pd 
import os 
import re
import csv
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("nucleosome_position", type=str,
                    help="The file contains the nucleosom dyad positions")
parser.add_argument("output_name", type=str,
                    help="output file name")

args = parser.parse_args()
nucleosome_position = args.nucleosome_position
output_name = args.output_name


# Reading hg 19 FASTA sequence of each chromosome
hg_files = ["chr1.fa", "chr2.fa", "chr3.fa", "chr4.fa", "chr5.fa", "chr6.fa", "chr7.fa", "chr8.fa", "chr9.fa",
             "chr10.fa", "chr11.fa", "chr12.fa", "chr13.fa" ,"chr14.fa" ,"chr15.fa" ,"chr16.fa", "chr17.fa", 
             "chr18.fa", "chr19.fa", "chr20.fa",  "chr21.fa", "chr22.fa", "chrX.fa", "chrY.fa"]

chr_seq = {}


for chr in hg_files:
    with open("hg19_GRCh37_UCSC/" + chr) as file:
        data = file.readlines()
        for i in range(0,len(data)):
            data[i] = data[i].strip("\n")
        seq = "".join(data)
        chr_seq[chr[0:-3]] = seq.upper()
        
#reading the nucleosome position 

nuc_pos = pd.read_csv(nucleosome_position,  sep = "\t", header=None)


def nuc_flank_CG(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1: position + flank_len +2]
    for m in re.finditer("CG", flank_seq):
        position = m.span()[0] - flank_len -1
        CG_count[position] += 1
        CG_count[position+1] += 1
# WW = AT, TA , AA , TT
# SS = CG, GC, CC, GG
# YY = CC, TT, CT, TC
# RR = AA, GG, AG, GA
def nuc_flank_WW(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1: position + flank_len +2]
    for m in re.finditer("AT", flank_seq):
        position = m.span()[0] - flank_len -1
        WW_count[position] += 1
        WW_count[position+1] += 1
    for m in re.finditer("TA", flank_seq):
        position = m.span()[0] - flank_len -1
        WW_count[position] += 1
        WW_count[position+1] += 1
    for m in re.finditer("AA", flank_seq):
        position = m.span()[0] - flank_len -1
        WW_count[position] += 1
        WW_count[position+1] += 1
    for m in re.finditer("TT", flank_seq):
        position = m.span()[0] - flank_len -1
        WW_count[position] += 1 
        WW_count[position+1] += 1
        
        
def nuc_flank_SS(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1: position + flank_len +2]
    for m in re.finditer("CG", flank_seq):
        position = m.span()[0] - flank_len -1
        SS_count[position] += 1
        SS_count[position+1] += 1
    for m in re.finditer("GC", flank_seq):
        position = m.span()[0] - flank_len -1
        SS_count[position] += 1
        SS_count[position+1] += 1
    for m in re.finditer("GG", flank_seq):
        position = m.span()[0] - flank_len -1
        SS_count[position] += 1
        SS_count[position+1] += 1
    for m in re.finditer("CC", flank_seq):
        position = m.span()[0] - flank_len -1
        SS_count[position] += 1 
        SS_count[position+1] += 1
        
def nuc_flank_YY(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1 : position + flank_len +2]
    for m in re.finditer("CC", flank_seq):
        position = m.span()[0] - flank_len -1
        YY_count[position] += 1
        YY_count[position+1] += 1
    for m in re.finditer("TT", flank_seq):
        position = m.span()[0] - flank_len -1
        YY_count[position] += 1
        YY_count[position+1] += 1
    for m in re.finditer("CT", flank_seq):
        position = m.span()[0] - flank_len -1
        YY_count[position] += 1
        YY_count[position+1] += 1
    for m in re.finditer("TC", flank_seq):
        position = m.span()[0] - flank_len -1
        YY_count[position] += 1 
        YY_count[position+1] += 1

def nuc_flank_RR(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1 : position + flank_len +2]
    for m in re.finditer("AA", flank_seq):
        position = m.span()[0] - flank_len -1
        RR_count[position] += 1
        RR_count[position+1] += 1
    for m in re.finditer("GG", flank_seq):
        position = m.span()[0] - flank_len -1
        RR_count[position] += 1
        RR_count[position+1] += 1
    for m in re.finditer("AG", flank_seq):
        position = m.span()[0] - flank_len -1
        RR_count[position] += 1
        RR_count[position+1] += 1
    for m in re.finditer("GA", flank_seq):
        position = m.span()[0] - flank_len -1
        RR_count[position] += 1 
        RR_count[position+1] += 1
#    print(CG_count)
        
def nuc_flank_GC(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1: position + flank_len + 2]
    for m in re.finditer("C", flank_seq):
        position = m.span()[0] - flank_len -1
        GC_count[position] += 1
        
    for m in re.finditer("G", flank_seq):
        position = m.span()[0] - flank_len -1
        GC_count[position] += 1
        
def nuc_flank_AT(chr_name, position, flank_len):
    flank_seq = chr_seq[chr_name][position - flank_len -1: position + flank_len + 2]
    for m in re.finditer("A", flank_seq):
        position = m.span()[0] - flank_len -1
        AT_count[position] += 1
        
    for m in re.finditer("T", flank_seq):
        position = m.span()[0] - flank_len -1
        AT_count[position] += 1
    

#CG_count = { i : 0 for i in range(-501,502) }
#GC_count = { i : 0 for i in range(-501,502) }
#AT_count = { i : 0 for i in range(-501,502) }
WW_count = { i : 0 for i in range(-501,502) }
SS_count = { i : 0 for i in range(-501,502) }
YY_count = { i : 0 for i in range(-501,502) }
RR_count = { i : 0 for i in range(-501,502) }

#for i in range(0,20000):
for i in range(0,len(nuc_pos)):
    if nuc_pos[0][i] + ".fa" in hg_files: 
#        nuc_flank_CG(nuc_pos[0][i], int(nuc_pos[1][i]),500)
#        nuc_flank_GC(nuc_pos[0][i], int(nuc_pos[1][i]),500)
#        nuc_flank_AT(nuc_pos[0][i], int(nuc_pos[1][i]),500)
        nuc_flank_WW(nuc_pos[0][i], int(nuc_pos[1][i]),500)
        nuc_flank_SS(nuc_pos[0][i], int(nuc_pos[1][i]),500)
        nuc_flank_YY(nuc_pos[0][i], int(nuc_pos[1][i]),500)
        nuc_flank_RR(nuc_pos[0][i], int(nuc_pos[1][i]),500)

#        print(i)
## mapping of the sequence into 200bp windows in the genome
#for chr in hg_files:
#    for pos in range(250,len(chr_seq[chr[0:-3]]),500):
 #   for pos in range(10100,20000,200):    
#        nuc_flank_CG(chr[0:-3], pos,500)
#        nuc_flank_GC(chr[0:-3], pos,500)
#        nuc_flank_AT(chr[0:-3], pos,500)
#        nuc_flank_WW(chr[0:-3], pos,500)
#        nuc_flank_SS(chr[0:-3], pos,500)
#        print(chr,pos)        


## remove -501 and 501 in tge output files
w = csv.writer(open("./patterns_different_celllines/WW_" + output_name + ".csv", "w"))
for key, val in WW_count.items():
    w.writerow([key, val])
    
w = csv.writer(open("./patterns_different_celllines/SS_" + output_name + ".csv", "w"))
for key, val in SS_count.items():
    w.writerow([key, val])

#w = csv.writer(open("./patterns_different_celllines/CG_" + output_name + ".csv", "w"))
#for key, val in CG_count.items():
#    w.writerow([key, val])
    
#w = csv.writer(open("./patterns_different_celllines/GC_" + output_name + ".csv", "w"))
#for key, val in GC_count.items():
#    w.writerow([key, val])
    
#w = csv.writer(open("./patterns_different_celllines/AT_" + output_name + ".csv", "w"))
#for key, val in AT_count.items():
#    w.writerow([key, val])
        
w = csv.writer(open("./patterns_different_celllines/YY_" + output_name + ".csv", "w"))
for key, val in YY_count.items():
    w.writerow([key, val])
    
w = csv.writer(open("./patterns_different_celllines/RR_" + output_name + ".csv", "w"))
for key, val in RR_count.items():
    w.writerow([key, val])
