#!/usr/bin/python
# This code is compatable with python3
# Developed by Ram Sai Nanduri
# Github: https://github.com/ramsainanduri
# Script to extract vcf records that falls in the given bed region

############################### Required Packages ###########################################

import os
import sys
import csv
import argparse
import warnings
import subprocess
import pandas as pd
from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

############################### Usage Block ###########################################
parser = argparse.ArgumentParser(description="Extract vcf records that falls in the given bed region")
parser.add_argument("-i", "--vcf", metavar="vcf", required=True, help="Input vcf file", type=str)
parser.add_argument("-b", "--bed", metavar="bed", required=True, help="Input Bed file", type=str)
parser.add_argument("-o", "--output", metavar="out_vcf", required=True, help="Output vcf file", type=str)
parser.add_argument("-r", "--reverse", metavar="reverse_output", default="N", help="Y, prints vcf records that doesnot fall in the given bed region, default N", type=str)
parser.add_argument("-v", "--version", help="Program's version", action='version', version='%(prog)s 1.0')
args = parser.parse_args()

########################## General Variable Defining ######################################

cwd = os.getcwd()
folder = os.path.basename(cwd)
vcf_ori = cwd+'/'+args.vcf
bed_ori = cwd+'/'+args.bed
final_vcf_out = cwd+'/'+args.output
reverse = args.reverse
vcf_head_name = vcf_ori+'.temp.header.vcf'
vcf_merged = cwd+'/All_merged.temp'
vcf_sorted = cwd+'/All_merged.sorted.temp'

########################### Deleteing if exsiting ##############################

if os.path.exists(final_vcf_out+'*'):
    os.remove(final_vcf_out+'*')
if os.path.exists(cwd+'/chr*.temp.vcf'):
    os.remove(cwd+'/chr*.temp.vcf')
if os.path.exists(cwd+'/chr*.temp.vcf'):
    os.remove(cwd+'/chr*.temp.vcf')
if os.path.exists(vcf_merged+'*'):
    os.remove(vcf_merged+'*')
if os.path.exists(vcf_sorted+'*'):
    os.remove(vcf_sorted+'*')
if os.path.exists(vcf_head_name):
    os.remove(vcf_head_name)

########################## Data Reading #####################################

vcf_raw = pd.read_csv(vcf_ori, sep='\t', header=None)
bed_raw = pd.read_csv(bed_ori, sep='\t', header=None)
vcf_head = vcf_raw[vcf_raw[0].str.contains('#')]
vcf_head.to_csv(vcf_head_name,sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
vcf_nohead = vcf_raw[~vcf_raw[0].str.contains('#')]
cols = len(vcf_nohead.columns)
chrs_bed = list(bed_raw[0].unique())
chrs_vcf = list(vcf_nohead[0].unique())

######################### Data Processing ##################################

for b in chrs_bed:
	bed_name = b+".temp.bed"
	vcf_name = b+".temp"
	chrs_bed_df = bed_raw[bed_raw[0] == b]
	bed_pos_srt = list(chrs_bed_df[1])
	bed_pos_end = list(chrs_bed_df[2])
	chrs_bed_df.to_csv(bed_name,sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
	vcf_nohead_chr = vcf_nohead[(vcf_nohead[0] == b)]
	vcf_nohead_passed = vcf_nohead_chr[(vcf_nohead_chr[6] == "PASS")]
	vcf_pos_fil_passed = list(vcf_nohead_passed[1])
	vcf_pos_fil_all = list(vcf_nohead_chr[1])
	
	vcf_pos_final_passed = []
	vcf_pos_final_all = []

	for vcf_pos in vcf_pos_fil_passed:
		for srt, end in zip(bed_pos_srt, bed_pos_end):
			if int(vcf_pos) >= int(srt) and int(vcf_pos) <= int(end):
				vcf_pos_final_passed.append(vcf_pos)

	for vcf_pos in vcf_pos_fil_all:
		for srt, end in zip(bed_pos_srt, bed_pos_end):
			if int(vcf_pos) >= int(srt) and int(vcf_pos) <= int(end):
				vcf_pos_final_all.append(vcf_pos)

	vcf_pos_final_out_passed = []
	vcf_pos_final_out_all = []

	if reverse == "N":
		for x in vcf_pos_final_passed:
			if x not in vcf_pos_final_out_passed:
				vcf_pos_final_out_passed.append(x)

		for x in vcf_pos_final_all:
			if x not in vcf_pos_final_out_all:
				vcf_pos_final_out_all.append(x)

		vcf_chr_final_passed = vcf_nohead_passed[vcf_nohead_passed[1].isin(vcf_pos_final_out_passed)]
		vcf_chr_final_all = vcf_nohead_chr[vcf_nohead_chr[1].isin(vcf_pos_final_out_all)]
		vcf_chr_final_passed.to_csv(vcf_name+'.passed.vcf',sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
		vcf_chr_final_all.to_csv(vcf_name+'.all.vcf',sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)

	elif reverse == "Y":

		for x in vcf_pos_fil_all:
			if x not in vcf_pos_final_all:
				vcf_pos_final_out_all.append(x)

		for x in vcf_pos_fil_passed:
			if x not in vcf_pos_final_passed:
				vcf_pos_final_out_passed.append(x)

		vcf_chr_final_passed = vcf_nohead_passed[vcf_nohead_passed[1].isin(vcf_pos_final_out_passed)]
		vcf_chr_final_all = vcf_nohead_chr[vcf_nohead_chr[1].isin(vcf_pos_final_out_all)]
		vcf_chr_final_passed.to_csv(vcf_name+'.reverse.passed.vcf',sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)
		vcf_chr_final_all.to_csv(vcf_name+'.reverse.all.vcf',sep='\t', index=False, header=False, quoting=csv.QUOTE_NONE)


######################### Writing Output Results ####################################
if reverse == "N":
	subprocess.call('cat '+cwd+'/chr*.temp.passed.vcf > '+vcf_merged+'.passed.vcf', shell=True)
	subprocess.call('cat '+cwd+'/chr*.temp.all.vcf > '+vcf_merged+'.all.vcf', shell=True)
	subprocess.call('sort -k1,1 -k2,2 -V -s '+vcf_merged+'.passed.vcf > '+vcf_sorted+'.passed.vcf', shell=True)
	subprocess.call('cat '+vcf_head_name+' '+vcf_sorted+'.passed.vcf > '+final_vcf_out+'.passed.vcf', shell=True)
	subprocess.call('sort -k1,1 -k2,2 -V -s '+vcf_merged+'.all.vcf > '+vcf_sorted+'.all.vcf', shell=True)
	subprocess.call('cat '+vcf_head_name+' '+vcf_sorted+'.all.vcf > '+final_vcf_out+'.all.vcf', shell=True)
elif reverse == "Y":
	subprocess.call('cat '+cwd+'/chr*.temp.reverse.passed.vcf > '+vcf_merged+'.reverse.passed.vcf', shell=True)
	subprocess.call('cat '+cwd+'/chr*.temp.reverse.all.vcf > '+vcf_merged+'.reverse.all.vcf', shell=True)
	subprocess.call('sort -k1,1 -k2,2 -V -s '+vcf_merged+'.reverse.passed.vcf > '+vcf_sorted+'.reverse.passed.vcf', shell=True)
	subprocess.call('cat '+vcf_head_name+' '+vcf_sorted+'.reverse.passed.vcf > '+final_vcf_out+'.reverse.passed.vcf', shell=True)
	subprocess.call('sort -k1,1 -k2,2 -V -s '+vcf_merged+'.reverse.all.vcf > '+vcf_sorted+'.reverse.all.vcf', shell=True)
	subprocess.call('cat '+vcf_head_name+' '+vcf_sorted+'.reverse.all.vcf > '+final_vcf_out+'.reverse.all.vcf', shell=True)

subprocess.call('rm '+vcf_head_name+'* '+vcf_merged+'* '+vcf_sorted+'* '+cwd+'/chr*.temp*.vcf '+cwd+'/chr*.temp*.bed ', shell=True)

##########################################################################