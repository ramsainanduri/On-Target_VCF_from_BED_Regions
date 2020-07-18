#!/usr/bin/python
# This code is compatable with python3
# Developed by Ram Sai Nanduri
# Github: https://github.com/ramsainanduri
# Script to extract vcf records that falls in the given bed region

############################### Required Packages ###########################################

import os
import sys
import argparse
import time

############################### Usage Block ###########################################
parser = argparse.ArgumentParser(description="Extract vcf records that falls in the given bed region, The current algorithm will hold all requested vcf records in memory prior to output. The user must ensure that there is adequate memory for this.")
parser.add_argument("-i", "--vcf", metavar="vcf", required=True, help="Input vcf file", type=str)
parser.add_argument("-b", "--bed", metavar="bed", required=True, help="Input Bed file", type=str)
parser.add_argument("-o", "--output", metavar="out_vcf", required=True, help="Output vcf file", type=str)
parser.add_argument("-r", "--offtarget", metavar="offtarget", default="N", help="Y, prints vcf records that doesnot fall in the given bed region (off-target), default N", type=str)
parser.add_argument("-v", "--version", help="Program's version", action='version', version='%(prog)s 1.0')
args = parser.parse_args()

########################## General Variable Defining ######################################
start_time = time.clock()
cwd = os.getcwd()
folder = os.path.basename(cwd)
vcf_ori = cwd+'/'+args.vcf
bed_ori = cwd+'/'+args.bed
offtarget = args.offtarget

if offtarget == "N":
	final_vcf_out = cwd+'/'+args.output+".on-target.vcf"
elif offtarget == "Y":
	final_vcf_out = cwd+'/'+args.output+".off-target.vcf"
else:
	print("Please provide Y or N for the on-target vcf output!")

########################### Deleteing if exsiting and output file ##############################

if os.path.exists(final_vcf_out):
	os.remove(final_vcf_out)
out_file=open(final_vcf_out, "a+")

########################### Reading vcf and Bed file ##############################

#vcf file reading
with open(vcf_ori, "r") as vcf:
	lines = vcf.read().split('\n')
	header = lines
	header = [vh for vh in lines if "#" in vh]
	vcf_records = [vr for vr in lines if "#" not in vr]

while("" in vcf_records):
	vcf_records.remove("")
vcf.close()

#bed file reading
with open(bed_ori, "r") as bed:
	reg = bed.read().split('\n')
	bed_records_tmp = [br for br in reg if "#" not in br]
	bed_records = [txt.replace("\r", "") for txt in bed_records_tmp]

while("" in vcf_records):
	vcf_records.remove("")
bed.close()

############################## Writing output after matching on-target and off-target ##############################
vcf_ontarget = []
vcf_offtarget = []
for i in header:
	out_file.write(i+'\n')

for vcf_line in vcf_records:
	vcf_elements = vcf_line.split('\t')

	for bed_region in bed_records:
		bed_elements = bed_region.split('\t')

		if offtarget == "N":	
			if (vcf_elements[0] == bed_elements[0]) and (int(vcf_elements[1]) > int(bed_elements[1]) and int(vcf_elements[1]) <= int(bed_elements[2])):
				#vcf_ontarget.append(vcf_line)
				out_file.write(vcf_line+"\n")
				break
	
		elif offtarget == "Y":
			if (vcf_elements[0] == bed_elements[0]) and (int(vcf_elements[1]) <= int(bed_elements[1]) or int(vcf_elements[1]) > int(bed_elements[2])):
				#vcf_offtarget.append(vcf_line)
				out_file.write(vcf_line+"\n")
				break		

out_file.close()

#Resetting all the larger variables
vcf_ontarget = vcf_offtarget = vcf_records = bed_records = lines = bed = bed_records_tmp = header = []
print("Program ran for "+ str(round((time.clock() - start_time)/60, 3))+" minutes")
########### Finished #####################3
