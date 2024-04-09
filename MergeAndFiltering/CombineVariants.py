#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
import re
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus')
sys.path.append('/home/zhujh/Codes/Pipelines/FamilyBasic')
from BasicWESorWGSLib import *
from FamilyBasicLib import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outfile for make')
parser.add_argument('-m', required=True, dest='module', action='store',help='module name')

if len(sys.argv) <= 3: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outpath4proj = args.outpath4proj
sample = args.sample
module=args.module

def Define_MergeStep_Output(require_setting):
	chr_group=['chr1 chr2',
		'chr3 chr4 chr5',
		'chr6 chr7 chr8',
		'chr9 chr10 chr11 chr12',
		'chr13 chr14 chr15 chr16 chr17',
		'chr18 chr19 chr20 chr21 chr22 chrX chrY chrM']


def Define_CombineVariants_Output():
	commands=[
		'java -jar ${GATK} -T CombineVariants \\\n',
		'\t-R ${ref}  \\\n',
		'\t--variant  ${vcf[@]} \\\n',   
		'\t-o ${outpath4step}/all.raw.vcf.gz -genotypeMergeOptions UNSORTED  --setKey null\n\n',
		]
	return(commands)

if __name__ == '__main__':
	##setting
	soft_opt=["GATK"]
	param_opt=["mergecalling","datatype"]
	ref_opt=["ref","bed","ref_chr"]
	in_opt=['inputtarget','vcf']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='CombineVariants'
	commands=Define_CombineVariants_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

