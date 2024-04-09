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

def Define_CombineGVCFs_Output():
	commands=[
		'chr=${name}\n',
		'chr=${chr##*_}\n',
		#'i=${i%%.sh}\n',
		#'for chrom in ${chr[i]}\n',
		#'do\n',
		'ref=${ref_chr}.${chr}.fasta\n',
		'java -jar ${GATK} -T CombineGVCFs \\\n',
		'\t-R ${ref} \\\n',
		'\t-V ${vcfs[@]} \\\n',
		'\t-o ${outpath4step}/${name}.raw.${chr}.gvcf.gz \n\n']
		#'java -jar ${GATK} -T GenotypeGVCFs \\\n',
		#'\t-R ${ref} -V ${outpath4step}/${name}.raw.${chr}.gvcf.gz\\\n',
		#'\t-o ${outpath4step}/${name}.raw.${chr}.vcf.gz\n',
		#'done\n\n'
	return(commands)
	#output=sum([require_setting,log1,commands,log2],[])
	#return(output)

if __name__ == '__main__':
	##setting
	soft_opt=["GATK"]
	param_opt=["chr"]
	ref_opt=["ref","bed","ref_chr"]
	in_opt=['inputtarget','vcfs']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
        step='CombineGVCFs'
	commands=Define_CombineGVCFs_Output()
	### Basic info, software, parameter,reference
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
