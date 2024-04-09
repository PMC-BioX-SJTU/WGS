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
module= args.module
def Define_MergeQCs_Output():
	commands=['paste ${Q20Q30[@]} > ${outpath4step}/Merge.QC.raw.Q20Q30\n\n',
		'paste ${Alig[@]} > ${outpath4step}/Merge.QC.raw.Alig\n\n',
		'paste ${CrossCont[@]} > ${outpath4step}/Merge.QC.raw.CrossCont\n\n',
		'paste ${Hs[@]} > ${outpath4step}/Merge.QC.raw.Hs\n\n',
		'perl ${module_path}/QCWriteExcel.pl  \\\n',
		'\t${outpath4step}/Merge.QC.raw.Q20Q30 \\\n',
		'\t${outpath4step}/Merge.QC.raw.Alig  \\\n',
		'\t${outpath4step}/Merge.QC.raw.CrossCont \\\n',
		'\t${outpath4step}/Merge.QC.raw.Hs ${outpath4step}/QC.merge.xls\n\n' ]
	return(commands)	

if __name__ == '__main__':
	##requirements preparation
	soft_opt=["GATK"]
	param_opt=["method"]
	ref_opt=["ref","bed","ref_chr"]
	in_opt=['inputtarget','Q20Q30','Alig','CrossCont','Hs']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	### Basic info, software, parameter,reference
	step='MergeQCs'
	commands=Define_MergeQCs_Output()
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
