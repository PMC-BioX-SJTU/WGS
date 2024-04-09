#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/home/zhujh/Codes/Pipelines/BasicWESorWGSPlus')
sys.path.append('/home/zhujh/Codes/Pipelines/FamilyBasic')
from BasicWESorWGSLib import *
from FamilyBasicLib import *

cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for obtain Variants that follow Inheritance Models!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='family ped info')
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
  

	
def Define_ADARXlinkFinderAndFiltering_Output():
	commands=[
                '#perl  ${module_path}/AD_AR_X_Y_selection.pl \\\n',
                '#\t${outpath4step}/../${name}.ped \\\n',
                '#\t${inputtarget} \\\n',
                '#\t${outpath4step}/${name}  1 \n',
		'for inm in `ls ${outpath4step}/*inm`\n',
		'do\n',
		'name=${inm%%.inm}\n',
		'name=${name##*/}\n',
		'for type in `cat ${inm}`\n',
		'do\n',
		#'inputtarget=${outpath4step}/${name}.primary.${type}.vcf\n',
		#'if [ ${sample} == multiple ] ; then\n',
                #'\tpara1="-format vcf4 -allsample -withfreq"\n',
                #'\tpara2=${name}.${type}.annovar\n',
                #'else\n',
                #'\tpara1="--format vcf4 --allsample"\n',
                #'\tpara2=annovar.${name}.avinput\n',
                #'fi\n\n',
		#'perl ${annovar}/convert2annovar.pl --includeinfo ${para1} \\\n',
		#'\t--outfile ${outpath4step}/${name}.${type}.annovar \\\n',
		#'\t${inputtarget}\n\n',
		#'perl ${annovar}/table_annovar.pl --buildver hg19 \\\n',
           	#'\t--thread 3 --remove  --otherinfo \\\n',
                #'\t--protocol ${database} \\\n',
                #'\t-operation ${operation}  \\\n',
                #'\t-nastring . ${outpath4step}/${para2} ${annovar}/humandb \\\n',
                #'\t--outfile ${outpath4step}/${name}.${type}.annovar > ${outpath4step}/${name}.${type}.annovar.log\n\n',
        	#'perl ${module_path}/add_title.pl  ${inputtarget}  \\\n',
                #'\t ${outpath4step}/${name}.${type}.annovar.hg19_multianno.txt   \\\n',
                #'\t>  ${outpath4step}/${name}.${type}.annovar.txt\n\n',
		
		'perl  ${module_path}/ADARXlinkFiltering.pl \\\n',
		'\t${outpath4step}/${name}.${type}.annovar.txt \\\n',
		'\t0.01 ${name} ${outpath4step}/${name}.${type}.annovar_filtered \n\n',
		
                #'perl  ${module_path}/filter_AF.pl  \\\n',
                #'\t${outpath4step}/${name}.${type}.annovar.txt \\\n',
                #'\t0.01 >${outpath4step}/${name}.${type}.annovar_filtered.txt\n\n',
		'done\n',
		'done\n\n'
		]
	return commands

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["annovar","GATK"]
	param_opt=["sample"]
	ref_opt=["dbsnp","ref","database","operation","ped"]
	in_opt=['inputtarget']
        out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='ADARXlinkFinderAndFiltering'
	commands=Define_ADARXlinkFinderAndFiltering_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

