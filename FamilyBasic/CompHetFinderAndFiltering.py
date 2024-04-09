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
  

	
def Define_CompHetFinderAndFiltering_Output():
	commands=[
                'if [ ${DEL} == NO ] ; then\n',
		'\tperl  ${module_path}/compound_het.pl \\\n',
                '\t${outpath4step}/../${name}.ped \\\n',
                '\t${inputtarget} \\\n',
                '\t${refGene} \\\n',
                '\t${outpath4step}/${name}.CompHet.vcf \n\n',
		'else\n',
		'\tperl ${module_path}/CompHetCNVSNVINDEL.pl \\\n',
        	'\t${outpath4step}/../${name}.ped \\\\n',
        	'\t${inputtarget} \\\n',
        	'\t${refGene} \\\n',
		'\t${outpath4step}/${name}.CompHet.vcf \\\n',
        	'\t${DEL}\n\n',
		'fi\n\n',
		'if [ ${sample} == multiple ] ; then\n',
                '\tpara1="-format vcf4 -allsample -withfreq"\n',
                '\tpara2=${name}.annovar\n',
                'else\n',
                '\tpara1="--format vcf4 --allsample"\n',
                '\tpara2=annovar.${name}.avinput\n',
                'fi\n\n',
		'perl ${annovar}/convert2annovar.pl --includeinfo ${para1} \\\n',
		'\t--outfile ${outpath4step}/${name}.annovar \\\n',
		'\t${outpath4step}/${name}.CompHet.vcf\n\n',
		'perl ${annovar}/table_annovar.pl --buildver hg19 \\\n',
           	'\t--thread 3 --remove  --otherinfo \\\n',
                '\t--protocol ${database} \\\n',
                '\t-operation ${operation}  \\\n',
                '\t-nastring . ${outpath4step}/${para2} ${annovar}/humandb \\\n',
                '\t--outfile ${outpath4step}/${name}.annovar > ${outpath4step}/${name}.annovar.log\n\n',
        	'perl ${module_path}/add_title.pl  ${outpath4step}/${name}.CompHet.vcf  \\\n',
                '\t ${outpath4step}/${name}.annovar.hg19_multianno.txt   \\\n',
                '\t>  ${outpath4step}/${name}.annovar.txt\n\n',
                'perl  ${module_path}/CompHetFilteringNoDEL.pl  \\\n',
                '\t${outpath4step}/${name}.annovar.txt \\\n',
                '\t0.01  ${name} ${outpath4step}/${name}.annovar_filtered\n\n'
		]
	return commands

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["annovar"]
	param_opt=["sample","DEL"]
	ref_opt=["database","refGene","operation"]
	in_opt=['inputtarget']
        out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='CompHetFinderAndFiltering'
	commands=Define_CompHetFinderAndFiltering_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

