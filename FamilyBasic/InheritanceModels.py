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
  

	
def Define_InheritanceModels_Output():
	commands=['function ADARXY(){\n',   
		'`ls ${outpath4step}/*ped | split -l 10 - ${outpath4step}/group`\n',
		'for group in `ls ${outpath4step}/group*` \n',
		'do\n',
		'\tfor fam in `cat ${group}`\n',
		'\tdo\n',
		'\tsample=${fam%%.ped}\n',
		'\tperl  ${module_path}/AD_AR_X_Y_selection.pl \\\n',
		'\t${fam} \\\n', 
		'\t${inputtarget} \\\n',
		'\t${sample}  1 &\n',
		'\tdone\n',
		'wait\n',
		'done\n',
		'}\n\n',
		'function CompHet(){\n',
		'`ls ${outpath4step}/*ped | split -l 10 - ${outpath4step}/group`\n',
		'for group in `ls ${outpath4step}/group*` \n',
		'do\n',
		'for fam in `cat ${group} `\n',
		'\tdo\n',
                '\tsample=${fam%%.ped}\n',
		'\tperl  ${module_path}/compound_het.pl \\\n',
		'\t${fam} \\\n',
		'\t${inputtarget} \\\n', 
		'\t${refGene} \\\n',
		'\t${sample}.CompHet.vcf &\n',
		'done\n',
		'wait\n',
		'done\n',
		'}\n\n',
		'function DeNovo(){\n',
		'java -jar ${GATK} -T PhaseByTransmission \\\n', 
		'\t-R ${ref} \\\n',
		'\t-V ${inputtarget}\\\n', 
		'\t-ped ${ped} \\\n',
		'\t-o ${outpath4step}/all.PH.denovo.vcf \n',

		'java -jar ${GATK} -T VariantAnnotator \\\n',
		'\t-R ${ref} \\\n',
		'\t-V ${outpath4step}/all.PH.denovo.vcf\\\n', 
		'\t-ped ${ped} -A PossibleDeNovo \\\n',
		'\t--dbsnp ${dbsnp} \\\n',
		'\t-o ${outpath4step}/all.denovo.raw.vcf\n',
		'egrep \'^#|hiConfDeNovo\' ${outpath4step}/all.denovo.raw.vcf \\\n',
		'\t   > ${outpath4step}/all.denovo_cand.vcf\n',
		'perl ${module_path}/separate_hiConfDeNovo.pl \\\n',
		'\t${outpath4step}/all.denovo_cand.vcf \\\n',
		'\t> ${outpath4step}/all.denovo.vcf\n',
		'}\n\n',
		'perl ${module_path}/split_ped.pl ${ped}  ${outpath4step}\n',	
		'if [ ${DeNovo} == YES ] ;then \n',
		'\tDeNovo\n',
		'fi\n',
		'if [ ${CompHet} == YES ] ;then \n',
		'\tCompHet\n',
		'fi\n',
		'if [ ${ALL} == YES ] ;then \n',
		'\tADARXY\n',
		'fi\n']
	return commands

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["GATK","plink"]
	param_opt=["DeNovo","CompHet","ALL"]
	ref_opt=["ref","dbsnp","refGene","ped"]
	in_opt=['inputtarget']
        out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='InheritanceModels'
	commands=Define_InheritanceModels_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

