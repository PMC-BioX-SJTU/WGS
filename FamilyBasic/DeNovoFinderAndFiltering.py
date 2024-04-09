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
  

	
def Define_DeNovoFinderAndFiltering_Output():
	commands=[
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
		
		'if [ ${sample} == multiple ] ; then\n',
                '\tpara1="-format vcf4 -allsample -withfreq"\n',
                '\tpara2=${name}.annovar\n',
                'else\n',
                '\tpara1="--format vcf4 --allsample"\n',
                '\tpara2=annovar.${name}.avinput\n',
                'fi\n\n',

		'perl ${annovar}/convert2annovar.pl --includeinfo ${para1} \\\n',
		'\t--outfile ${outpath4step}/${name}.annovar \\\n',
		'\t${outpath4step}/${name}.denovo.vcf\n\n',

		'perl ${annovar}/table_annovar.pl --buildver hg19 \\\n',
           	'\t--thread 3 --remove  --otherinfo \\\n',
                '\t--protocol ${database} \\\n',
                '\t-operation ${operation}  \\\n',
                '\t-nastring . ${outpath4step}/${para2} ${annovar}/humandb \\\n',
                '\t--outfile ${outpath4step}/${name}.annovar > ${outpath4step}/${name}.annovar.log\n\n',

        	'perl ${module_path}/add_title.pl  ${outpath4step}/${name}.denovo.vcf  \\\n',
                '\t ${outpath4step}/${name}.annovar.hg19_multianno.txt   \\\n',
                '\t>  ${outpath4step}/${name}.annovar.txt\n\n',

                'perl  ${module_path}/filter_AF.pl  \\\n',
                '\t${outpath4step}/${name}.annovar.txt \\\n',
                '\t0.01 >${outpath4step}/${name}.filtered.per1.txt\n\n'

		'perl ${module_path}/count_deonvo.pl \\\n',
		'\t${outpath4step}/${name}.annovar.txt  1.07 without \\\n',
		'\t${outpath4step}/${name}.filtered.txt  \\\n', 
		'\t>${outpath4step}/${name}.count.stat.txt\n\n',

		'perl  ${module_path}/MPC/missense_anno.pl  \\\n',
		'\t${module_path}/MPC/miss_baddness.txt   \\\n',
		'\t${module_path}/MPC/gamma.txt   \\\n',
		'\t${module_path}/MPC/segment.txt \\\n',
		'\t${outpath4step}/${all}.filtered.txt > ${outpath4step}/MPC_cand.txt\n\n'
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
	step='DeNovoFinderAndFiltering'
	commands=Define_DeNovoFinderAndFiltering_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

