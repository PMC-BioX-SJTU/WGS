#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
sys.path.append('/opt/NfsDir/UserDir/sunfl/Codes/Pipelines/BasicWESorWGS')
from PreProcess import *
cf = ConfigParser.ConfigParser()

parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outpath4proj for the project')

if len(sys.argv) <= 3: 
	parser.print_help() 
	sys.exit(1) 
else: 
	args = parser.parse_args()
config = args.config
outpath4proj = args.outpath4proj
sample = args.sample


def Define_Annotation_Output(requirements):
	sample=['name=',requirements[0][0],'\n','outpath4proj=',requirements[0][2],'\n',
		'outpath4step=',requirements[0][3],'\n','inputtarget=',requirements[0][1],'\n\n']
	software=['annovar=',requirements[1][1],'\n','module_path=',requirements[1][0],'\n\n']
	#parameter=['vcftype=',requirements[2][0],'\n\n']
	reference=['database=',requirements[3][0],'\n','operation=',requirements[3][1],'\n\n']
	if(requirements[2][0]=='multiple'):
		para1='-format vcf4 -allsample -withfreq'
		para2='annovar'
	else:
		para1='--format vcf4 --allsample'
                para2='annovar.${name}.avinput'
	log1=['echo -e The Step Annotation started at `date` "\\n" >>${outpath4step}/run.log;\n']
	log2=['echo -e The Step Annotation ended at `date` "\\n" >>${outpath4step}/run.log;\n',
		'mv ${outpath4step}/run.log ${outpath4step}/Annotation.log;\n\n']
	### for singe sample analysis
	anno_prep=['perl ${annovar}/convert2annovar.pl --includeinfo '+para1+' --outfile ${outpath4step}/annovar ${inputtarget}\n\n']
	annotation=['perl ${annovar}/table_annovar.pl --buildver hg19 \\\n',
		'\t--thread 3 --remove  --otherinfo \\\n',
		'\t--protocol ${database} \\\n',
		'\t-operation ${operation}  \\\n',
		'\t-nastring . ${outpath4step}/'+para2+' ${annovar}/humandb \\\n',
		'\t--outfile ${outpath4step}/${name}.annovar > ${outpath4step}/${name}.annovar.log\n\n']
	stat=['perl ${module_path}/annovar_statistic.pl  -a ${outpath4step}/${name}.annovar.hg19_multianno.txt -s ${name}\n']
	other=['perl ${module_path}/add_title.pl  ${inputtarget}  \\\n',
		'\t ${outpath4step}/${name}.annovar.hg19_multianno.txt   \\\n',
		'\t>  ${outpath4step}/${name}.annovar.txt\n\n',
		'perl  ${module_path}/filter_AF.pl  \\\n',
		'\t${outpath4step}/${name}.annovar.txt \\\n',
		'\t0.01  ${type} >${name}.annovar_filtered.txt\n\n']
	output1=sum([sample,software,reference,log1,anno_prep,annotation,stat,other,log2],[])
	### for family inheritance model annotation
	fam_anno=['for inm in `ls ${outpath4proj}/InheritanceModels/*inm`\n',
		'do\n',
		'name=${inm%%.inm}\n',
		'name=${name##*/}\n',
		'for type in `cat ${inm}`\n',
		'do\n',
		'inputtarget=${outpath4proj}/InheritanceModels/${name}.primary.${type}.vcf\n',
		'perl ${annovar}/convert2annovar.pl --includeinfo '+para1+' --outfile ${outpath4step}/annovar ${inputtarget}\n\n',
		'perl ${annovar}/table_annovar.pl --buildver hg19  \\\n',
		'\t--thread 3 --remove  --otherinfo \\\n',
		'\t--protocol ${database} \\\n',
		'\t-operation ${operation} \\\n',
		'\t-nastring . ${outpath4step}/'+para2+' ${annovar}/humandb \\\n',
		'\t--outfile ${outpath4step}/${name}.${type}.annovar \\\n',
		'\t> ${outpath4step}/${name}.annovar.log\n\n',
		'perl ${module_path}/annovar_statistic.pl  \\\n',
		'\t-a ${outpath4step}/${name}.${type}.annovar.hg19_multianno.txt -s ${name}.${type}\n',
		'perl ${module_path}/add_title.pl  ${inputtarget}  \\\n',
		'\t ${outpath4step}/${name}.${type}.annovar.hg19_multianno.txt \\\n',
		'\t>  ${outpath4step}/${name}.${type}.annovar.txt\n\n',
		'perl  ${module_path}/filter_AF.pl \\\n',
		'\t${outpath4step}/${name}.${type}.annovar.txt \\\n',
		'\t0.01  ${type} >${name}.${type}.annovar_filtered.txt\n',
		'done\n','done\n']
	output2=sum([sample,software,reference,log1,anno_prep,annotation,stat,other,fam_anno,log2],[])
	output3=sum([sample,software,reference,log1,fam_anno,log2],[])
	if(requirements[2][1]=='NO' and requirements[2][2]=='NO'):
		return (output1)
	elif(requirements[2][1]=='YES' and requirements[2][2]=='YES'):
		return (output2)
	elif(requirements[2][1]=='YES' and requirements[2][2]=='NO'):
		return (output3)

	
if __name__ == '__main__':
	#requirements preparation
	soft_opt=["module_path","annovar"]
	param_opt=["sample","InheritanceModels","DeNovo"]
	ref_opt=["database","operation"]
	opt=[soft_opt,param_opt,ref_opt]
	sec=['software','Annotation.parameter','reference']
	### SampleAndOutpat, software, parameter,reference
	requirements=Input_Requiremens(config,sec,opt,sample,"Annotation",outpath4proj)
	### write makefile
	makefilename=requirements[0][3]+'/Annotation.sh'
	f=open(makefilename,"w")
	output=Define_Annotation_Output(requirements)
	f.writelines(output)
	f.close
