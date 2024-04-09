#/usr/bin/python

import ConfigParser
import string
import os
import argparse
import sys
from BasicWESorWGSLib import *
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

def Define_MappingAndCleanUp_Output():
	BQSR=[	'java -jar ${GATK} -T BaseRecalibrator \\\n',
		'\t-R ${ref} -I ${outpath4step}/${name}.align.reorder.sorted.markdup.bam \\\n',
		'\t-knownSites ${dbsnp} -nct ${thread} -knownSites ${known1} -knownSites ${known2} \\\n',
		'\t-o ${outpath4step}/${name}.recal.table\n\n',

		'java -jar ${GATK} -T PrintReads \\\n',
		'\t-R ${ref} -I ${outpath4step}/${name}.align.reorder.sorted.markdup.bam \\\n', 
		'\t-BQSR ${outpath4step}/${name}.recal.table -nct  ${thread}\\\n',
		'\t-o ${outpath4step}/${name}.recal.bam;\n']

	IndelRealign=['java -jar ${GATK} -T IndelRealigner \\\n',
		'\t-R ${ref} -I ${outpath4step}/${name}.align.reorder.sorted.markdup.bam  \\\n',
		'\t-known ${known1} -known ${known2} \\\n',
		'\t-targetIntervals ${outpath4step}/${name}.realigner.intervals \\\n',
		'\t-o ${outpath4step}/${name}.realigner.bam\n\n']

	RealignerTargetCreator=['echo -e The Step IndelRealignAndBQSR started at `date` "\\n" >>${outpath4step}/run.log;\n',
		'java -jar ${GATK} -T RealignerTargetCreator \\\n',
		'\t-R ${ref} -I ${outpath4step}/${name}.align.reorder.sorted.markdup.bam \\\n',
		'\t-known ${known1} -known ${known2} -nt ${thread}\\\n',
		'\t-o ${outpath4step}/${name}.realigner.intervals\n\n']
	make_dup=['java -jar ${picard} MarkDuplicates \\\n',
		'\tI= ${outpath4step}/${name}.align.reorder.sorted.bam \\\n',
		'\tO= ${outpath4step}/${name}.align.reorder.sorted.markdup.bam \\\n',
		'\tASSUME_SORTED=true METRICS_FILE= ${outpath4step}/${name}.align.metrics \\\n',
		'\tVALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true;\n',
		'\techo -e The Step MappingAndMarkDup ended at `date` "\\n" >>${outpath4step}/run.log;\n\n']

	samtools_index=['${samtools} index ${outpath4step}/${name}.align.reorder.sorted.bam\n\n']

	sort=[	'${samtools} sort  -@  ${thread} -m 2G\\\n',
		'\t${outpath4step}/${name}.align.reorder.bam \\\n',
		'\t-o ${outpath4step}/${name}.align.reorder.sorted.bam \n\n']

	reorder=['java -jar ${picard} ReorderSam \\\n',
		'\tI= ${outpath4step}/${name}.align.bam \\\n',
		'\tO= ${outpath4step}/${name}.align.reorder.bam \\\n',
		'\tR=${ref} VALIDATION_STRINGENCY=LENIENT\n\n']
	make_clean=[ 'if [ ${clean} == "YES" ] ;then\n',
		'\trm ${outpath4step}/${name}*.align.*;\n',
                '\trm ${outpath4step}/${name}.realigner.*;\n',
                '\trm ${outpath4step}/${name}*recal.table\n',
		'fi\n\n']
		
	

	mapping=['echo -e The Step MappingAndMarkDup started at `date` "\\n" >>${outpath4step}/run.log;\n',
		'bam_list=""\n',
		'for(( i=0;i<${#fq1[@]};i++))\n',
		'do\n',
		'${bwa} mem -t ${thread}  -M \\\n',
		'\t-R "@RG\\tID:QSY${name}\\tLB:QSY${name}\\tSM:${name}\\tPL:ILLUMINA" \\\n',
		'\t${ref}  ${fq1[i]} ${fq2[i]}  \\\n',
		'\t|gzip -3 > ${outpath4step}/${name}fq${i}.align.sam\n\n',

		'${samtools} view  -@ ${thread} -bS \\\n',
		'\t${outpath4step}/${name}fq${i}.align.sam \\\n',
		'\t-o ${outpath4step}/${name}fq${i}.align.bam\n\n',
		'\tbam_list="${bam_list} ${outpath4step}/${name}fq${i}.align.bam"\n',
		'done\n',	
		'if [ ${#fq1[@]} -gt  1 ] ;then\n',
		'\t${samtools} merge -f  ${outpath4step}/${name}.align.bam ${bam_list}\n\n',
		'else\n',
		'\tmv ${outpath4step}/${name}fq0.align.bam ${outpath4step}/${name}.align.bam\n',
		'fi\n\n']
	#output=sum([mapping,reorder,sort,samtools_index,make_dup,RealignerTargetCreator,IndelRealign, BQSR,make_clean],[])
	output=sum([mapping,reorder,sort,samtools_index,make_dup, BQSR,make_clean],[])
	return output

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["bwa","samtools",'picard','GATK']
	param_opt=["bwa.mem","thread","clean"]
	ref_opt=["ref","known1","known2","dbsnp"]
	in_opt=['fq1','fq2','inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='MappingAndCleanUp'
	commands=Define_MappingAndCleanUp_Output()
        
	### write shell scripts
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
	
	
