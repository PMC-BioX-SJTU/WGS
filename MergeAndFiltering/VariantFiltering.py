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

	
def Define_VariantFiltering_Output():
	commands=['if [ ${method} != VQSR ] ;then \n',
		'\t## harder filering of variant sites, Extract the SNPs from the call set\n',
		'echo -e The Step VariantFiltering by hard filtering started at `date` "\\n" >>${outpath4step}/run.log\n',
		'java -jar ${GATK}  -T SelectVariants \\\n',
		'\t-R ${ref}  -V ${inputtarget}  \\\n',
		'\t-selectType SNP -o ${outpath4step}/${name}.raw_snp.vcf;\n',
		'\t## Extract the Indels from the call set\n',
		'java -jar ${GATK}  -T SelectVariants \\\n',
		'\t-R ${ref}  -V ${inputtarget} \\\n',
		'\t-selectType INDEL -o ${outpath4step}/${name}.raw_indel.vcf\n\n',
		'\t## Determine parameters for filtering SNPs \n',
		'\t## QD 2.0; FS 60.0; MQ 40.0; MQRankSum -12.5; ReadPosRankSum -8.0; DP 10\n',
		'\t## Apply the filter to the SNP call set\n',
		'java -jar  ${GATK} -T VariantFiltration \\\n',
		'\t-R ${ref} -V ${outpath4step}/${name}.raw_snp.vcf  \\\n',
		'\t--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 ||DP < ${DP} " \\\n',
		'\t--filterName "my_filter" -o  ${outpath4step}/${name}.filtered_snps.vcf\n\n',
		'\t### Determine parameters for filtering Indels.\n',
		'\t### QD 2.0; FS 200.0; ReadPosRankSum 20.0; DP 10\n',
		'\t### Apply the filter to the Indel call set\n',
		'java -jar ${GATK} -T VariantFiltration \\\n',
		'\t-R ${ref} -V  ${outpath4step}/${name}.raw_indel.vcf  \\\n',
		'\t--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || DP < ${DP}"  \\\n',
		'\t--filterName "my_filter"  -o ${outpath4step}/${name}.filtered_indels.vcf\n\n',
		'grep -v \'^#\' ${outpath4step}/${name}.filtered_indels.vcf >> ${outpath4step}/${name}.filtered_snps.vcf;\n',
		'grep -v \'my_filter\' ${outpath4step}/${name}.filtered_snps.vcf >${outpath4step}/${name}.final1.vcf\n',
		'perl ${module_path}/DepthFiltering.pl all.final1.vcf ${DP} >all.final2.vcf',
		'perl ${module_path}/sort_chr_vcf.pl ${outpath4step}/${name}.final1.vcf >${outpath4step}/${name}.final.vcf\n',
		'perl ${module_path}/sort_chr_vcf.pl ${outpath4step}/${name}.final2.vcf >${outpath4step}/${name}.filtered.vcf\n',
		'${igvtools} index ${outpath4step}/${name}.filtered.vcf\n',
		'${igvtools} index ${outpath4step}/${name}.final.vcf\n',
		'echo -e The Step VariantCalling by hard filter ended at `date` "\\n" >>${outpath4step}/run.log\n',
		'mv ${outpath4step}/run.log ${outpath4step}/VariantFiltering.log\n',
		'else\n',
		'echo -e The Step VariantFiltering started at `date` "\\n" >>${outpath4step}/run.log\n',
		'java -jar ${picard} MakeSitesOnlyVcf \\\n',
		'\tINPUT=${inputtarget} \\\n', 
		'\tOUTPUT=${outpath4step}/output.sites_only.unfiltered.vcf.gz\n\n',
		#perform SNP VQSR on sites only file
		'java -jar ${GATK} -T VariantRecalibrator \\\n',
		'\t--disable_auto_index_creation_and_locking_when_reading_rods \\\n',
		'\t-R ${ref} -input ${outpath4step}/output.sites_only.unfiltered.vcf.gz \\\n',
		'\t--num_threads 2 -recalFile ${outpath4step}/output.snps.recal -tranchesFile ${outpath4step}/output.snps.tranches -allPoly \\\n',
		'\t-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 \\\n',
		'\t-tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \\\n',
		'\t-an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an InbreedingCoeff \\\n',
		'\t-resource:hapmap,known=false,training=true,truth=true,prior=15 ${hapmap} \\\n',
		'\t-resource:omni,known=false,training=true,truth=true,prior=12 ${omni} \\\n',
		'\t-resource:1000G,known=false,training=true,truth=false,prior=10 ${known3} \\\n',
		'\t--resource:dbsnp137,known=false,training=false,truth=false,prior=7 ${dbsnp} \\\n',
		'\t--resource:dbsnp129,known=true,training=false,truth=false,prior=3 ${dbsnp1} \\\n',
		'\t--maxGaussians 6 -mode SNP -rscriptFile ${outpath4step}/output.snps.recalibration_plots.rscript >> ${outpath4step}/vqsr.log\n\n',
		#perform INDEL VQSR on sites only file
		'java -jar ${GATK} -T VariantRecalibrator \\\n',
		'\t--disable_auto_index_creation_and_locking_when_reading_rods \\\n',
		'\t-R ${ref} -input ${outpath4step}/output.sites_only.unfiltered.vcf.gz \\\n',
		'\t--num_threads 1 -recalFile ${outpath4step}/output.indels.recal -tranchesFile ${outpath4step}/output.indels.tranches -allPoly \\\n',
		'\t-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 \\\n',
		'\t-tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 \\\n',
		'\t-tranche 91.0 -tranche 90.0 \\\n',
		'\t-an FS -an ReadPosRankSum -an InbreedingCoeff -an MQRankSum -an QD \\\n',
		'\t-resource:mills,known=false,training=true,truth=true,prior=12 ${known2} \\\n', 
		'\t-resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiom}  \\\n',
		'\t-resource:dbsnp137,known=true,training=false,truth=false,prior=2 ${dbsnp} \\\n',
		'\t--maxGaussians 6 -mode INDEL -rscriptFile ${outpath4step}/output.indels.recalibration_plots.rscript >> ${outpath4step}/vqsr.log\n\n',
		#Apply SNP VQSR on genotypes file (SNP sensitivity 99.6)
		'java -jar ${GATK} -T ApplyRecalibration \\\n',
		'\t--disable_auto_index_creation_and_locking_when_reading_rods \\\n',
		'\t-R ${ref} -o ${outpath4step}/output.snps.filtered.vcf.gz  \\\n',
		'\t-input ${inputtarget} -recalFile ${outpath4step}/output.snps.recal \\\n',
		'\t-tranchesFile ${outpath4step}/output.snps.tranches \\\n',
		'\t-ts_filter_level 99.6 -mode SNP >> ${outpath4step}/vqsr.log \n\n',
		#Apply Indel VQSR on genotypes file (INDEL sensitivity 95.0) # snp file as input file
		'java -jar ${GATK} -T ApplyRecalibration \\\n',
		'\t--disable_auto_index_creation_and_locking_when_reading_rods \\\n',
		'\t-R ${ref} -o ${outpath4step}/all.final.vcf  \\\n',
		'\t-input ${outpath4step}/output.snps.filtered.vcf.gz \\\n',
		'\t-recalFile ${outpath4step}/output.indels.recal \\\n',
		'\t-tranchesFile ${outpath4step}/output.indels.tranches \\\n',
		'\t-ts_filter_level 95.0 -mode INDEL >> ${outpath4step}/vqsr.log \n\n',
		'\tjava -jar  ${GATK} -T SelectVariants \\\n',
		'\t-R ${ref} \\\n',
		'\t-V ${outpath4step}/all.final.vcf \\\n',
		'\t-o ${outpath4step}/all.filtered.vcf \\\n',
		'\t--setFilteredGtToNocall -ef -select "QD > 10.0 "\n\n',
		'echo -e The Step VariantFiltering by VQSR ended at `date` "\\n" >>${outpath4step}/run.log\n',
		'fi\n\n']
	return(commands)


def Makeclean_VariantFilgering(path,name):
        clean=['clean:\n',
                '\trm'+path+name+'.raw*.vcf;\n',
                '\trm'+path+name+'.filtered*.vcf;\n']
        retrun (clean)

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["GATK","picard","igvtools"]
	param_opt=["method","DP"]
	ref_opt=["ref","dbsnp","dbsnp1","known3","omni","axiom","hapmap","known2"]
	in_opt=['inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	### Basic info, software, parameter,reference
	step='VariantFiltering'
	commands=Define_VariantFiltering_Output()
	Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
