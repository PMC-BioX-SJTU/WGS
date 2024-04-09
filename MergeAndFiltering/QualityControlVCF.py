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

	
def Define_QualityControlVCF_Output():
	dbsnp=['###the ped file in the part must be in concordance with the sample names in the header of the vcf file\n',
		'java -jar ${GATK} -T VariantAnnotator \\\n', 
		'\t-R ${ref} -V ${inputtarget}  \\\n',
		'\t-o ${outpath4step}/${name}.annotated.vcf  \\\n',
		'\t--dbsnp ${dbsnp}\n']
	#filter_snp=['awk -F "\t" \'{sum=0;for(i=10;i<=NF;i++) \\\n',
	#	'\t{split($i,a_i,":");if(a_i[4]>=20)sum=sum+1}; \\\n',
	#	'\tif(($3!="." && sum>=NF-9)||($1~/^\#/))  print $0}\' \\\n',
	#	'\t ${outpath4step}/${name}.annotated.vcf\\\n',
	#	'\t > ${outpath4step}/${name}.dbsnp.vcf\n']
	binary_ped=['${vcftools} --vcf ${outpath4step}/${name}.annotated.vcf \\\n', 
		'\t--plink  --out ${outpath4step}/all.plink\n\n',
		'${plink} --file ${outpath4step}/all.plink --missing --out ${outpath4step}/all.plink \n\n',
		'awk \'NR < 2 { next } $5 > 0\' ${outpath4step}/all.plink.lmiss \\\n',
		'\t > ${outpath4step}/all.plink.lmiss.exclude \n',
		'${plink} --file ${outpath4step}/all.plink \\\n',
		'\t--exclude ${outpath4step}/all.plink.lmiss.exclude \\\n',
		'\t--make-bed --out all.plink_QC1 \n\n',
		'cp ${ped} ${outpath4step}/${name}.plink_QC1.fam \n\n']
		#'java -jar ${GATK} -T VariantsToBinaryPed \\\n',
		#'\t-R ${ref} -V ${outpath4step}/${name}.dbsnp.vcf  \\\n',
		#'\t-m  ${outpath4step}/${name}.data.fam  -mgq 0 \\\n',
		#'\t-bed ${outpath4step}/${name}.plink.bed \\\n',
		#'\t-bim ${outpath4step}/${name}.plink.bim \\\n',
		#'\t-fam ${outpath4step}/${name}.plink.fam\n'


	ibd=['${plink} --bfile \\\n',
		'\t${outpath4step}/${name}.plink_QC1 \\\n',
		'\t--maf 0.01 --genome  \\\n',
		'\t--out ${outpath4step}/${name}.ibd\n\n']

	ibd_sex=['${plink} --bfile \\\n',
		'\t${outpath4step}/${name}.plink_QC1 \\\n',
		'\t--check-sex \\\n',
		'\t--out ${outpath4step}/${name}.sex\n\n']

	merged_vcf_qc=['java -jar ${GATK} -T VariantEval  \\\n',
		'\t-R ${ref} \\\n',
		'\t--eval ${inputtarget} \\\n', 
		'\t-D ${dbsnp} -noST -ST Sample \\\n', 
		'\t-noEV -EV CompOverlap -EV IndelSummary \\\n',
		'\t-EV TiTvVariantEvaluator -EV CountVariants \\\n',
		'\t-EV MultiallelicSummary -o ${outpath4step}/vcf.rawqc.txt \n']

	narr_qc=['perl ${module_path}/QCvcf_table.pl \\\n',
		'\t${outpath4step}/vcf.rawqc.txt \\\n',
		'\t>${outpath4step}/vcf.qc.xls\n']

	output=sum([dbsnp,binary_ped,ibd,ibd_sex,merged_vcf_qc,narr_qc],[])
	return output

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["GATK","plink","vcftools"]
	param_opt=["method"]
	ref_opt=["ref","dbsnp","ped"]
	in_opt=['inputtarget']
	out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='QualityControlVCF'
        commands=Define_QualityControlVCF_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)
