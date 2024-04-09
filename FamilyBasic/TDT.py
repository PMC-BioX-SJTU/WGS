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

parser = argparse.ArgumentParser(description='the Step for the Transmission Disequilibrium Test!')
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
  

	
def Define_TDT_Output():
	commands=[
		'###vcf2map\n',
		'java -jar ${GATK} -T SelectVariants \\\n',
		'\t-R ${ref} \\\n',
		'\t-V ${inputtarget} \\\n',
		'\t-o ${outpath4step}/${name}.SNP.final.vcf \\\n',
		'\t-selectType SNP --setFilteredGtToNocall \\\n',
		'\t-ef -select "QD > 10.0 " \n',
		'java -jar ${GATK} -T VariantAnnotator \\\n',
		'\t-R ${ref} -V ${inputtarget}  \\\n',
		'\t-o ${outpath4step}/${name}.annotated.vcf  \\\n',
		'\t--dbsnp ${dbsnp} \n\n',
		'${vcftools} --vcf ${outpath4step}/${name}.annotated.vcf \\\n',
		'\t--plink  \\\n',
		'\t--out ${outpath4step}/${name}.plink \n',

		'##QC1: Remove SNPs with missing genotype rate >5% \n',
		'${plink} --file ${outpath4step}/${name}.plink --missing --out ${outpath4step}/${name}.plink \n',
		'awk \'NR < 2 { next } $5 > 0.05\' ${outpath4step}/${name}.plink.lmiss \\\n',
		'\t> ${outpath4step}/${name}.plink.lmiss.exclude \n',
		'${plink} --file ${outpath4step}/${name}.plink \\\n',
		'\t--exclude ${outpath4step}/${name}.plink.lmiss.exclude \\\n',
		'\t--make-bed --out ${outpath4step}/${name}.plink_QC1 \n',
		'cp ${ped} ${outpath4step}/${name}.plink_QC1.fam\n',
		
		'##QC2: Remove Individuals with missing genotype rate >5% and heterozygosity rate > 3*sigma \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC1 \\\n',
		'\t--missing --out ${outpath4step}/${name}.plink_QC1 \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC1 \\\n',
		'\t--het --out ${outpath4step}/${name}.plink_QC1  \n',
		'Rscript ${module_path}/missing_genotype_heterogygosity.R \\\n',
		'\t${outpath4step}/${name}.plink_QC1.imiss \\\n',
		'\t${outpath4step}/${name}.plink_QC1.het \\\n',
		'\t${outpath4step}/${name}.plink.imiss.het.pdf \\\n',
		'\t${outpath4step}/${name}.plink.imiss.het.remove \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC1 \\\n',
		'\t--remove ${outpath4step}/${name}.plink.imiss.het.remove \\\n',
		'\t--make-bed --out ${outpath4step}/${name}.plink_QC2 \n',

		'##QC3: Remove SNPs missing genotype rate >1%\n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC2  \\\n',
		'\t--missing --out ${outpath4step}/${name}.plink_QC2 \n',
		'awk \'NR < 2 { next } $5 > 0.01\' ${outpath4step}/${name}.plink_QC2.lmiss \\\n',
		'\t> ${outpath4step}/${name}.plink_QC2.lmiss.exclude  \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC2   \\\n',
		'\t--exclude ${outpath4step}/${name}.plink_QC2.lmiss.exclude   \\\n',
		'\t--make-bed --out ${outpath4step}/${name}.plink_QC3 \n\n',

		##QC4: Remove SNPs with differential missing genotype rate between cases and controls (+-5%)
		#plink --bfile all.plink_QC3 --test-missing --out all.plink_QC3
		#awk '$3-$4 > 0.05 || $3-$4 < -0.05 {print $2}' all.plink_QC3.missing > all.plink_QC3.missing.exclude
		#plink --bfile all.plink_QC3 --exclude all.plink_QC3.missing.exclude --make-bed --out all.plink_QC4

		'##QC5 :Remove SNPs deviated from HWE (10e-6)\n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC3 \\\n', 
		'\t--hardy --out ${outpath4step}/${name}.plink_QC3 \n',
		'awk \'$3=="UNAFF" && $9 <0.000001 {print $2}\' ${outpath4step}/${name}.plink_QC3.hwe  \\\n', 
		'\t> ${outpath4step}/${name}.plink_QC3.hwe.exclude \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC3 \\\n',
		'\t--exclude ${outpath4step}/${name}.plink_QC3.hwe.exclude \\\n', 
		'\t--make-bed --out ${outpath4step}/${name}.plink_QC4 \n',

		'##QC6: Principal component analysis (PCA) \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC4  \\\n',
		'\t--indep-pairwise 200 5 0.2  \\\n',
		'\t--out ${outpath4step}/${name}.plink_QC4  \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC4   \\\n',
		'\t--extract ${outpath4step}/${name}.plink_QC4.prune.in     \\\n',
		'\t--make-bed --out ${outpath4step}/${name}.plink_QC4_indep \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC4_indep \\\n',
		'\t--pca --out ${outpath4step}/${name}.plink_QC4_indep \n',
		'##Plot\n',
		'Rscript ${module_path}/PCA.R ${outpath4step}/${name}.plink_QC4_indep.eigenvec \\\n',
		'\t${outpath4step}/${name}.plink_QC4_indep.pca.pdf \\\n',
		'\t${outpath4step}/${name}.plink_QC4_indep.remove \n\n',
		'${plink} --bfile ${outpath4step}/${name}.plink_QC4  \\\n',
		'\t--remove ${outpath4step}/${name}.plink_QC4_indep.remove \\\n',
		'\t--make-bed --out ${outpath4step}/${name}.plink_clean \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_clean \\\n',
		'\t--pca --out ${outpath4step}/${name}.plink_clean \\\n',
		'##Plot\n',
		'Rscript ${module_path}/PCA1.R ${outpath4step}/${name}.plink_clean.eigenvec \\\n',
		'\t${outpath4step}/${name}.plink_clean.pca.pdf \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_clean \\\n',
		'\t--freq --out ${outpath4step}/${name}.plink_clean\n',
		'##Plot\n',
		'Rscript  ${module_path}/Plot_alleleFQ.R ${outpath4step}/${name}.plink_clean.frq \\\n',
		'\t${outpath4step}/${name}.plink_clean.maf.pdf \n',

		'###TDT\n',
		'###1)select common variations with MAF>1% \n',
		'${plink} --bfile ${outpath4step}/${name}.plink_clean \\\n', 
		'\t--maf 0.01 --make-bed \\\n',
		'\t--out ${outpath4step}/${name}.plink_clean_common_variants \n\n',
	
		'###2)TDT\n',
		'${plink} --bfile ${outpath4step}/${name}.plink_clean_common_variants \\\n',
		'\t--tdt --ci 0.95 --out ${outpath4step}/${name}.plink_clean_common_variants \n\n',
		'perl ${module_path}/sort_TDT.pl ${outpath4step}/${name}.plink_clean_common_variants.tdt \\\n',
		'\t>${outpath4step}/${name}.plink_clean_common_variants.tdt.sorted \n',
		'head -500 ${outpath4step}/${name}.plink_clean_common_variants.tdt.sorted \\\n', 
		'\t>${outpath4step}/${name}.plink_clean_common_variants.tdt.sort_top500 \n\n',
		'## Plot\n',
		'Rscript ${module_path}/Plot.manhattan_qq.R \\\n',
		'\t${outpath4step}/${name}.plink_clean_common_variants.tdt \\\n',
		'\t${outpath4step}/${name}.plink_clean_common_variants.tdt.manhattan.jpg \\\n',
		'\t${outpath4step}/${name}.plink_clean_common_variants.tdt.QQ.jpg\n\n']
	return commands

if __name__ == '__main__':
	#requirements preparation
	soft_opt=["GATK","plink","vcftools"]
	param_opt=["DeNovo"]
	ref_opt=["ref","dbsnp","ped"]
	in_opt=['inputtarget']
        out_opt=["outtarget"]
	inhouse_opt=["module_path"]
	opt=[soft_opt,param_opt,ref_opt,inhouse_opt,in_opt,out_opt]
	step='TDT'
	commands=Define_TDT_Output()
        ### Basic info, software, parameter,reference
        Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands)

