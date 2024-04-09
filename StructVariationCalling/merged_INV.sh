###Basic
outpath4proj=/home/zhujh/Projects
name=all
outpath4step=/home/zhujh/Projects/StructVariationCalling
InputPath=${outpath4proj}/StructVariationCalling
OutputPath=${outpath4step}
inputtarget=(${InputPath}/F203/DellyCalling/F203.INV.bcf \
${InputPath}/F213/DellyCalling/F213.INV.bcf \
${InputPath}/F225/DellyCalling/F225.INV.bcf \
${InputPath}/F258/DellyCalling/F258.INV.bcf \
${InputPath}/F260/DellyCalling/F260.INV.bcf \
${InputPath}/F264/DellyCalling/F264.INV.bcf \
${InputPath}/F271/DellyCalling/F271.INV.bcf \
${InputPath}/F273/DellyCalling/F273.INV.bcf \
${InputPath}/F280/DellyCalling/F280.INV.bcf \
${InputPath}/F283/DellyCalling/F283.INV.bcf \
${InputPath}/F304/DellyCalling/F304.INV.bcf \
${InputPath}/F308/DellyCalling/F308.INV.bcf \
${InputPath}/F310/DellyCalling/F310.INV.bcf \
${InputPath}/F311/DellyCalling/F311.INV.bcf \
${InputPath}/F314/DellyCalling/F314.INV.bcf \
${InputPath}/F319/DellyCalling/F319.INV.bcf \
${InputPath}/F320/DellyCalling/F320.INV.bcf \
${InputPath}/F325/DellyCalling/F325.INV.bcf \
${InputPath}/F327/DellyCalling/F327.INV.bcf \
${InputPath}/F328/DellyCalling/F328.INV.bcf \
${InputPath}/F329/DellyCalling/F329.INV.bcf \
${InputPath}/F330/DellyCalling/F330.INV.bcf \
${InputPath}/F334/DellyCalling/F334.INV.bcf \
${InputPath}/F335/DellyCalling/F335.INV.bcf \
${InputPath}/F336/DellyCalling/F336.INV.bcf \
${InputPath}/F342/DellyCalling/F342.INV.bcf \
${InputPath}/F343/DellyCalling/F343.INV.bcf \
${InputPath}/F347/DellyCalling/F347.INV.bcf \
${InputPath}/F348/DellyCalling/F348.INV.bcf \
${InputPath}/F349/DellyCalling/F349.INV.bcf \
${InputPath}/F352/DellyCalling/F352.INV.bcf \
${InputPath}/F354/DellyCalling/F354.INV.bcf \
${InputPath}/F358/DellyCalling/F358.INV.bcf \
${InputPath}/F362/DellyCalling/F362.INV.bcf \
${InputPath}/F363/DellyCalling/F363.INV.bcf \
${InputPath}/F364/DellyCalling/F364.INV.bcf \
${InputPath}/F366/DellyCalling/F366.INV.bcf \
${InputPath}/F370/DellyCalling/F370.INV.bcf \
${InputPath}/F371/DellyCalling/F371.INV.bcf \
${InputPath}/F373/DellyCalling/F373.INV.bcf \
${InputPath}/F374/DellyCalling/F374.INV.bcf \
${InputPath}/F375/DellyCalling/F375.INV.bcf \
${InputPath}/F379/DellyCalling/F379.INV.bcf \
${InputPath}/F380/DellyCalling/F380.INV.bcf \
${InputPath}/F382/DellyCalling/F382.INV.bcf \
${InputPath}/F383/DellyCalling/F383.INV.bcf \
${InputPath}/F384/DellyCalling/F384.INV.bcf \
${InputPath}/F385/DellyCalling/F385.INV.bcf \
${InputPath}/F387/DellyCalling/F387.INV.bcf \
${InputPath}/F388/DellyCalling/F388.INV.bcf \
${InputPath}/F393/DellyCalling/F393.INV.bcf \
${InputPath}/F395/DellyCalling/F395.INV.bcf \
${InputPath}/F396/DellyCalling/F396.INV.bcf \
${InputPath}/F398/DellyCalling/F398.INV.bcf \
${InputPath}/F399/DellyCalling/F399.INV.bcf \
${InputPath}/F403/DellyCalling/F403.INV.bcf \
${InputPath}/F404/DellyCalling/F404.INV.bcf \
${InputPath}/F407/DellyCalling/F407.INV.bcf \
${InputPath}/F408/DellyCalling/F408.INV.bcf \
${InputPath}/F409/DellyCalling/F409.INV.bcf \
${InputPath}/F415/DellyCalling/F415.INV.bcf \
${InputPath}/F421/DellyCalling/F421.INV.bcf \
${InputPath}/F422/DellyCalling/F422.INV.bcf \
${InputPath}/F423/DellyCalling/F423.INV.bcf \
${InputPath}/F425/DellyCalling/F425.INV.bcf \
${InputPath}/F427/DellyCalling/F427.INV.bcf \
${InputPath}/F434/DellyCalling/F434.INV.bcf \
${InputPath}/F435/DellyCalling/F435.INV.bcf \
${InputPath}/F436/DellyCalling/F436.INV.bcf \
${InputPath}/F437/DellyCalling/F437.INV.bcf \
${InputPath}/F438/DellyCalling/F438.INV.bcf \
${InputPath}/F439/DellyCalling/F439.INV.bcf \
${InputPath}/F441/DellyCalling/F441.INV.bcf \
${InputPath}/F443/DellyCalling/F443.INV.bcf \
${InputPath}/F444/DellyCalling/F444.INV.bcf \
${InputPath}/F450/DellyCalling/F450.INV.bcf \
${InputPath}/F454/DellyCalling/F454.INV.bcf \
${InputPath}/F455/DellyCalling/F455.INV.bcf \
${InputPath}/F456/DellyCalling/F456.INV.bcf \
${InputPath}/F457/DellyCalling/F457.INV.bcf \
${InputPath}/F560/DellyCalling/F560.INV.bcf \
${InputPath}/F561/DellyCalling/F561.INV.bcf \
${InputPath}/F562/DellyCalling/F562.INV.bcf \
${InputPath}/F563/DellyCalling/F563.INV.bcf \
${InputPath}/F567/DellyCalling/F567.INV.bcf \
${InputPath}/F584/DellyCalling/F584.INV.bcf \
${InputPath}/F595/DellyCalling/F595.INV.bcf \
${InputPath}/F612/DellyCalling/F612.INV.bcf \
${InputPath}/F613/DellyCalling/F613.INV.bcf \
${InputPath}/F614/DellyCalling/F614.INV.bcf \
${InputPath}/F617/DellyCalling/F617.INV.bcf \
${InputPath}/F618/DellyCalling/F618.INV.bcf \
${InputPath}/F623/DellyCalling/F623.INV.bcf \
${InputPath}/F636/DellyCalling/F636.INV.bcf \
${InputPath}/F641/DellyCalling/F641.INV.bcf \
${InputPath}/F643/DellyCalling/F643.INV.bcf \
${InputPath}/F647/DellyCalling/F647.INV.bcf \
${InputPath}/F712/DellyCalling/F712.INV.bcf \
${InputPath}/F727/DellyCalling/F727.INV.bcf \
${InputPath}/F743/DellyCalling/F743.INV.bcf \
${InputPath}/F745/DellyCalling/F745.INV.bcf \
${InputPath}/F750/DellyCalling/F750.INV.bcf \
${InputPath}/F1001/DellyCalling/F1001.INV.bcf \
${InputPath}/F1002/DellyCalling/F1002.INV.bcf \
${InputPath}/F1003/DellyCalling/F1003.INV.bcf )


outtarget=(${OutputPath}/${name}.INV.filtered.merged.tab ${OutputPath}/${name}.DUP.filtered.merged.tab ${OutputPath}/${name}.INV.filtered.merged.tab)

###software
delly=/home/zhujh/software/delly/src/delly
bcftools=/home/zhujh/software/delly/src/bcftools/bcftools

###DellyCalling.parameter
Filtering_para="centelBuffer=1000:gapOverlap=0.5:gapsBuffer=600"
sv_type="INV DUP INV"

###reference
ref=/home/zhujh/reference/ucsc.hg19.fasta
hg19_excl=/home/zhujh/software/delly/excludeTemplates/human.hg19.excl.tsv
gaps=/home/zhujh/Codes/Pipelines/StructVariationCalling/hg19_gap.txt
centel=/home/zhujh/Codes/Pipelines/StructVariationCalling/hg19_cen_tel.txt

###inhouse_scripts
module_path=/home/zhujh/Codes/Pipelines/StructVariationCalling

echo ${inputtarget[@]}
${delly} merge -t INV -o ${outpath4step}/all.INV.bcf ${inputtarget[@]}

