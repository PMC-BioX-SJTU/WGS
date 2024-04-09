#name=F1001
outpath4proj=/home/zhujh/Projects
for name in F1003  F264  F308  F325  F335  F349  F364  F375  F385  F398  F409  F427  F439  F455  F563  F614  F643  F750 F203   F271  F310  F327  F336  F352  F366  F379  F387  F399  F415  F434  F441  F456  F567  F617  F647 F213   F273  F311  F328  F342  F354  F370  F380  F388  F403  F421  F435  F443  F457  F584  F618  F712 F225   F280  F314  F329  F343  F358  F371  F382  F393  F404  F422  F436  F444  F560  F595  F623  F727 F1001            F258   F283  F319  F330  F347  F362  F373  F383  F395  F407  F423  F437  F450  F561  F612  F636  F743  F1002            F260   F304  F320  F334  F348  F363  F374  F384  F396  F408  F425  F438  F454  F562  F613  F641  F745  
do
outpath4step=/home/zhujh/Projects/StructVariationCalling/${name}

cat ${outpath4step}/BreakDancerCalling/${name}.DEL.filtered.merged.tab \
	${outpath4step}/CNVnatorCalling/${name}.DEL.filtered.merged.tab \
	${outpath4step}/DellyCalling/${name}.DEL.filtered.merged.tab >${outpath4step}/${name}.DEL.tab

cat ${outpath4step}/CNVnatorCalling/${name}.DUP.filtered.merged.tab \
	${outpath4step}/DellyCalling/${name}.DUP.filtered.merged.tab >${outpath4step}/${name}.DUP.tab

cat ${outpath4step}/DellyCalling/${name}.INV.filtered.merged.tab \
	${outpath4step}/BreakDancerCalling/${name}.INV.filtered.merged.tab >${outpath4step}/${name}.INV.tab

perl ~/Codes/Pipelines/StructVariationCalling/MergeBySample.pl -n ${name} -o ${outpath4step}
done
