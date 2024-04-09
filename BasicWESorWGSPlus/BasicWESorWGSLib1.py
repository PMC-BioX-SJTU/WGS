#/usr/bin/python
import ConfigParser
import string
import os
import argparse
import sys
import re
cf = ConfigParser.ConfigParser()
class Read_ConfigFile:
	def __init__(self,FileName):
		cf.read(FileName)
	
	def Get_Options(self,sec_name):
		info=cf.options(sec_name)
		return info

	def Get_OptionItems(self,sec_name,opt_name):
		item=cf.get(sec_name, opt_name)
		return item
	
	def Get_OptionLists(self,sec,opt):
		lists=[]
		for k in range(0,len(opt)):
			lists.append(cf.get(sec, opt[k]))
		return lists
	
	def Get_OptionLists4shell(self,sec,opt,i):
		setting=['###'+sec[i]+'\n']
		for k in range(0,len(opt[i])):
			setting.append(opt[i][k]+'='+cf.get(sec[i], opt[i][k])+'\n')
		setting.append('\n')
		return setting


class Obtain_SampleInfo(Read_ConfigFile):
	def __init__(self,sample): 
		self.sample=sample
		Read_ConfigFile.__init__(self,self.sample)
	
	def Get_Info(self,sec_name):
		info=self.Get_Options(sec_name)
		Length=len(info)
		return(info,Length)
        
	def Check_Length(self):
		(self.SampleName,len1)=self.Get_Info("Sample.name")
		(self.SampleData,len2)=self.Get_Info("Sample.data")
		(self.SampleDir, len3)=self.Get_Info("Sample.dir")
		if (len1 != len2 or len1 != len3 ): ### check sample.ini input
        	        print "Sample count is not matched with fq count!!!"
                	sys.exit(1)

        def Define_Names(self):
                self.names=self.Get_OptionLists("Sample.name",self.SampleName)
	
	def Define_Dirs(self):
		self.dirs=self.Get_OptionLists("Sample.dir",self.SampleDir)

	def Define_FQs(self):
		self.fqs=self.Get_OptionLists("Sample.data",self.SampleData)

class Input_Output_Settings(Read_ConfigFile):
	def __init__(self, outpath4proj,module,step,sample,config):
        	self.outpath4proj = outpath4proj
        	self.sample = sample 
		self.step = step
		self.module=module
		self.famGroup=""
		self.fq=""
		self.chrGroup=""
		self.chrfamGroup=""
		self.allGroup=""
		Read_ConfigFile.__init__(self,config)
	
	def Process_SampleInfo(self):
		if (re.findall('fq.gz|fastq.gz',self.sample)):
			(self.name,self.fq)=self.sample.split('|')
		elif(re.search('^Fam',self.sample)):
			(self.name,self.mem)=self.sample.split('|')
		elif(re.search('^chr',self.sample)):
			(self.name,self.mem)=self.sample.split('|')
		elif(re.search('^all\|',self.sample)):
                        (self.name,self.mem)=self.sample.split('|')
		else:
			self.name=self.sample
	
	def Define_Outdirs(self):
		#if(re.search('^FamGroup',self.name)):
		#self.outpath4step=self.outpath4proj+'/'+self.step
		#else:
		self.outpath4step=self.outpath4proj+'/'+self.module+'/'+self.name+'/'+self.step
		if  not os.path.exists(self.outpath4step):
			os.makedirs(self.outpath4step)
	
	def Define_Softwares(self,sec,opt):
		self.software=self.Get_OptionLists4shell(sec,opt,0)

	def Define_Parameters(self,sec,opt):
		self.parameter=self.Get_OptionLists4shell(sec,opt,1)	
	
	def Define_References(self,sec,opt):
		self.reference=self.Get_OptionLists4shell(sec,opt,2)
	
	def Define_InHouse_Scripts(self,sec,opt):
		self.InHouse_Scripts=self.Get_OptionLists4shell(sec,opt,3)

	def Define_InputTargets(self):
		self.InputTarget=self.Get_OptionItems(self.step + ".input","input.target")
		if (self.InputTarget=="from_sample.ini"):
			self.InputTarget=self.fq
			(R1,R2,self.InputPath)=self.fq.split(":")
			R1_lst=R1.split()
			R2_lst=R2.split()
			if(len(R1_lst) != len(R2_lst)):
				print name,"left Read cannot match right Read completely!!"
				sys.exit(1)
			else:
				self.InputTarget=[['${InputPath}/'+ x for x in R1_lst],['${InputPath}/'+ y for y in R2_lst]]
		elif(self.InputTarget=="from_fam.ped"):
			famMemLst=self.mem.split(":")
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
			self.InputTarget=['${InputPath}/BasicWESorWGS/' + x + '/MappingAndCleanUp/' + x + '.recal.bam'  for x in famMemLst]
		elif(self.InputTarget=="from_fam.group"):
			famMemLst=self.mem.split(":")
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
			Target1=['${InputPath}/BasicWESorWGS/'+ x + '/MappingAndCleanUp/' + x + '.recal.bam'  for x in famMemLst]
			Target2=['${InputPath}/BasicWESorWGS/'+ x + '/GermlineVariantCalling/'+ x + '.raw.gvcf.gz'  for x in famMemLst]
			tar=' '.join(sum([Target1,Target2],[]))
			tar1=' -I '.join(Target1)
			tar2=' -V '.join(Target2)
			self.InputTarget=[tar,tar1,tar2]
			self.famGroup='from_fam.group'
		elif(self.InputTarget=="from_chr.group"):
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")			
			chrLst=self.mem.split(":")
			tar=' -V '.join(chrLst)
			tar1=' '.join(chrLst)
			self.InputTarget=[tar,tar1]
			self.chrGroup='from_chr.group'
		elif(self.InputTarget=="from_chrfam.group"):
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
			chrfamLst=self.mem.split(":")
			tar=' --variant '.join(chrfamLst)
			tar1=' '.join(chrfamLst)
			self.InputTarget=[tar,tar1]
			self.chrfamGroup='from_chrfam.group'
		elif(self.InputTarget == 'from_all.group'):
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
			allLst=self.mem.split(":")
			tar1=['${InputPath}/BasicWESorWGS/'+ x + '/QualityControlBam/' + x + '.Q20Q30.final.txt'  for x in allLst]
			tar2=['${InputPath}/BasicWESorWGS/'+ x + '/QualityControlBam/' + x + '.Alig.final.txt'  for x in allLst]
			tar3=['${InputPath}/BasicWESorWGS/'+ x + '/QualityControlBam/' + x + '.CrossCont.final.txt'  for x in allLst]
			tar4=['${InputPath}/BasicWESorWGS/'+ x + '/QualityControlBam/' + x + '.Hs.final.txt'  for x in allLst]
			tar=' '.join(sum([tar1,tar2,tar3,tar4],[]))
			self.InputTarget=[tar,tar1,tar2,tar3,tar4]
			self.allGroup='from_all.group'
		else:
			self.InputPath=self.Get_OptionItems(self.step + ".input","input.path")
			out=self.InputTarget.split(":")
			self.InputTarget=['${InputPath}/'+ x for x in out]

	def Define_OutputTargets(self):
		self.OutputTarget=self.Get_OptionItems(self.step + ".output","output.target")
		self.OutputPath=self.Get_OptionItems(self.step + ".output","output.path")
		out=self.OutputTarget.split(":")
		self.OutputTarget=['${OutputPath}/' + x for x in out]
	 	
	def Define_Basic(self):
		self.basic=['###Basic\n','outpath4proj='+self.outpath4proj+'\n']
		self.basic.append('name='+self.name+'\n')
		self.basic.append('outpath4step='+self.outpath4step+'\n')
		self.basic.append('InputPath='+self.InputPath+'\n')
		self.basic.append('OutputPath='+self.OutputPath+'\n')
		if(self.fq != ""):
			self.basic.append('fq1=('+ ' '.join(self.InputTarget[0]) +')\n')
			self.basic.append('fq2=('+ ' '.join(self.InputTarget[1]) +')\n')
			self.basic.append('inputtarget=(${fq1[@]} ${fq2[@]})\n')
		elif(self.famGroup == 'from_fam.group'):
			self.basic.append('inputtarget=('+self.InputTarget[0]+')\n')
			self.basic.append('bams=('+self.InputTarget[1]+')\n')
			self.basic.append('vcfs=('+self.InputTarget[2]+')\n')
		elif(self.chrGroup == 'from_chr.group' ):
			self.basic.append('inputtarget=('+self.InputTarget[1]+')\n')
			self.basic.append('gvcf=('+self.InputTarget[0]+')\n')
		elif(self.chrfamGroup == 'from_chrfam.group'):
			self.basic.append('inputtarget=('+self.InputTarget[1]+')\n')
			self.basic.append('vcf=('+self.InputTarget[0]+')\n')
		elif(self.allGroup == 'from_all.group'):
			self.basic.append('inputtarget=('+self.InputTarget[0]+')\n')
			self.basic.append('Q20Q30=('+' '.join(self.InputTarget[1])+')\n')
			self.basic.append('Alig=('+' '.join(self.InputTarget[2])+')\n')
			self.basic.append('CrossCont=('+' '.join(self.InputTarget[3])+')\n')
			self.basic.append('Hs=('+' '.join(self.InputTarget[4])+')\n')
		else:
			self.basic.append('inputtarget=('+' '.join(self.InputTarget)+')\n')
		self.basic.append('outtarget=('+' '.join(self.OutputTarget)+')\n\n')

def Define_Requirements(outpath4proj,module,step,sample,config,sec,opt):	
	requirements=Input_Output_Settings(outpath4proj,module,step,sample,config)
	requirements.Process_SampleInfo()
	requirements.Define_Outdirs()
	requirements.Define_Softwares(sec,opt)
	requirements.Define_Parameters(sec,opt)
	requirements.Define_References(sec,opt)
	requirements.Define_InHouse_Scripts(sec,opt)
	requirements.Define_InputTargets()
	requirements.Define_OutputTargets()
	requirements.Define_Basic()
	return requirements

def Check_InputOutputTarget(step):
	log1=['for file in ${inputtarget[@]}\n',
		'do\n',
        	'if [ ! -f "${file}" ]; then\n',
        	'\techo -e the target file ${file} cannot be generated successfully!!\n'
        	'\texit 0\n',
                'fi\n',
                'done\n',
		'echo -e The Step '+step+' started at `date` "\\n" >>${outpath4step}/run.log;\n\n']
	
        log2=[	'echo -e The Step '+step+' ended at `date` "\\n" >>${outpath4step}/run.log;\n',
		'for file in ${outtarget[@]}\n',
        	'do\n',
        	'if [ ! -f "${file}" ]; then\n',
        	'\techo -e the target file ${file} cannot bed generated successfully!!\n',
        	'\texit 0\n',
        	'fi\n',
        	'done\n',
        	'\tmv ${outpath4step}/run.log ${outpath4step}/'+step+'.log\n\n']
	return(log1,log2)

def Check_SampleRawFq(name,fq):
	(r1,r2,inputdir)=fq.split(":")
	raw_r1=r1.split()
	raw_r2=r2.split()
	if(len(raw_r1) != len(raw_r2)):
		print name,"left Read cannot match right Read completely!!"
		sys.exit(1)
	else:
		return(raw_r1,raw_r2,r1,r2,inputdir)

def Check_Exists(file):
        for i in range(0,len(file)):
                if not os.path.exists(file[i]):
                        print "We cannot find the file "+file[i]+"!!"
                        sys.exit(1)

def Write_makefile(outfile,outpath,tar_out,ext_out,actStep):
	step_tar=['outpath=',outpath,'\n',
		'\nstep_tar=',tar_out,
		'\n.PHONY:all\n',
		'all:$(step_tar)\n']
	step0=['define steps\n']
	step1=['$2:'+'$(outpath)/$1/'+actStep+'/'+actStep+'.sh\n',
		'\tsh $(outpath)/$1/'+actStep+'/'+actStep+'.sh\n']
	step2=['endef\n']
	output=sum([step_tar,step0,step1,step2,ext_out],[])
	f=open(outfile+actStep+'.mk',"w")
	f.writelines(output)
	f.close

def Write_eachSteps(outpath4proj,module,step,sample,config,opt,commands): 
        ### Basic info, software, parameter,reference
        sec=['software',step+'.parameter','reference','inhouse_scripts']
	requirements=Define_Requirements(outpath4proj,module,step,sample,config,sec,opt)
        require=requirements
        require_setting=[]
        require_setting=sum([require.basic,require.software,require.parameter,require.reference,require.InHouse_Scripts],[])
        (log1,log2)=Check_InputOutputTarget(step)
	### write makefile
	output=sum([require_setting,log1,commands,log2],[])
	makefilename=require.outpath4step+'/'+step+'.sh'
	f=open(makefilename,"w")
	output=sum([require_setting,log1,commands,log2],[])
	f.writelines(output)
	f.close

if __name__ == '__main__':
	import ConfigParser
	import string
	import os
	import argparse
	import sys
	import re
	cf = ConfigParser.ConfigParser()

	parser = argparse.ArgumentParser(description='the Step for Mapping and Mark Duplication reads!')
	parser.add_argument('-c', required=True, dest='config',action='store',help='config file for the module')
	parser.add_argument('-s', required=True, dest='sample',action='store',help='sample info')
	parser.add_argument('-o', required=True, dest='outpath4proj', action='store',help='outpath for the project')

	if len(sys.argv) <= 3:
        	parser.print_help()
        	sys.exit(1)
	else:
		args = parser.parse_args()
	config = args.config
	outpath4proj = args.outpath4proj
	sample = args.sample
	
	#requirements preparation
	actSteps=[]
	ConfigInfo=Read_ConfigFile(config)
	module_path=ConfigInfo.Get_OptionItems('software','module_path')
	steps=ConfigInfo.Get_OptionItems('actSteps','StepstoRun')
	actSteps=steps.split(',')
	waitSteps=ConfigInfo.Get_OptionItems('actSteps','StepstoWait')
	print actSteps, waitSteps
