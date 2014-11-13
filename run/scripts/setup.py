#!/usr/bin/python
import os



def mkMainDirs(TargetOrtho_path):
	if not os.path.exists('%s/run/output' %TargetOrtho_path):os.makedirs('%s/run/output' %TargetOrtho_path)
	if not os.path.exists('%s/run/fimo_out' %TargetOrtho_path):os.makedirs('%s/run/fimo_out' %TargetOrtho_path)
	if not os.path.exists('%s/run/associate_genes_out' %TargetOrtho_path):os.makedirs('%s/run/associate_genes_out' %TargetOrtho_path)
	if not os.path.exists('%s/run/input' %TargetOrtho_path):os.makedirs('%s/run/input' %TargetOrtho_path)
	if not os.path.exists('%s/run/sqlite_tmp_files' %TargetOrtho_path):os.makedirs('%s/run/sqlite_tmp_files' %TargetOrtho_path)

def mkOutputDir(jobID,output_dir,PSSMnamesList,PSSMinfile,queryFile):
	print 'running setup.mkOutputDir()'	
	os.system("rm -r %s" %(output_dir))	
	os.makedirs(output_dir)
	for n in range(len(PSSMnamesList)):
		print PSSMnamesList[n]
		os.makedirs("%s/%s" %(output_dir,PSSMnamesList[n].split()[0]))
		os.makedirs("%s/%s/Results_plots" %(output_dir,PSSMnamesList[n].split()[0]))
		os.makedirs("%s/%s/Results_text" %(output_dir,PSSMnamesList[n].split()[0]))
		os.makedirs("%s/%s/Results_html" %(output_dir,PSSMnamesList[n].split()[0]))
		os.makedirs("%s/%s/inputFiles" %(output_dir,PSSMnamesList[n].split()[0]))
		#cp inputFile to inputFiles directory
		os.system("cp %s %s/%s/inputFiles" %(PSSMinfile,output_dir,PSSMnamesList[n].split()[0]))
		os.system("cp %s %s/%s/inputFiles" %(queryFile,output_dir,PSSMnamesList[n].split()[0]))
	
def writeInputSum(jobID,jobTag,paramDic,PSSMnamesList,output_dir):
	weightDic={'scale_A':'scale factor for average species conserved site scor',
		  'scale_B':'scale factor for average species region score',
		  'scale_C':'scale factor for average species region site count',
		  'scale_D':'scale factor for average species average gene score',
		  'scale_E':'scale factor for average species gene site count',
		  'scale_F':'scale factor for conservation',
		  'scale_G':'scale factor for average species offset variance'}
			  
	print 'running setup.writeInputSum()'
	inputSum=file('%s/%s_inputSummary.txt' %(output_dir,jobID), 'w')
	inputSum.write('jobID: %s\njob name: %s\n' %(jobID,jobTag))
	inputSum.write('motifs run:\t')
	for name in PSSMnamesList:
		inputSum.write('%s;' %(name))
	inputSum.write('\n')
	paramList=[]
	for k,v in paramDic.iteritems():
		paramList.append(k)
	paramList.sort()
	for k in paramList:
		print k
		if 'order' in k or 'dist' in k:
			if len(PSSMnamesList)>1:inputSum.write("%s:\t%s\n" %(k,paramDic[k]))
		elif 'total_hits' in k:
			total_hits=str(paramDic[k].split()[:len(PSSMnamesList)])[1:-1]
			inputSum.write("%s:\t%s\n" %(k,total_hits))
		elif 'd_' in k:
			if paramDic['species']=='melanogaster5':
				inputSum.write("%s:\t%s\n" %(k,paramDic[k]))
		elif 'scale' in k:
				inputSum.write("%s:\t%s\n" %(weightDic[k],paramDic[k]))
			
		else:inputSum.write("%s:\t%s\n" %(k,paramDic[k]))
	
