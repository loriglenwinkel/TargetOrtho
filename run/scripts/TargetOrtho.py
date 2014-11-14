#!/usr/bin/python
import os
import sys
import datetime
from optparse import OptionParser
import scipy.stats
import setup
import validate_infile
import threadScripts
#command line version (or run through galaxy at alice.cumc.columbia.edu:8080)    
parser = OptionParser()
parser.add_option( '-f', '--file', dest='infile', help='-f filepath/motif_input.txt' )
parser.add_option('-o', '--output_dir',dest='output_dir',help='-o my_results_dir, default output is to www.alice.cumc.columabie.du/TargetOrtho/jobID. A link to this webpage is printed at the end of the job')
parser.add_option( '-q', '--file2', dest='infile2', help='-q filepath.queryList.txt' )
parser.add_option( '-x', '--max_upstream', dest='max_upstream', help='-x None #maximum distnace to search for motif match upstream of ATG')
parser.add_option( '-i', '--max_downstream', dest='max_downstream', help='-i 500')
parser.add_option( '-v', '--max_offset_variance', dest='max_offset_variance', help='-v 1')
parser.add_option( '-e', '--filter_exons', dest='filter_exons', help='-e False')
parser.add_option( '-P', '--region_overlap', dest='region_overlap', help='-p True or -p False #True if motif genomic region should be shared in all orthologs. i.e. all motif matches in exon, intron, intergenic-downstream, or intergenic-upstream')
parser.add_option( '-k', '--top_ranked_per_gene', dest='top_ranked_per_gene', type="int",help='-k 1')
parser.add_option( '-d', '--distance_between', dest='bp_dist', type="int",help='-d 1000 or -d None#maximum distance in nucleotides between motif matches (optional when searching for combination of 2 or more motifs matches')
parser.add_option( '-y', '--ordered', dest='motif_order', help='True or False: -o True means motif matches must be same order as input file on DNA')
parser.add_option( '-c', '--config' , dest='configFilePath', help='-c /TargetOrtho_path/run/scripts/config.txt')
parser.add_option( '-z', '--genes_between' , dest='genes_between', help='number of genes allowed between motif match and associategene: -z 0')
parser.add_option( '-Z', '--max_dist_if_genes_between' , dest='max_dist_if_genes_between', help='absolute nucleotide distance to search if genes are present between site and associated gene -z 2000')
parser.add_option('-M','--gene_name',dest='gene_name',help='-a unc-3')
parser.add_option('-w','--QueryOnly',dest='QueryOnly',help='-w True to only run TargetOrtho on query list genes')
parser.add_option('-p',dest='p_cutoff',help='--p_cutoff 1e-6')
parser.add_option('-g','--jobID',dest='jobID',help='-g 123')
parser.add_option('-J','--jobTag',dest='jobTag',help='-J ASE_known_genes')
parser.add_option('-a','--scale_A',dest='scale_A',help='-a 1 ..factor to scale the avg species site score before ranking final results')
parser.add_option('-b','--scale_B',dest='scale_B',help='-b 1 ..factor to scale the avg species region score before ranking final results')
parser.add_option('-H','--scale_C',dest='scale_C',help='-c 1 ..factor to scale the avg species region site count before ranking final results')
parser.add_option('-I','--scale_D',dest='scale_D',help='-d 1 ..factor to scale the avg species total gene score before ranking final results')
parser.add_option('-K','--scale_E',dest='scale_E',help='-e 1 ..factor to scale the avg species gene site counts  before ranking final results')
parser.add_option('-L','--scale_F',dest='scale_F',help='-f 1 ..factor to scale the conservation score before ranking final results')
parser.add_option('-N','--scale_G',dest='scale_G',help='-g 1 ..factor to scale the avg species offset variance  before ranking final results')
parser.add_option('-O','--species',dest='species',help='-O melagnogaster or -O nematodes -O nem_mel')
( options, args ) = parser.parse_args()

print options.jobID,'jobID in TargetOrtho.py'

def fetchConfig(configFile):
    infile=file(configFile,'r')
    lines=infile.readlines()
    configDic={}
    for line in lines:
        sep=line.split('\'')
        if len(sep)>1:configDic[sep[0]]=sep[1]
    return configDic
    
#get TargetOrtho path
p=os.path.abspath(__file__).split('/')[1:-3]
s=''
for i in p:
    s=s +'/%s'%i
TargetOrtho_path=s
print TargetOrtho_path,'TargetOrtho_path'

#find confit.txt
#try to use the config.txt file given by the user, next use the default config file in scripts/
try:configDic=fetchConfig(options.configFilePath)
except:configDic=fetchConfig('%s/run/scripts/config.txt' %TargetOrtho_path)

#define species list
#species1 is reference species
speciesList=['c_eleg','c_brig','c_bren','c_rema','c_japo']
speciesList2=['d_mel','d_sec','d_sim','d_yak','d_ere']
if options.species=='melanogaster5':speciesList=speciesList2
speciesStr=''
for n in speciesList:
    speciesStr=speciesStr+ n + '-'
sys.path.append('%s/run/scripts' %(TargetOrtho_path))
gene_name=options.gene_name

def main():
    def assign_jobID():
        dt=datetime.datetime.now()
        if options.jobID==None:jobID="j%s%s%s%s%s%s" %(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)        
        else:jobID=options.jobID
        if options.jobTag==None:jobTag=jobID
        else:jobTag=options.jobTag
        print jobID,jobTag
        return jobID,jobTag

    def load_defaults():
        """unique to cmd line version"""
        print 'loading default dic: load_defaults()'

        defaults_dic={     
        'species':'nematodes',
        'c_eleg_genomeVersion':configDic["c_eleg_genome_version="],
        'c_brig_genomeVersion':configDic["c_brig_genome_version="],
        'c_bren_genomeVersion':configDic["c_bren_genome_version="],
        'c_rema_genomeVersion':configDic["c_rema_genome_version="],
        'c_japo_genomeVersion':configDic["c_japo_genome_version="],
        'd_mel_genomeVersion':configDic["d_mel_genome_version="],
        'd_sec_genomeVersion':configDic["d_sec_genome_version="],
        'd_sim_genomeVersion':configDic["d_sim_genome_version="],
        'd_yak_genomeVersion':configDic["d_yak_genome_version="],
        'd_ere_genomeVersion':configDic["d_ere_genome_version="],
        'QueryOnly':'False',
        'p_cutoff':'1e-4',
        'genes_between':20,
        'max_dist_if_genes_between':6000,
        'max_upstream': 6000,
        'max_downstream': 0,
        'filter_exons': 'False',
        'region_overlap':'True',
        'top_ranked_per_gene':1,
        'max_offset_variance':1.0,
        'bp_dist':6000,
        'motif_order':'False',
        'scale_A':1,
        'scale_B':1,
        'scale_C':1,
        'scale_D':1,
        'scale_E':1,
        'scale_F':1,
        'scale_G':1,
        }    
        return defaults_dic
        
    def load_cmd_line_params():
        """unique to cmd line version"""
        print 'running use_command_line_params()'
        cmd_line_params_dic={
        'species':options.species,
        'QueryOnly':options.QueryOnly,
        'motif_infile':options.infile,
        'QueryFile':options.infile2,
        'p_cutoff':options.p_cutoff,
        'max_upstream': options.max_upstream,
        'max_downstream': options.max_downstream,
        'filter_exons': options.filter_exons,
        'region_overlap':options.region_overlap,
        'top_ranked_per_gene':options.top_ranked_per_gene,
        'max_offset_variance':options.max_offset_variance,
        'scale_A':options.scale_A,
        'scale_B':options.scale_B,
        'scale_C':options.scale_C,
        'scale_D':options.scale_D,
        'scale_E':options.scale_E,
        'scale_F':options.scale_F,
        'scale_G':options.scale_G,
        'bp_dist':options.bp_dist,
        'motif_order':options.motif_order,
        'genes_between':options.genes_between,
        'max_dist_if_genes_between':options.max_dist_if_genes_between}
        #fill in missing command line parameters with default values
       
        
        #validate input motif file
        matrix_count,matrixNamesList=validate_infile.validate(options.infile,TargetOrtho_path,jobID)
        print 'filling in missing parameters using default values'
        print jobID,'jobID'
        print matrix_count,matrixNamesList
        
        #load parameters from command line or fill in default values. If run through galaxy, galaxy will pass parameters to command line options.
        defaults_dic=load_defaults()
        for k,v in defaults_dic.iteritems():            
            if k in cmd_line_params_dic.keys():
                if cmd_line_params_dic[k]==None:cmd_line_params_dic[k]=v#if cmd line  param is not used, set equal to default value, v
            else:cmd_line_params_dic[k]=v                        
        return cmd_line_params_dic,matrix_count,matrixNamesList
        
    def print_parameters(dict):
        print 'printing parameters'
        for k,v in dict.iteritems():
            print k,v,type(v)    
        
    def build_queues(matrix_count,speciesList,param_dic,matrixNamesList,jobID,gene_name):
        """sets up queues for threadScripts module to execute fimo,associate_genes.py, and findOrthoAssoc.py for each species, for each matrix in parallel"""
        queue1=[]
        queue2=[]
        queue3=[]
        queue4=[]
        queue5=[]        
        for species in speciesList:
            for n in range(matrix_count):
                print matrix_count,'matrix_count',n
                if options.p_cutoff!=None:
                    p=float(options.p_cutoff.split()[n])
                    assert p>0 and p<=1, 'P value must be between 0 and 1'
                    pthresh='--thresh %s' %(p)

                else:pthresh=''
                print "%s/run/scripts/fimo --text %s %s/run/input/%s_meme%s.txt %s/genomes/%s_%s/%s_%s.fa > %s/run/fimo_out/%s_%s_%s_fimo.out" %(TargetOrtho_path,pthresh,TargetOrtho_path,jobID,n+1,TargetOrtho_path,species,param_dic['%s_genomeVersion' %(species)],species,param_dic['%s_genomeVersion' %(species)],TargetOrtho_path,jobID,species,n+1)
                queue1.append("%s/run/scripts/fimo --text %s %s/run/input/%s_meme%s.txt %s/genomes/%s_%s/%s_%s.fa > %s/run/fimo_out/%s_%s_%s_fimo.out" %(TargetOrtho_path,pthresh,TargetOrtho_path,jobID,n+1,TargetOrtho_path,species,param_dic['%s_genomeVersion' %(species)],species,param_dic['%s_genomeVersion' %(species)],TargetOrtho_path,jobID,species,n+1))
                queue2.append("python %s/run/scripts/associate_genes.py -j %s -n %s -m %s  -c %s -v %s -z %s -Z %s -x %s -p %s -w %s" %(TargetOrtho_path,jobID,species,n+1,TargetOrtho_path,param_dic['%s_genomeVersion' %(species)],param_dic["genes_between"],param_dic["max_dist_if_genes_between"],param_dic['filter_exons'],param_dic['max_upstream'],param_dic['max_downstream'])                )
                queue3.append("python %s/run/scripts/findOrthoAssoc.py -j %s -n %s -m %s -v %s -g %s -r %s -f %s -Z %s -s %s -c %s" %(TargetOrtho_path,jobID,species,n+1,param_dic['%s_genomeVersion' %(species)],gene_name,param_dic['QueryOnly'],param_dic['QueryFile'],param_dic["max_dist_if_genes_between"],speciesStr[:-1],TargetOrtho_path))
        for n in range(matrix_count):
            queue4.append("python %s/run/scripts/conserved_ortho.py -j %s -m %s -s %s  -p %s -x %s -i %s -k %s -c %s" %(TargetOrtho_path,jobID,n+1,speciesStr[:-1],param_dic['region_overlap'],param_dic['max_upstream'],param_dic['max_downstream'],param_dic["top_ranked_per_gene"],TargetOrtho_path))
            queue5.append("python %s/run/scripts/getResults.py -j %s -m %s -v %s -i %s -l %s -n \'%s\' -t %s  -a %s -b %s -c %s -d %s -e %s -f %s -g %s -k %s " %(TargetOrtho_path,jobID,n+1,param_dic['max_offset_variance'],output_dir,speciesStr,matrixNamesList[n],param_dic["top_ranked_per_gene"],param_dic['scale_A'],param_dic['scale_B'],param_dic['scale_C'],param_dic['scale_D'],param_dic['scale_E'],param_dic['scale_F'],param_dic['scale_G'],TargetOrtho_path))
        return queue1,queue2,queue3,queue4,queue5
    
    #assign a jobID, make output dir,load parameters
    jobID,jobTag=assign_jobID()    
    output_dir="%s/run/output/%s" %(TargetOrtho_path,jobID)
    setup.mkMainDirs(TargetOrtho_path)
    param_dic,matrix_count,matrixNamesList=load_cmd_line_params()
    print matrixNamesList        
    print_parameters(param_dic),'final params'

    try:
        setup.mkOutputDir(jobID,output_dir,matrixNamesList,param_dic['motif_infile'],param_dic['QueryFile'])
        setup.writeInputSum(jobID,jobTag,param_dic,matrixNamesList,output_dir)
        
        #build queues for parallel execution of scripts
        queue1,queue2,queue3,queue4,queue5=build_queues(matrix_count,speciesList,param_dic,matrixNamesList,jobID,gene_name)
        
        #execute fimo for each matrix for each species using parallel threads
        threadScripts.execute_queue(queue1) 

        #check fimo output
        for species in speciesList:
            for PSSMnum in range(matrix_count):
                filesize=float(os.path.getsize('%s/run/fimo_out/%s_%s_%s_fimo.out' %(TargetOrtho_path,jobID,species,PSSMnum+1)))
                if filesize < 86:
                    message="FIMO p value thresshold too stringent. No motif matches found for %s,motif %s" %(species,matrixNamesList[PSSMnum])
                    raise Exception(message)

        #execute associate_genes.py for each matrix for each species using parallel theads
        threadScripts.execute_queue(queue2)
        
        #execute findOrthoAssoc.py for each matrix for each species using parallel threads
        threadScripts.execute_queue(queue3)    
        
        #execute conserved_ortho.py
        threadScripts.execute_queue(queue4)
        
        #execute getResults1.py
        threadScripts.execute_queue(queue5)
        
        #execute cis_regulatory_hits.py if matrix_count > 1
        if matrix_count > 1:
            os.system("python %s/run/scripts/cis_reg_hits.py -j %s -p %s -b %s -r %s  -m %s -l %s -c %s" %(TargetOrtho_path,jobID,param_dic['region_overlap'],param_dic['bp_dist'],param_dic['motif_order'],matrix_count,speciesStr,TargetOrtho_path))
        
        #execute checkQueryList.py (if queryFile specified)     
        if (param_dic["QueryFile"]!=None) and (param_dic["QueryFile"] !='None'):
            os.system("python %s/run/scripts/checkQueryList.py -j %s -q %s -m %s -r %s -v %s -s %s -c %s" %(TargetOrtho_path,jobID,param_dic["QueryFile"],matrix_count, param_dic["QueryOnly"],param_dic['%s_genomeVersion' %(speciesList[0])],speciesList[0],TargetOrtho_path))
        
        #execute mkHTMLtables.py
        matrixNames=''
        for n in matrixNamesList:
            print n,'n'
            print n.split()[0]
            matrixNames=matrixNames + n.split()[0] + '*'
        os.system("python %s/run/scripts/mkHTMLtables.py -j %s -g \'%s\' -i %s -c %s  -m %s -n \'%s\' -q %s -w %s -S %s -c %s" %(TargetOrtho_path,jobID,jobTag,output_dir,TargetOrtho_path,matrix_count,matrixNames,param_dic["QueryFile"],param_dic["QueryOnly"],speciesStr,TargetOrtho_path))
    
        print jobID,'jobID'        
        #clear all job data (except output directory)
        os.system("python %s/run/scripts/clear_job_data.py -j %s -c %s" %(TargetOrtho_path,jobID,TargetOrtho_path))            
        print 'results at: %s' %(output_dir)        
        
    except:
        #clear all job data (except output directory)
        print "Unexpected error:", sys.exc_info()[0]        
        os.system("python %s/run/scripts/clear_job_data.py -j %s -c %s" %(TargetOrtho_path,jobID,TargetOrtho_path))    
        raise    
    
main()

