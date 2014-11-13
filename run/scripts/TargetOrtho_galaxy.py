import os,sys,datetime
#run TargetOrtho tool via Galaxy with this script. Also requires TargetOrtho.xml
p=os.path.abspath(__file__).split('/')[1:-3]
s=''
for i in p:
    s=s +'/%s'%i
TargetOrtho_path=s

def assign_jobID():
    dt=datetime.datetime.now()
    jobID="g%s%s%s%s%s%s" %(dt.year,dt.month,dt.day,dt.hour,dt.minute,dt.second)
    return jobID

def get_galaxy_input():
    dic={
        'jobTag':sys.argv[1],
        'PSSMinfile':sys.argv[2],
        'QueryFile':sys.argv[3],
        'QueryOnly':sys.argv[5],
        'reference_species':sys.argv[6]}
    dic["jobID"]=assign_jobID()
    if dic['jobTag']=='_':dic['jobTag']=dic['jobID']
    if len(sys.argv)>8:
        assert (sys.argv[8]==('None') or (int(sys.argv[8])==0 or int(sys.argv[8])>0)),'check job name, no spaces allowed. upstream distance must be None, 0, or positive integer'
        dic['max_upstream']=sys.argv[8]
        assert (sys.argv[9]==('None') or (int(sys.argv[9])==0 or int(sys.argv[9])>0)),'downstream distance must be None, 0, or positive integer'
        dic['max_downstream']=sys.argv[9]        
        dic['filter_exons']=sys.argv[10]
        dic['region_overlap']=sys.argv[11]    
        assert int(sys.argv[12])>-1,'genes between must be postive integer or 0'
        dic['genes_between']=int(sys.argv[12])        
        assert int(sys.argv[13])>-1,'distance allowed if genes between must be postive integer or 0'
        dic['dist_if_genes_between']=int(sys.argv[13])
        assert int(sys.argv[14])>0, 'top ranked per gene must positive integer'
        dic['top_ranked_per_gene']=int(sys.argv[14])        
        assert float(sys.argv[15])>0, "max offset variance must be greater than 0"
        dic['max_offset_variance']=float(sys.argv[15])
        dic['scale_A']=float(sys.argv[16])
        dic['scale_B']=float(sys.argv[17])
        dic['scale_C']=float(sys.argv[18])
        dic['scale_D']=float(sys.argv[19])
        dic['scale_E']=float(sys.argv[20])
        dic['scale_F']=float(sys.argv[21])
        dic['scale_G']=float(sys.argv[22])                
        for n in [16,17,18,19,20,21,22]:
            assert float(sys.argv[n])>=0,"species weight must be greater than or equal to 0"
        dic['scale_D']=float(sys.argv[19])
        dic['scale_E']=float(sys.argv[20])
        dic['scale_F']=float(sys.argv[21])
        dic['scale_G']=float(sys.argv[22])                
        for n in [16,17,18,19,20,21,22]:
            assert float(sys.argv[n])>=0,"species weight must be greater than or equal to 0"
        if int(sys.argv[23])>0:            
            assert(float(sys.argv[24])>0 and float(sys.argv[24])<=1),'P value must be greater than 0 and less than or equal to 1. (1e-4 or 0.0001)'
            dic['p_thresh']='%s' %sys.argv[24]
        if int(sys.argv[23])>1:     
            dic['bp_dist']=int(sys.argv[25])
            dic['motif_order']=sys.argv[26]
            assert(float(sys.argv[27])>0 and float(sys.argv[27])<=1),'P value must be greater than 0 and less than or equal to 1. (1e-4 or 0.0001)'
            dic['p_thresh']='%s %s' %(sys.argv[24], sys.argv[27])
        if int(sys.argv[23])>2:
            assert(float(sys.argv[28])>0 and float(sys.argv[28])<=1),'P value must be greater than 0 and less than or equal to 1. (1e-4 or 0.0001)'
            dic['p_thresh']='%s %s %s' %(sys.argv[24], sys.argv[27], sys.argv[28])
        if int(sys.argv[23])>3:
            assert(float(sys.argv[29])>0 and float(sys.argv[29])<=1),'P value must be greater than 0 and less than or equal to 1. (1e-4 or 0.0001)'
            dic['p_thresh']='%s %s %s %s' %(sys.argv[24], sys.argv[27], sys.argv[28], sys.argv[29])
        if int(sys.argv[23])>4:
            assert(float(sys.argv[30])>0 and float(sys.argv[30])<=1),'P value must be greater than 0 and less than or equal to 1. (1e-4 or 0.0001)'
            dic['p_thresh']='%s %s %s %s %s' %(sys.argv[24], sys.argv[27], sys.argv[28], sys.argv[29], sys.argv[30])
        bp_dist=''
        motif_order=''
        if sys.argv[23]>'1':
            bp_dist='-d %s'%dic["bp_dist"]
            motif_order='-y %s'%dic["motif_order"]
        command="python %s/run/scripts/TargetOrtho.py -g %s -J \'%s\' -O %s -f %s -q %s -x %s -i %s -v %s -e %s -P %s -k %s %s %s -a %s -b %s  -H %s  -I %s  -K %s  -L %s  -N  %s -Z %s -z %s -w %s -p \'%s\' >%s/run/%sTargetOrthoStdout.out 2>&1" %(TargetOrtho_path,dic["jobID"],dic["jobTag"],dic['reference_species'],dic["PSSMinfile"],dic["QueryFile"],dic["max_upstream"],dic["max_downstream"],dic["max_offset_variance"],dic["filter_exons"],dic["region_overlap"],dic["top_ranked_per_gene"],bp_dist,motif_order,dic["scale_A"],dic["scale_B"],dic["scale_C"],dic["scale_D"],dic["scale_E"],dic["scale_F"],dic["scale_G"],dic["dist_if_genes_between"],dic["genes_between"],dic["QueryOnly"],dic['p_thresh'],TargetOrtho_path,dic["jobID"])
        os.system(command)                    
    else:
        command="python %s/run/scripts/TargetOrtho.py -g %s -J \'%s\' -O %s -f %s -q %s -w %s >%s/run/%sTargetOrthoStdout.out 2>&1" %(TargetOrtho_path,dic["jobID"],dic["jobTag"],dic['reference_species'],dic["PSSMinfile"],dic["QueryFile"],dic["QueryOnly"],TargetOrtho_path,dic["jobID"])
        os.system(command)        
    return dic

def main():
    dic=get_galaxy_input()
    output_dir = "%s/run/output/%s" %(TargetOrtho_path,dic["jobID"])    
    #print the stdout logg, in case errors have occurs (force galaxy to display stdout from tool)
    stdout_file="%s/run/%sTargetOrthoStdout.out" %(TargetOrtho_path,dic["jobID"])
    for l in file(stdout_file,'r').readlines():
        print l
    #move the stdout logg to job output directory
    os.system("mv %s/run/%sTargetOrthoStdout.out %s/" %(TargetOrtho_path,dic["jobID"],output_dir))
    #move job output directory to galaxy output directory 
    os.system("mv  %s %s/" %(output_dir,sys.argv[7]))
    #copy index.html with main results to galaxy output
    os.system("cp %s/index.html %s" %(sys.argv[7],sys.argv[4]))

main()
