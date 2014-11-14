#!/usr/bin/python
#Author: Lori Glenwinkel
# CREATE DATE: 2014
# DESCRIPTION: finds co-occurances of motifs by cross referencing results in All_conserved_hits tables from each input motif

import sys
from optparse import OptionParser
import sqlite3 as lite

parser = OptionParser()
parser.add_option( '-p', '--feature_overlap', dest='feature_overlap', help='-p True #True if motif genomic feature should be shared in all orthologs. i.e. all motif matches in exon, intron, intergenic-downstream, or intergenic-upstream')
parser.add_option( '-b', '--distance_between', dest='bp_dist', type="int",help='-b 1000  or -b None #maximum distance in nucleotides between motif matches (optional when searching for combination of 2 or more motifs matches. If None chosen, then any distance between motifs is allowed ')
parser.add_option( '-r', '--ordered', dest='PSSMorder', help='True or False: -o True means motif matches must be same order as input file on DNA')
parser.add_option('-m','--matrix_count',help='matrix_count: -m 3')
parser.add_option('-j','--jobID',dest='jobID',help='-j 123345')
parser.add_option('-l', '--speciesStr', dest='speciesStr',help='-l \'eleg brig bren rema japo\'')
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/newTargetOrtho/run' )

( options, args ) = parser.parse_args()
TargetOrtho_path=options.TargetOrtho_path
speciesStr=options.speciesStr.split('-')
speciesList=[n for n in speciesStr]
jobID=options.jobID
feature_overlap=options.feature_overlap
bp_dist=int(options.bp_dist)
#bp_dist=50#remove
PSSMorder=options.PSSMorder
matrix_count=int(options.matrix_count)

conn = lite.connect('%s/run/sqlite_tmp_files/%s_TargetOrtho.db' %(TargetOrtho_path,jobID),isolation_level=None)
cursor = conn.cursor()
for m in range(matrix_count):
    cursor.execute("attach database '%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' as %s%s_TargetOrtho" %(TargetOrtho_path,jobID,m+1,jobID,m+1))

def create_CRM_hits_table(matrix_count):
    print matrix_count,'matrix_count'
    print 'creating CRM table'
    if matrix_count <2:sys.exit(1)
    cursor.execute("drop table if exists  %s_CRM_results" %(jobID)) 
    if matrix_count==2:
        cursor.execute("create table %s_CRM_results (associated_gene1 varchar(100),motif1_site_id1 int(10),motif2_site_id1 int(10),motif1_site_rank int(10),motif2_site_rank int(10),motif1_conservation1 int(10),motif2_conservation1 int(10),motif1_site_seq1 varchar(30), motif2_site_seq1 varchar(30),motif1_site_position1 varchar(50),motif2_site_position1 varchar(50),motif1_site_offset1 int(10),motif2_site_offset1 int(10),motif1_site_region1 varchar(30),motif2_site_region1 varchar(30),motif1_site_count1 int(10),motif2_site_count1 int(10),motif1_offset_variance double(10,7),motif2_offset_variance double(10,7))" %(jobID))
    if matrix_count==3:
        cursor.execute("create table %s_CRM_results (associated_gene1 varchar(100),motif1_site_id1 int(10),motif2_site_id1 int(10),motif3_site_id1 int(10),motif1_site_rank int(10),motif2_site_rank int(10),motif3_site_rank int(10),motif1_conservation1 int(10),motif2_conservation1 int(10),motif3_conservation1 int(10),motif1_site_seq1 varchar(30), motif2_site_seq1 varchar(30),motif3_site_seq1 varchar(30),motif1_site_position1 varchar(50),motif2_site_position1 varchar(50),motif3_site_position1 varchar(50),motif1_site_offset1 int(10),motif2_site_offset1 int(10),motif3_site_offset1 int(10),motif1_site_region1 varchar(30),motif2_site_region1 varchar(30),motif3_site_region1 varchar(30),motif1_site_count1 int(10),motif2_site_count1 int(10),motif3_site_count1 int(10),motif1_offset_variance double(10,7),motif2_offset_variance double(10,7),motif3_offset_variance double(10,7))" %(jobID))
    if matrix_count==4:
        cursor.execute("create table %s_CRM_results (associated_gene1 varchar(100),motif1_site_id1 int(10),motif2_site_id1 int(10),motif3_site_id1 int(10),motif4_site_id1,motif1_site_rank int(10),motif2_site_rank int(10),motif3_site_rank int(10),motif4_site_rank int(10),motif1_conservation1 int(10),motif2_conservation1 int(10),motif3_conservation1 int(10),motif4_conservation1 int(10),motif1_site_seq1 varchar(30), motif2_site_seq1 varchar(30),motif3_site_seq1 varchar(30),motif4_site_seq1 varchar(30),motif1_site_position1 varchar(50),motif2_site_position1 varchar(50),motif3_site_position1 varchar(50),motif4_site_position1 varchar(50),motif1_site_offset1 int(10),motif2_site_offset1 int(10),motif3_site_offset1 int(10),motif4_site_offset1 int(10),motif1_site_region1 varchar(30),motif2_site_region1 varchar(30),motif3_site_region1 varchar(30),motif4_site_region1,motif1_site_count1 int(10),motif2_site_count1 int(10),motif3_site_count1 int(10),motif4_site_count1 int(10),motif1_offset_variance double(10,7),motif2_offset_variance double(10,7),motif3_offset_variance double(10,7),motif4_offset_variance double(10,7))" %(jobID))
    if matrix_count==5:
        cursor.execute("create table %s_CRM_results (associated_gene1 varchar(100),motif1_site_id1 int(10),motif2_site_id1 int(10),motif3_site_id1 int(10),motif4_site_id1,motif5_site_id1,motif1_site_rank int(10),motif2_site_rank int(10),motif3_site_rank int(10),motif4_site_rank int(10),motif5_site_rank int(10),motif1_conservation1 int(10),motif2_conservation1 int(10),motif3_conservation1 int(10),motif4_conservation1 int(10),motif5_conservation1 int(10),motif1_site_seq1 varchar(30), motif2_site_seq1 varchar(30),motif3_site_seq1 varchar(30),motif4_site_seq1 varchar(30),motif5_site_seq1 varchar(30),motif1_site_position1 varchar(50),motif2_site_position1 varchar(50),motif3_site_position1 varchar(50),motif4_site_position1 varchar(50),motif5_site_position1 varchar(50),motif1_site_offset1 int(10),motif2_site_offset1 int(10),motif3_site_offset1 int(10),motif4_site_offset1 int(10),motif5_site_offset1 int(10),motif1_site_region1 varchar(30),motif2_site_region1 varchar(30),motif3_site_region1 varchar(30),motif4_site_region1,motif5_site_region1 varchar(30),motif1_site_count1 int(10),motif2_site_count1 int(10),motif3_site_count1 int(10),motif4_site_count1 int(10),motif5_site_count1 int(10),motif1_offset_variance double(10,7),motif2_offset_variance double(10,7),motif3_offset_variance double(10,7),motif4_offset_variance double(10,7),motif5_offset_variance double(10,7))" %(jobID))

    cursor.execute("drop table if exists  %s_top_ranked_per_gene_CRM" %(jobID))
    print "creating top ranked per gene CRM table"
    cursor.execute("create table %s_top_ranked_per_gene_CRM as select * from  %s_CRM_results where 0" %(jobID,jobID))    
    cursor.execute("drop table if exists  %s_temp" %(jobID))
    cursor.execute("create table %s_temp as select * from  %s_CRM_results where 0" %(jobID,jobID))

def get_All_cons_ids(matrix_count,conn):
    cursor=conn.cursor()
    common_ids=[]
    for n in range(matrix_count):
        cursor.execute("select distinct associated_gene1 from %s_All_conserved_hits%s_ranked" %(jobID,n+1))
        if n==0:hits1=cursor.fetchall()
        elif n==1:
            hits2=cursor.fetchall()
            common_ids=(set(hits1) & set(hits2))
        elif n==2:
            hits3=cursor.fetchall()
            common_ids=(set(hits1) & set(hits2) & set(hits3))
        elif n==3:
            hits4=cursor.fetchall()
            common_ids=(set(hits1) & set(hits2) & set(hits3) & set(hits4))
        elif n==4:
            hits5=cursor.fetchall()
            common_ids=(set(hits1) & set(hits2) & set(hits3) & set(hits4) & set(hits5))            
    return common_ids

def check_dist(offsetsList,PSSMorder):
    if bp_dist=='None':return True
    if PSSMorder=='True':
        if matrix_count>1:
            dist1=abs(offsetsList[0]-offsetsList[1])
            if dist1>bp_dist:return False
            elif matrix_count==2:return True
            if matrix_count>2:
                dist2=abs(offsetsList[1]-offsetsList[2])
                if dist2>bp_dist:return False
                elif matrix_count==3:return True
                if matrix_count>3:
                    dist3=abs(offsetsList[2]-offsetsList[3])
                    if dist3>bp_dist:return False
                    elif matrix_count==4:return True
                    if matrix_count>4:
                        dist4=abs(offsetsList[3]-offsetsList[4])
                        if dist4>bp_dist:return False
                        elif matrix_count==4:return True
    else:                  
        max_pos=max(offsetsList)
        min_pos=min(offsetsList)
        mid_pos=(max_pos+min_pos)/2
        check_offsets=True
        for i in offsetsList:
            if abs(i-mid_pos)>bp_dist:check_offsets=False
        return check_offsets

def check_overlap(featuresList,feature_overlap):
    if feature_overlap=='False':return True
    else:
        if len(set(featuresList))==1:return True
        else:return False
        
def fill_CRM_hits(gids):
    if matrix_count==1:
        print 'CRM for unconserved combos coming soon'
    elif matrix_count==2:
        cursor.execute("insert into %s_CRM_results (associated_gene1,motif1_conservation1,motif2_conservation1,motif1_site_id1,motif2_site_id1,motif1_site_seq1,motif2_site_seq1,motif1_site_position1,motif2_site_position1,motif1_site_offset1,motif2_site_offset1,motif1_site_region1,motif2_site_region1,motif1_site_count1,motif2_site_count1,motif1_site_rank,motif2_site_rank,motif1_offset_variance,motif2_offset_variance) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,gids[0][0],gids[0][3],gids[1][3],gids[0][4],gids[1][4],gids[0][5],gids[1][5],gids[0][6],gids[1][6],gids[0][1],gids[1][1],gids[0][2],gids[1][2],gids[0][7],gids[1][7],gids[0][8],gids[1][8],gids[0][9],gids[1][9]))
    elif matrix_count==3:
        cursor.execute("insert into %s_CRM_results (associated_gene1,motif1_conservation1,motif2_conservation1,motif3_conservation1,motif1_site_id1,motif2_site_id1,motif3_site_id1,motif1_site_seq1,motif2_site_seq1,motif3_site_seq1,motif1_site_position1,motif2_site_position1,motif3_site_position1,motif1_site_offset1,motif2_site_offset1,motif3_site_offset1,motif1_site_region1,motif2_site_region1,motif3_site_region1,motif1_site_count1,motif2_site_count1,motif3_site_count1,motif1_site_rank,motif2_site_rank,motif3_site_rank,motif1_offset_variance,motif2_offset_variance,motif3_offset_variance) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,gids[0][0],gids[0][3],gids[1][3],gids[2][3],gids[0][4],gids[1][4],gids[2][4],gids[0][5],gids[1][5],gids[2][5],gids[0][6],gids[1][6],gids[2][6],gids[0][1],gids[1][1],gids[2][1],gids[0][2],gids[1][2],gids[2][2],gids[0][7],gids[1][7],gids[2][7],gids[0][8],gids[1][8],gids[2][8],gids[0][9],gids[1][9],gids[2][9]))
    elif matrix_count==4:
        cursor.execute("insert into %s_CRM_results (associated_gene1,motif1_conservation1,motif2_conservation1,motif3_conservation1,motif4_conservation1,motif1_site_id1,motif2_site_id1,motif3_site_id1,motif4_site_id1,motif1_site_seq1,motif2_site_seq1,motif3_site_seq1,motif4_site_seq1,motif1_site_position1,motif2_site_position1,motif3_site_position1,motif4_site_position1,motif1_site_offset1,motif2_site_offset1,motif3_site_offset1,motif4_site_offset1,motif1_site_region1,motif2_site_region1,motif3_site_region1,motif4_site_region1,motif1_site_count1,motif2_site_count1,motif3_site_count1,motif4_site_count1,motif1_site_rank,motif2_site_rank,motif3_site_rank,motif4_site_rank,motif1_offset_variance,motif2_offset_variance,motif3_offset_variance,motif4_offset_variance) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,gids[0][0],gids[0][3],gids[1][3],gids[2][3],gids[3][3],gids[0][4],gids[1][4],gids[2][4],gids[3][4],gids[0][5],gids[1][5],gids[2][5],gids[3][5],gids[0][6],gids[1][6],gids[2][6],gids[3][6],gids[0][1],gids[1][1],gids[2][1],gids[3][1],gids[0][2],gids[1][2],gids[2][2],gids[3][2],gids[0][7],gids[1][7],gids[2][7],gids[3][7],gids[0][8],gids[1][8],gids[2][8],gids[3][8],gids[0][9],gids[1][9],gids[2][9],gids[3][9]))
    elif matrix_count==5:
        cursor.execute("insert into %s_CRM_results (associated_gene1,motif1_conservation1,motif2_conservation1,motif3_conservation1,motif4_conservation1,motif5_conservation1,motif1_site_id1,motif2_site_id1,motif3_site_id1,motif4_site_id1,motif5_site_id1,motif1_site_seq1,motif2_site_seq1,motif3_site_seq1,motif4_site_seq1,motif5_site_seq1,motif1_site_position1,motif2_site_position1,motif3_site_position1,motif4_site_position1,motif5_site_position1,motif1_site_offset1,motif2_site_offset1,motif3_site_offset1,motif4_site_offset1,motif5_site_offset1,motif1_site_region1,motif2_site_region1,motif3_site_region1,motif4_site_region1,motif5_site_region1,motif1_site_count1,motif2_site_count1,motif3_site_count1,motif4_site_count1,motif5_site_count1,motif1_site_rank,motif2_site_rank,motif3_site_rank,motif4_site_rank,motif5_site_rank,motif1_offset_variance,motif2_offset_variance,motif3_offset_variance,motif4_offset_variance,motif5_offset_variance) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,gids[0][0],gids[0][3],gids[1][3],gids[2][3],gids[3][3],gids[4][3],gids[0][4],gids[1][4],gids[2][4],gids[3][4],gids[4][4],gids[0][5],gids[1][5],gids[2][5],gids[3][5],gids[4][5],gids[0][6],gids[1][6],gids[2][6],gids[3][6],gids[4][6],gids[0][1],gids[1][1],gids[2][1],gids[3][1],gids[4][1],gids[0][2],gids[1][2],gids[2][2],gids[3][2],gids[4][2],gids[0][7],gids[1][7],gids[2][7],gids[3][7],gids[4][7],gids[0][8],gids[1][8],gids[2][8],gids[3][8],gids[4][8],gids[0][9],gids[1][9],gids[2][9],gids[3][9],gids[4][9]))

def find_CRM_hits(common_ids,conn):
    for id in common_ids:
        cursor.execute("select associated_gene1,site_offset1,site_region1,conservation,site_id1,site_sequence1,site_position1,region_site_count1,rank,offset_variance from %s_All_conserved_hits1_ranked where associated_gene1=\'%s\'" %(jobID,id[0]))
        gids1=cursor.fetchall()
        for gid1 in gids1:
            cursor.execute("select associated_gene1,site_offset1,site_region1,conservation,site_id1,site_sequence1,site_position1,region_site_count1,rank,offset_variance from %s_All_conserved_hits2_ranked where associated_gene1=\'%s\'" %(jobID,gid1[0]))
            gids2=cursor.fetchall()
            for gid2 in gids2:
                if matrix_count==2:
                    offsetsList=[int(gid1[1]),int(gid2[1])]
                    featuresList=[gid1[2],gid2[2]]
                    if check_dist(offsetsList,PSSMorder)==True:
                        if check_overlap(featuresList,feature_overlap)==True:
                            fill_CRM_hits([gid1,gid2])
                elif matrix_count > 2:
                    cursor.execute("select associated_gene1,site_offset1,site_region1,conservation,site_id1,site_sequence1,site_position1,region_site_count1,rank,offset_variance from %s_All_conserved_hits3_ranked where associated_gene1=\'%s\'" %(jobID,gid2[0]))
                    gids3=cursor.fetchall()
                    for gid3 in gids3:
                        if matrix_count==3:
                            offsetsList=[int(gid1[1]),int(gid2[1]),int(gid3[1])]
                            featuresList=[gid1[2],gid2[2],gid3[2]]
                            if check_dist(offsetsList,PSSMorder)==True:
                                if check_overlap(featuresList,feature_overlap)==True:
                                    fill_CRM_hits([gid1,gid2,gid3])
                                
                        elif matrix_count > 3:
                            cursor.execute("select associated_gene1,site_offset1,site_region1,conservation,site_id1,site_sequence1,site_position1,region_site_count1,rank,offset_variance  from %s_All_conserved_hits4_ranked where associated_gene1=\'%s\'" %(jobID,gid3[0]))
                            gids4=cursor.fetchall()
                            for gid4 in gids4:
                                if matrix_count==4:
                                    offsetsList=[int(gid1[1]),int(gid2[1]),int(gid3[1]),int(gid4[1])]
                                    featuresList=[gid1[2],gid2[2],gid3[2],gid4[2]]
                                    if check_dist(offsetsList,PSSMorder)==True:
                                        if check_overlap(featuresList,feature_overlap)==True:
                                            fill_CRM_hits([gid1,gid2,gid3,gid4])

                                elif matrix_count > 4:
                                    cursor.execute("select associated_gene1,site_offset1,site_region1,conservation,site_id1,site_sequence1,site_position1,region_site_count1,rank,offset_variance  from %s_All_conserved_hits5_ranked where associated_gene1=\'%s\'" %(jobID,gid4[0]))
                                    gids5=cursor.fetchall()
                                    for gid5 in gids5:
                                        if matrix_count==5:
                                            offsetsList=[int(gid1[1]),int(gid2[1]),int(gid3[1]),int(gid4[1]),int(gid5[1])]
                                            featuresList=[gid1[2],gid2[2],gid3[2],gid4[2]]
                                            if check_dist(offsetsList,PSSMorder)==True:
                                                if check_overlap(featuresList,feature_overlap)==True:
                                                    fill_CRM_hits([gid1,gid2,gid3,gid4,gid5])


def remove_dups():
    cursor.execute("drop table if exists  %s_CRM" %(jobID))
    cursor.execute("create table %s_CRM as select * from %s_CRM_results where 0" %(jobID,jobID))
    cursor.execute("insert into %s_CRM select distinct* from %s_CRM_results" %(jobID,jobID))
    cursor.execute("delete from %s_CRM_results" %(jobID))
    cursor.execute("insert into %s_CRM_results select * from %s_CRM" %(jobID,jobID))

def getTopRanked(matrix_count):
    print 'getting top ranked per gene'
    cursor.execute("select distinct associated_gene1 from %s_CRM_results" %(jobID))
    eleg_gids=cursor.fetchall()            
    rank3=''
    rank4=''
    rank5=''
    if matrix_count >2:rank3=' +motif3_site_rank'
    if matrix_count >3:rank4=' +motif4_site_rank'
    if matrix_count >4:rank5=' +motif5_site_rank'
    for eleg_gid in eleg_gids:         
        
        print eleg_gid[0]
        cursor.execute("insert into %s_temp  select * from %s_CRM_results where associated_gene1=\'%s\' limit 0,1" %(jobID,jobID,eleg_gid[0]))        
    cursor.execute("insert into %s_top_ranked_per_gene_CRM select * from %s_temp order by (motif1_site_rank + motif2_site_rank%s%s%s) ASC" %(jobID,jobID,rank3,rank4,rank5))        

def main():    
    create_CRM_hits_table(matrix_count)
    common_ids=get_All_cons_ids(matrix_count,conn)
    find_CRM_hits(common_ids,conn)
    remove_dups()
    getTopRanked(matrix_count)
    conn.close()
main()
