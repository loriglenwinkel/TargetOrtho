#!/usr/bin/python
#Author: Lori Glenwinkel
# CREATE DATE: 2014
# DESCRIPTION: filters TargetOrtho results for user defined genes of interest (from query list input)

import os,re
import sys
from optparse import OptionParser
import sqlite3 as lite

parser = OptionParser()
parser.add_option('-q','--file',dest='filepath',help='-f queryFiles.txt')
parser.add_option('-m','--matrix_count',dest='matrix_count',help='-m 3')
parser.add_option('-j','--jobID',dest='jobID',help='-j 123345')
parser.add_option('-r', '--QueryOnly', dest='QueryOnly',help='-r True')
parser.add_option('-v','--version',dest='version',help='-v WS233')
parser.add_option('-s','--ref_species',dest='ref_species',help='-s c_eleg or -s d_mel')
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/newTargetOrtho/run' )
( options, args ) = parser.parse_args()
TargetOrtho_path=options.TargetOrtho_path
ref_species=options.ref_species
jobID=options.jobID
matrix_count=int(options.matrix_count)
filepath=options.filepath
filename=filepath.split('/')[-1]
QueryOnly=options.QueryOnly
version=options.version
conn = lite.connect('%s/run/sqlite_tmp_files/%s_TargetOrtho.db' %(TargetOrtho_path,jobID),isolation_level=None)
cursor = conn.cursor()

for m in range(matrix_count):
    cursor.execute("attach database '%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' as %s%s_TargetOrtho" %(TargetOrtho_path,jobID,m+1,jobID,m+1))
cursor.execute("attach database '%s/run/TargetOrtho.db' as c" %(TargetOrtho_path))
def prepare_sql_tables(jobID,matrix_count):
    #cursor=conn.cursor()
    print matrix_count
    for n in range(matrix_count):
        cursor.execute("drop table if exists %s_QueryListResults%s" %(jobID,n+1))
        cursor.execute("create table %s_QueryListResults%s as select * from  %s_All_conserved_hits%s_ranked where 0" %(jobID,n+1,jobID,matrix_count))
        if n==1:
            print 'creating QueryListResults_CRM table'
            cursor.execute("drop table if exists %s_QueryListResults_CRM" %(jobID))
            cursor.execute("create table %s_QueryListResults_CRM as select * from %s_CRM_results where 0" %(jobID,jobID))

def transfer_result(gene_name,PSSMnum):
    #cursor=conn.cursor()
    print 'transferring to query table ', PSSMnum
    cursor.execute("insert into %s_QueryListResults%s select * from %s_All_conserved_hits%s_ranked where associated_gene1='%s'" %(jobID,PSSMnum,jobID,PSSMnum,gene_name[0]))

def transfer_result_CRM(gene_name):
    print gene_name[0],'gene_name[gid]'
    #cursor=conn.cursor()
    cursor.execute("insert into %s_QueryListResults_CRM select * from %s_CRM_results where associated_gene1='%s'" %(jobID,jobID,gene_name[0]))

def lookup_query(name_dic,PSSMnum):
    cursor.execute("select distinct associated_gene1 from %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum))
    gids=cursor.fetchall()
    foundList=[]
    for gid in gids:
        name=gid[0]
        if name in name_dic.keys():
            transfer_result(gid,PSSMnum)
            foundList.append(name)
            print name, 'found'
        elif name in name_dic.values():
            transfer_result(gid,PSSMnum)
            foundList.append(name)
            print name,'found'
    return foundList
def lookup_query_CRM(name_dic):
    cursor.execute("select distinct associated_gene1 from %s_CRM_results" %(jobID))
    gids=cursor.fetchall()
    foundList=[]
    for gid in gids:
        name=gid[0]
        if name in name_dic.keys():
            transfer_result_CRM(gid)
            foundList.append(name)
        elif name in name_dic.values():
            transfer_result_CRM(gid)
            foundList.append(name)
    return foundList

def getSeqID(gene_nameList):
    #cursor=conn.cursor()
    name_dic={}
    for gene_name in gene_nameList:
        cursor.execute("select seq_id from %s_gene_info%s where gene_name = '%s' limit 0,1" %(ref_species,version,gene_name))
        seq_id=cursor.fetchone()
        if seq_id !=None:name_dic[gene_name]=seq_id[0]
    return name_dic

def readQueryList(filepath):
    infile=file('%s' %(filepath),'r')
    lines=infile.readlines()
    gene_nameList=[]
    for line in lines:
        print line
        l=line.split()
        if '#' in l[0]:
            queryName=l[0][1:]
        else:gene_nameList.append(l[0])
    return gene_nameList,queryName

def writeQueryResults(foundList,total,PWMnum,Query_name):
    #cursor=conn.cursor()
    foundString=''
    for n in foundList:
        foundString = foundString + n + ''
    if len(foundList)<1:percentFound=0
    else:percentFound=float(len(foundList))/total*100
    cursor.execute("insert into %s_ResultsSummary%s (result,count) values('percent of Query Genes with hits from %s (out of %s total)','%s')" %(jobID,PWMnum,Query_name,total,percentFound))
    cursor.execute("insert into %s_ResultsSummary%s (result,count) values('genes with Hits from %s','%s')" %(jobID,PWMnum,Query_name,foundString))


def main():
    print 'running checkQueryList.py'
    prepare_sql_tables(jobID,matrix_count)
    for n in range(matrix_count):
        PSSMnum=n+1
        if QueryOnly=='True':
            print 'inserting top ranked hits results into querylist results'
            cursor.execute("insert into %s_QueryListResults%s select * from %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum,jobID,PSSMnum))
            if PSSMnum==2:cursor.execute("insert into %s_QueryListResults_CRM select * from %s_CRM_results" %(jobID,jobID))
        else:
            print ''
            gene_nameList,queryName=readQueryList(filepath)
            name_dic=getSeqID(gene_nameList)
            foundList=lookup_query(name_dic,PSSMnum)
            total=len(gene_nameList)
            writeQueryResults(foundList,total,PSSMnum,queryName)
            if PSSMnum==2:
                foundList=lookup_query_CRM(name_dic)
    conn.close()
main()




