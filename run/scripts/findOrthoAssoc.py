#!/usr/bin/python
#Author: Lori Glenwinkel
import sys
from optparse import OptionParser
import sqlite3 as lite

parser = OptionParser()
parser.add_option( '-n', '--name', dest='name', help='name (eleg,japo,brig,bren,rema)' )
parser.add_option( '-m', '--matrix', dest='matrix', help='martrix (-m 1 -m 2 -m 3)' )
parser.add_option('-v', '--version', dest='version',help='genome version: -v WS220')
parser.add_option('-j','--jobID',dest='jobID',help='-j 123345')
parser.add_option('-g','--gene_name',dest='gene_name',help='-g unc-3')
parser.add_option('-r', '--QueryOnly', dest='QueryOnly',help='-r True')
parser.add_option('-f', '--QueryFilePath', dest='QueryFilePath',help='-f QueryFilePath')
parser.add_option( '-Z', '--max_dist_if_genes_between' , dest='max_dist_if_genes_between', help='absolute nucleotide distance to search if genes are present between site and associated gene -z 2000')
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/newTargetOrtho/run' )
parser.add_option('-s','--speciesList',help='speciesList: -s [c_eleg,c_brig,c_bren,c_japo,c_rema]')
( options, args ) = parser.parse_args()
TargetOrtho_path=options.TargetOrtho_path


gene_name=options.gene_name
max_dist_if_genes_between=options.max_dist_if_genes_between
jobID=options.jobID
species = options.name 
PSSMnum = options.matrix
version=options.version
speciesStr=options.speciesList.split('-')
speciesList=[n for n in speciesStr]

conn = lite.connect('%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,species,PSSMnum))
cursor = conn.cursor()
print "%s_%s%s_TargetOrtho.db" %(jobID,species,PSSMnum)
def getQueryListIDs(QueryFilePath):
    infile=file('%s' %(QueryFilePath),'r')
    lines=infile.readlines()
    gene_nameList=[]
    for line in lines:
        l=line.split()
        if len(l)>0:
            if '#' in l[0]:
                queryName=l[0][1:]
            else:
                queryName='queryList'
                gene_nameList.append(l[0])
    return gene_nameList,queryName

def prepare_sql_tables(jobID,species,PSSMnum):
    try:
        cursor.execute("drop table if exists %s_%s_hit_ortho%s" %(jobID,species,PSSMnum))
        cursor.execute("create table %s_%s_hit_ortho%s (ortho_id int(10),association_id int(10), site_id int(10),species varchar(30),site_sequence varchar(100),dna varchar(40),strand varchar(10),gene_name varchar(100),offset int(15),genes_between int(10),ortho_in_%s varchar(100),score double(10,4),p_value double(10,4),start int(10),end int(10),site_region varchar(200))" %(jobID,species,PSSMnum,speciesList[0]))
        cursor.execute("drop table if exists %s_%s_temp%s" %(jobID,species,PSSMnum))
        cursor.execute("create table %s_%s_temp%s as select * from  %s_%s_associate_genes%s where 0" %(jobID,species,PSSMnum,jobID,species,PSSMnum))
    except lite.Error, e:
        if conn:
            conn.rollback()
            print "Error %s:" % e.args[0]
            sys.exit(1)
    
def associate_gene(gene_names):
    cursor.execute("attach database '%s/run/TargetOrtho.db' as c" %(TargetOrtho_path))
    exonID_list=[]
    l=float(len(set(gene_names)))
    assocCount=0
    for gene_name in set(gene_names):
        assocCount+=1
        n=(int(PSSMnum)-1)*5
        sys.stdout.write("%s%.2f \r" % ('\t'*(n+speciesList.index(species)),assocCount/l*100) )
        sys.stdout.flush()        
        if species!=speciesList[0]:
            ortho_ids=[row[0] for row in cursor.execute("select distinct ortholog_gid from %s_%s_ortho where %s_gid='%s'" %(speciesList[0],species,speciesList[0],gene_name))]            
            for ortho_id in ortho_ids:                
                exonIDs=list([row[0] for row in cursor.execute("select exonID from %s_exons%s where gene_name='%s'" %(species,version,ortho_id.strip()))])
                exonID_list=exonID_list + exonIDs
        else:
            exonIDs=list([row[0] for row in cursor.execute("select exonID from %s_exons%s where gene_name=\'%s\'" %(species,version,gene_name))])
            exonID_list=exonID_list+exonIDs
    assocCount=0
    cursor.execute("select * from %s_%s_associate_genes%s" %(jobID,species,PSSMnum))
    assocs=cursor.fetchall()
    l=float(len(assocs))
    for assoc in assocs:
        assocCount+=1
        #print assoc,'assoc from associate_genes'
        #(0, 286, u'HEAD', u'SAME', 1050, 0, 46288) assoc from associate_genes
        #0|association_id|int(10)|0||0
        #1|site_id|int(10)|0||0
        #2|site_region|varchar(30)|0||0
        #3|sameStrand|varchar(10)|0||0
        #4|assoc_distance|int(10)|0||0
        #5|GeneBetween|int(10)|0||0
        #6|exonID|int(10)|0||0
        #sys.exit(1)
        if int(assoc[6]) in exonID_list:
            cursor.execute("insert into %s_%s_temp%s (association_id,site_id,site_region,GeneBetween,sameStrand,assoc_distance,exonID) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')"%(jobID,species,PSSMnum,assoc[0],assoc[1],assoc[2],assoc[5],assoc[3],assoc[4],assoc[6]))
        n=(int(PSSMnum)-1)*5
        sys.stdout.write("%s%.2f \r" % ('\t'*(n+speciesList.index(species)),assocCount/l*100) )
        sys.stdout.flush()
    conn.commit()
    
def associate_ortho(QueryOnly):
    try:cursor.execute("attach database '%s/run/TargetOrtho.db' as c" %(TargetOrtho_path))#change back to TargetOrtho.db
    except:None
    if QueryOnly=='True':cursor.execute("select * from %s_%s_temp%s" %(jobID,species,PSSMnum))
    else:cursor.execute("select * from %s_%s_associate_genes%s" %(jobID,species,PSSMnum))
    assocs=cursor.fetchall()
    ortho_id=0
    assocCount=0
    l=float(len(assocs))
    conn.commit()
    for assoc in assocs:
        assocCount+=1
        n=(int(PSSMnum)-1)*5
        sys.stdout.write("%s%.2f \r" % ('\t'*(n+speciesList.index(species)),assocCount/l*100) )
        sys.stdout.flush()
        #print assoc,'assoc',assoc[6],'assoc[6]'
        #print "select gene_name from %s_exons%s where exonID = %s limit 1" %(species,version,assoc[6]),'match fetch'
        
        cursor.execute("select gene_name from %s_exons%s where exonID = %s limit 1" %(species,version,assoc[6]))
        match=cursor.fetchone()
        #print match,'match'
        #print "select site_sequence,score,start,end,dna,strand from %s_%s_hits%s where id = %s limit 1" %(jobID,species,PSSMnum,int(assoc[1]))
        cursor.execute("select site_sequence,score,p_value,start,end,dna,strand from %s_%s_hits%s where id = %s limit 1" %(jobID,species,PSSMnum,int(assoc[1])))
        hit=cursor.fetchone()                
        cursor.execute("select strand,start,end from %s_gene_info%s where gene_name= '%s' limit 1" %(species,version,match[0].strip()))
        gene_info=cursor.fetchone()   

        #print hit,'hit'
        #print gene_info,'gene_info'
        #print match,'match'
        #print assoc,'assoc'
     
        strand=int(gene_info[0])*-1        
        site_region=assoc[2]            
        x1=int(gene_info[0])*-1
        offset=x1*(int(gene_info[1])-int(hit[3]))       
        #print 'offset = gene_info[1]-int(hit[2]))'
        #print strand,'strand'
        #print site_region,'site_region'
        #print offset,'offset'
        
        if site_region== 'HEAD' or site_region=='TAIL':
            if offset > 0:
                site_region='downstream'
                offsetD=x1*(int(gene_info[2])-int(hit[3])) 
                offset=offsetD#uncomment to allow separate downstream offset measured from stop codon
            elif offset < 0 or offset==0:site_region='upstream'
        if species!=speciesList[0]:
            command="select distinct %s_gid from %s_%s_ortho where ortholog_gid= '%s'" %(speciesList[0],speciesList[0],species,match[0])
            cursor.execute(command)
            eleg_ids=cursor.fetchall()
            conn.commit()
            for eleg_id in eleg_ids:  
                ortho_id+=1
                command="insert into %s_%s_hit_ortho%s (ortho_id,association_id, site_id,species,site_region,ortho_in_%s, offset,genes_between,gene_name,site_sequence,score,p_value,start,end,dna,strand) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,species,PSSMnum,speciesList[0],ortho_id,assoc[0],assoc[1],species, site_region, eleg_id[0], offset, assoc[5], match[0], hit[0],hit[1],hit[2],hit[3],hit[4],hit[5],hit[6])
                #print assoc,'assoc'
                #print hit,'hit'
                #print command, 'command'
                #print

                cursor.execute(command)
                conn.commit()
        else:
            if 'X' in hit[5]:
                dna=hit[5][hit[5].index('X'):]
            elif 'IV' in hit[5]:
                dna=hit[5][hit[5].index('IV'):]
            elif 'V' in hit[5]:
                dna=hit[5][hit[5].index('V'):]
            elif 'III' in hit[5]:
                dna=hit[5][hit[5].index('III'):]
            elif 'II' in hit[5]:
                dna=hit[5][hit[5].index('II'):]
            elif 'I' in hit[5]:
                dna=hit[5][hit[5].index('I'):]                        
            else:
                dna=hit[5]                    
            ortho_id+=1
            ortho_in_ref=match[0]
            command="insert into %s_%s_hit_ortho%s (ortho_id,association_id, site_id,species,site_region,ortho_in_%s, offset,genes_between,gene_name,site_sequence,score,p_value,start,end,dna,strand) values(\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,species,PSSMnum,speciesList[0],ortho_id,assoc[0],assoc[1],species,site_region,ortho_in_ref,offset,assoc[5],match[0],hit[0],hit[1],hit[2],hit[3],hit[4],dna,hit[6])
            #print assoc,'assoc'
            #print hit,'hit'
            #print command, 'command'
            #print

            cursor.execute(command)
    conn.commit()
        
def main():
    print 'starting associated_ortho for ', species
    prepare_sql_tables(jobID,species,PSSMnum)   
    conn.commit()
    if options.QueryOnly=='True':
        gene_names,queryName=getQueryListIDs(options.QueryFilePath)
        conn.commit()
        associate_gene(gene_names)
    associate_ortho(options.QueryOnly)
    cursor.execute("delete from %s_%s_hit_ortho%s where genes_between > 0 and abs(offset) > %s " %(jobID,species,PSSMnum,max_dist_if_genes_between))
    conn.close()
main()
