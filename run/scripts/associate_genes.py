#!/usr/bin/python
import sys
import os
from optparse import OptionParser
import sqlite3 as lite
import fileinput

parser = OptionParser()
parser.add_option( '-n', '--name', dest='name', help='name (eleg,japo,brig,bren,rema): -n eleg' )
parser.add_option( '-m', '--matrix', dest='matrix', help='matrix (PSSM number): -m 1 -m 2 -m 3)' )
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/TargetOrtho/run' )
parser.add_option( '-v', '--genome_version', dest='version', help=': genome release version: -v WS220' )
parser.add_option( '-z', '--genes_between', dest='genes_between', help='number of genes allowed between motif match and associated gene: -z 1' )
parser.add_option( '-Z', '--max_dist_if_genes_between', dest='max_dist_if_genes_between', help='maximum absolute distance  allowed between motif match and associated gene if there are intervening genes: -Z 2000' )
parser.add_option( '-x', '--filter_exons', dest='filter_exons', help='set to True to remove exons from analysis: -x True (otherwise do not include option' )
parser.add_option('-j','--jobID',dest='jobID',help='-j 2012420')
parser.add_option('-p','--max_upstream',dest='max_upstream',help='-u 2000')
parser.add_option('-w','--max_downstream',dest='max_downstream',help='=w 500')
( options, args ) = parser.parse_args()

max_upstream=options.max_upstream
max_downstream=options.max_downstream
jobID=options.jobID
TargetOrtho_path =options.TargetOrtho_path
version=options.version
genes_between=options.genes_between
max_dist_if_genes_between=options.max_dist_if_genes_between
filter_exons=options.filter_exons

species = options.name 

if filter_exons=='True':filter_exons='-x True'
else: filter_exons=''

PSSMnum = options.matrix

conn = lite.connect('%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,species,PSSMnum))
#conn.isolation_level = None
cursor = conn.cursor()

def prepare_sql_tables(jobID,species,TargetOrtho_path):

    try:
        cursor.execute("drop table if exists %s_%s_hits%s" %(jobID,species,PSSMnum))
        cursor.execute("create table %s_%s_hits%s (id int(10),dna varchar(40),start int(10),end int(10),strand varchar(10),score double(10,4),p_value double(10,4), site_sequence varchar(100))" %(jobID,species,PSSMnum))    
        conn.commit()
    except lite.Error, e:
        if conn:
            conn.rollback()
            print "Error %s:" % e.args[0]
            sys.exit(1)    
    
    try:
        cursor.execute("drop table if exists  %s_%s_associate_genes%s" %(jobID,species,PSSMnum))
        cursor.execute("create table %s_%s_associate_genes%s (association_id int(10),site_id int(10),site_region varchar(30),sameStrand varchar(10),assoc_distance int(10),GeneBetween int(10),exonID int(10))" %(jobID,species,PSSMnum))        
        conn.commit()
    except lite.Error, e:
        if conn:
            conn.rollback()
            print "Error %s:" % e.args[0]
            sys.exit(1)
    
def mv_hits_output_to_sql():
    print 'writing hits file for associate genes'    
    print '%s/run/fimo_out/%s_%s_hits%s.txt' %(TargetOrtho_path,jobID,species,PSSMnum),'outfile'
    outfile=file('%s/run/fimo_out/%s_%s_hits%s.txt' %(TargetOrtho_path,jobID,species,PSSMnum),'w')
    filesize=float(os.path.getsize('%s/run/fimo_out/%s_%s_%s_fimo.out' %(TargetOrtho_path,jobID,species,PSSMnum)))
    print filesize
    if filesize < 86:
        message="FIMO p value thresshold too stringent. No motif matches found for %s" %species
        raise Exception(message)
 
    if filesize > 50000000:#20000000:
        print 'using slow method to preserve memory'
        hit_id=0
        line_count=0
        for line in fileinput.input('%s/run/fimo_out/%s_%s_%s_fimo.out' %(TargetOrtho_path,jobID,species,PSSMnum)):
            l=line.split('\t')
            print line, 'fimo line'
            if line_count>0:
                hit_id+=1
                if l[4].strip()=='+':strand='POS'
                else:strand='NEG'
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(hit_id,species,l[1],l[2],l[3],strand))
                j=[hit_id]+['%s' %j.strip() for j in l[1:]]
                k=j[:-2]+j[-1:]
                cursor.execute("INSERT INTO %s_%s_hits%s (id,dna,start,end,strand,score,p_value,site_sequence) VALUES (?, ?, ?, ?, ?, ?, ?, ?);" %(jobID,species,PSSMnum), k)
                conn.commit()
            line_count+=1

    else:
        print 'using fast method'
        data=[row.split('\t') for row in fileinput.input('%s/run/fimo_out/%s_%s_%s_fimo.out' %(TargetOrtho_path,jobID,species,PSSMnum),bufsize=1000)]
        data2=[]
        hit_id=0
        for i in data[1:]:
            hit_id+=1
            if i[4].strip()=='+':strand='POS'
            else:strand='NEG'
            outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' %(hit_id,species,i[1],i[2],i[3],strand))
            j=[hit_id]+['%s' %j.strip() for j in i[1:]]
            k=j[:-2]+j[-1:]
            data2.append(k)
        #upload hits table to sqlite db
        cursor.executemany("INSERT INTO %s_%s_hits%s (id,dna,start,end,strand,score,p_value,site_sequence) VALUES (?, ?, ?, ?, ?, ?, ?, ?);" %(jobID,species,PSSMnum), data2)
        conn.commit()
    
def associate_genes():

    #execute assocate_genes script using fimo output (hits file) as input
    print 'running associate genes command'
    command="%s/run/scripts/associate_genes -b 10000000 -m %s -o %s/run/associate_genes_out/%s_%s_associate_genes%s -d %s/genomes/%s_%s/%s.dnas -p %s/genomes/%s_%s/ -q %s/run/fimo_out/%s_%s_hits%s.txt -e %s/run/exon_files/%s_exons%s.txt %s" %(TargetOrtho_path,genes_between,TargetOrtho_path,jobID,species,PSSMnum,TargetOrtho_path,species,version,species,TargetOrtho_path,species,version,TargetOrtho_path,jobID,species,PSSMnum,TargetOrtho_path,species,version,filter_exons)
    print command,'associate_genes command'
    os.system(command)

def upload_assoc_genes():    
    #determine greatest distance restrction for filtering associate_genes results.
    if max_upstream =='None':p =1000000
    else:p=max_upstream
    if max_downstream =='None':w=1000000
    else:w=max_downstream
    m=max([int(p),int(w)])
    print m, 'max'
    cursor.execute('PRAGMA synchronous = OFF')
    conn.isolation_level = None
    print "removing hits beyond max association distance"
    
    print 'upload hit assocation to mysql'
    def upload(x):
        #print x
        try:cursor.execute("INSERT INTO %s_%s_associate_genes%s (association_id,site_id,site_region,sameStrand,assoc_distance,GeneBetween,exonID) VALUES (?, ?, ?, ?, ?, ?, ?);" %(jobID,species,PSSMnum), x)
        except lite.Error, e:                
            if conn:
                conn.rollback()
                print "Error %s:" % e.args[0]
                sys.exit(1)

    #data=[d.split('\t') for d in fileinput.input('%s/run/associate_genes_out/%s_%s_associate_genes%s' %(TargetOrtho_path,jobID,species,PSSMnum))]
    for d in fileinput.input('%s/run/associate_genes_out/%s_%s_associate_genes%s' %(TargetOrtho_path,jobID,species,PSSMnum)):
        i=d.split('\t')
        if int(i[4].strip())<m:
            if int(i[5].strip())<int(genes_between):
                if int(i[5].strip())>0:
                    if abs(int(i[4].strip()))< int(max_dist_if_genes_between):upload(['%s' %j.strip() for j in i])
                else:upload(['%s' %j.strip() for j in i])

def main():        
    prepare_sql_tables(jobID,species,TargetOrtho_path)
    mv_hits_output_to_sql()
    associate_genes()    
    upload_assoc_genes()    
    conn.close()
main()              


