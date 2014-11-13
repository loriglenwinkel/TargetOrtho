#!/usr/bin/python
import itertools,sys
import sqlite3 as lite
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-m','--PSSMnum',help='PSSMnum: --PSSMnum 1')
parser.add_option('-s','--speciesList',help='speciesList: -s [c_eleg,c_brig,c_bren,c_japo,c_rema]')
parser.add_option( '-p', '--region_overlap', dest='region_overlap', help='-p True #True if motif genomic region should be shared in all orthologs i.e. all motif matches in exon, intron, intergenic-downstream, or intergenic-upstream')
parser.add_option( '-x', '--max_upstream', dest='max_upstream', help='-x None or -x 1000')
parser.add_option( '-i', '--max_downstream', dest='max_downstream', help='-i None or -x 1000')
parser.add_option( '-k', '--top_ranked_per_gene', dest='top_ranked_per_gene', type="int",help='-k 1')
parser.add_option('-j','--jobID',dest='jobID',help='-j 2012420')
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/newTargetOrtho' )
( options, args ) = parser.parse_args()
TargetOrtho_path=options.TargetOrtho_path
top_ranked_per_gene=options.top_ranked_per_gene
global PSSMnum
jobID=options.jobID
PSSMnum=int(options.PSSMnum)
region_overlap=options.region_overlap
max_upstream=options.max_upstream
max_downstream=options.max_downstream

if max_upstream !='None':max_upstream="offset > %s" %(-1*int(max_upstream))
else:max_upstream ='offset > -1000000'
if max_downstream !='None':max_downstream="and offset < %s" %(max_downstream)
else:max_downstream='and offset < 1000000'

speciesStr=options.speciesList.split('-')
speciesList=[n for n in speciesStr]
print TargetOrtho_path
conn = lite.connect('%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,PSSMnum),isolation_level=None)
cursor = conn.cursor()
for species in speciesList:
    cursor.execute("attach database '%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' as %s_%s%s_db" %(TargetOrtho_path,jobID,species,PSSMnum,jobID,species,PSSMnum))

def prepare_sql_tables():
    for i in range(0,5):        
        def addCols(table,cons):
            try:
                j=['varchar(100)','varchar(100)','int(10)','varchar(30)','double(10,4)','int(10)','varchar(40)','varchar(10)','int(10)','varchar(200)']
                factors=['species','gene_name','site_id','site_sequence','score','start','dna','strand','offset','site_region']
                for k in factors:
                    cursor.execute("alter table %s_%s_conserved_hits%s add column %s%s %s" %(jobID,table,PSSMnum,k,cons,j[factors.index(k)]))
            except lite.Error, e:
                if conn:
                    conn.rollback()
                    print "Error %s:" % e.args[0]
                    sys.exit(1)
        try:
            cursor.execute("drop table if exists %s_%s_conserved_hits%s" %(jobID,i+1,PSSMnum))
            cursor.execute("create table %s_%s_conserved_hits%s (species1 varchar(100),gene_name1 varchar(100),site_id1 int(10),site_sequence1 varchar(30),score1 int(10),start1 int(10), dna1 varchar(40),strand1 varchar(10),offset1 int(10),site_region1 varchar(200))" %(jobID,i+1,PSSMnum))
            print i+1,PSSMnum,'species number, motifnum'
        except lite.Error, e:
                if conn:
                    conn.rollback()
                    print "Error %s:" % e.args[0]
                    sys.exit(1)
    for i in range(1,5):
        addCols(i+1,2)
    for i in range(2,5):
        addCols(i+1,3)
    for i in range(3,5):
        addCols(i+1,4)
    addCols(5,5)

def fill_conserved_hit_table(combo,foundList):
    if conservation ==1:cursor.execute("insert into %s_1_conserved_hits%s select species,gene_name, site_id,site_sequence, score,start,dna,strand,offset,site_region from %s_%s_hit_ortho%s where ortho_id = \'%s\'" %(jobID,PSSMnum,jobID,speciesList[0],PSSMnum,combo[0]))
    
    if conservation ==2:cursor.execute("insert into %s_%s_conserved_hits%s select eleg.species as species1,eleg.gene_name as gene_name1,eleg.site_id as site_id1, eleg.site_sequence as site_sequence1, eleg.score as score1, eleg.start as start1, eleg.dna as dna1, eleg.strand as strand1, eleg.offset as offset1,eleg.site_region as site_region1, spc2.species as species2,spc2.gene_name as gene_name2,spc2.site_id as site_id2, spc2.site_sequence as site_sequence2, spc2.score as score2, spc2.start as start2, spc2.dna as dna2, spc2.strand as strand2, spc2.offset as offset2,spc2.site_region as site_region2 from %s_%s_hit_ortho%s as eleg, %s_%s_hit_ortho%s as spc2  where eleg.ortho_id = \'%s\' and spc2.ortho_id = \'%s\' " %(jobID,conservation,PSSMnum,jobID,speciesList[0],PSSMnum,jobID,foundList[1],PSSMnum,combo[0],combo[1]))

    if conservation ==3:cursor.execute("insert into %s_%s_conserved_hits%s select eleg.species as species1,eleg.gene_name as gene_name1,eleg.site_id as site_id1, eleg.site_sequence as site_sequence1, eleg.score as score1, eleg.start as start1, eleg.dna as dna1, eleg.strand as strand1, eleg.offset as offset1,eleg.site_region as site_region1, spc2.species as species2,spc2.gene_name as gene_name2,spc2.site_id as site_id2, spc2.site_sequence as site_sequence2, spc2.score as score2, spc2.start as start2, spc2.dna as dna2, spc2.strand as strand2, spc2.offset as offset2,spc2.site_region as site_region2, spc3.species as species3,spc3.gene_name as gene_name3,spc3.site_id as site_id3, spc3.site_sequence as site_sequence3, spc3.score as score3, spc3.start as start3, spc3.dna as dna3, spc3.strand as strand3, spc3.offset as offset3,spc3.site_region as site_region3 from %s_%s_hit_ortho%s as eleg, %s_%s_hit_ortho%s as spc2, %s_%s_hit_ortho%s as spc3 where eleg.ortho_id = \'%s\' and spc2.ortho_id = \'%s\' and spc3.ortho_id = \'%s\' " %(jobID,conservation,PSSMnum,jobID,speciesList[0],PSSMnum,jobID,foundList[1],PSSMnum,jobID,foundList[2],PSSMnum,combo[0],combo[1],combo[2]))
    if conservation ==4:cursor.execute("insert into %s_%s_conserved_hits%s select eleg.species as species1,eleg.gene_name as gene_name1,eleg.site_id as site_id1, eleg.site_sequence as site_sequence1, eleg.score as score1, eleg.start as start1, eleg.dna as dna1, eleg.strand as strand1, eleg.offset as offset1,eleg.site_region as site_region1, spc2.species as species2,spc2.gene_name as gene_name2,spc2.site_id as site_id2, spc2.site_sequence as site_sequence2, spc2.score as score2, spc2.start as start2, spc2.dna as dna2, spc2.strand as strand2, spc2.offset as offset2,spc2.site_region as site_region2, spc3.species as species3,spc3.gene_name as gene_name3,spc3.site_id as site_id3, spc3.site_sequence as site_sequence3, spc3.score as score3, spc3.start as start3, spc3.dna as dna3, spc3.strand as strand3, spc3.offset as offset3,spc3.site_region as site_region3, spc4.species as species4,spc4.gene_name as gene_name4,spc4.site_id as site_id4, spc4.site_sequence as site_sequence4, spc4.score as score4, spc4.start as start4, spc4.dna as dna4, spc4.strand as strand4, spc4.offset as offset4,spc4.site_region as site_region4 from %s_%s_hit_ortho%s as eleg, %s_%s_hit_ortho%s as spc2, %s_%s_hit_ortho%s as spc3, %s_%s_hit_ortho%s as spc4  where eleg.ortho_id = \'%s\' and spc2.ortho_id = \'%s\' and spc3.ortho_id = \'%s\' and spc4.ortho_id = \'%s\' "%(jobID,conservation,PSSMnum,jobID,speciesList[0],PSSMnum,jobID,foundList[1],PSSMnum,jobID,foundList[2],PSSMnum,jobID,foundList[3],PSSMnum,combo[0],combo[1],combo[2],combo[3]))
    if conservation ==5:cursor.execute("insert into %s_%s_conserved_hits%s select eleg.species as species1,eleg.gene_name as gene_name1,eleg.site_id as site_id1, eleg.site_sequence as site_sequence1, eleg.score as score1, eleg.start as start1, eleg.dna as dna1, eleg.strand as strand1, eleg.offset as offset1,eleg.site_region as site_region1,spc2.species as species2,spc2.gene_name as gene_name2,spc2.site_id as site_id2, spc2.site_sequence as site_sequence2, spc2.score as score2, spc2.start as start2, spc2.dna as dna2, spc2.strand as strand2, spc2.offset as offset2,spc2.site_region as site_region2, spc3.species as species3,spc3.gene_name as gene_name3,spc3.site_id as site_id3, spc3.site_sequence as site_sequence3, spc3.score as score3, spc3.start as start3, spc3.dna as dna3, spc3.strand as strand3, spc3.offset as offset3,spc3.site_region as site_region3, spc4.species as species4,spc4.gene_name as gene_name4,spc4.site_id as site_id4, spc4.site_sequence as site_sequence4, spc4.score as score4, spc4.start as start4, spc4.dna as dna4, spc4.strand as strand4, spc4.offset as offset4,spc4.site_region as site_region4, spc5.species as species5,spc5.gene_name as gene_name5,spc5.site_id as site_id5, spc5.site_sequence as site_sequence5, spc5.score as score5, spc5.start as start5, spc5.dna as dna5, spc5.strand as strand5, spc5.offset as offset5,spc5.site_region as site_region5 from %s_%s_hit_ortho%s as eleg, %s_%s_hit_ortho%s as spc2, %s_%s_hit_ortho%s as spc3, %s_%s_hit_ortho%s as spc4,%s_%s_hit_ortho%s as spc5 where eleg.ortho_id = \'%s\' and spc2.ortho_id = \'%s\' and spc3.ortho_id = \'%s\' and spc4.ortho_id = \'%s\' and spc5.ortho_id = \'%s\'"   %(jobID,conservation,PSSMnum,jobID,speciesList[0],PSSMnum,jobID,foundList[1],PSSMnum,jobID,foundList[2],PSSMnum,jobID,foundList[3],PSSMnum,jobID,foundList[4],PSSMnum,combo[0],combo[1],combo[2],combo[3],combo[4]))    
def get_all_gene_combos(all_species_ids_List,foundList):
    all_ids=[]
    for n in all_species_ids_List: 
        for i in n: all_ids.append(i)
    all_combos= set(itertools.combinations(all_ids,conservation))
    for combo in all_combos:
        if combo[0] in (all_species_ids_List[0]):
            if conservation == 1: fill_conserved_hit_table(combo,foundList)
            elif conservation > 1:            
                if combo[1] in species2_ids:
                    if conservation ==2: fill_conserved_hit_table(combo,foundList)
                    elif conservation > 2:
                        if combo[2] in species3_ids:
                            if conservation ==3: fill_conserved_hit_table(combo,foundList)
                            elif conservation > 3:
                                if combo[3] in species4_ids:
                                    if conservation ==4: fill_conserved_hit_table(combo,foundList)
                                    elif conservation > 4:
                                        if combo[4] in species5_ids: fill_conserved_hit_table(combo,foundList)                                        
                                    
def get_eleg_gid_list(eleg_gid,species,region_overlap,eleg_region,offset):
    if region_overlap=='True':region_overlap="and site_region=\'%s\' " %(eleg_region)
    else:region_overlap = ''
    cursor.execute("select ortho_id from %s_%s_hit_ortho%s where ortho_in_%s = \'%s\' %s and %s %s order by abs(offset - %s) limit 1" %(jobID,species,PSSMnum,speciesList[0],eleg_gid,region_overlap,max_upstream,max_downstream,offset))
    ortho_id=cursor.fetchone()
    if ortho_id==None:
        cursor.execute("select ortho_id from %s_%s_hit_ortho%s where ortho_in_%s = \'%s\' %s and %s and site_region='INTRON' order by abs(offset - %s) limit 1" %(jobID,species,PSSMnum,speciesList[0],eleg_gid,region_overlap,max_upstream,offset))
        ortho_id=cursor.fetchone()
    if ortho_id!=None:
        return [ortho_id[0]]
    else: return []

def progress(Rcount,l):
    if int(PSSMnum)==1:sys.stdout.write("%sprogress_%s: %s   \r" % ('',PSSMnum,Rcount/l*100) )
    else:sys.stdout.write("%sprogress_%s: %s   \r" % ('\t\t\t'*int(PSSMnum-1),PSSMnum,Rcount/l*100) )        
    sys.stdout.flush()            
def main():
    print 'running conserved_ortho.py'
    prepare_sql_tables()
    Rcount=0            
    command= "select offset,site_region as eleg_region,ortho_id,gene_name,site_id from %s_%s_hit_ortho%s where %s %s order by score DESC" %(jobID,speciesList[0],PSSMnum,max_upstream,max_downstream)
    command2= "select offset,site_region as eleg_region,ortho_id,gene_name,site_id from %s_%s_hit_ortho%s where %s and site_region='INTRON' order by score DESC" %(jobID,speciesList[0],PSSMnum,max_upstream)
    cursor.execute(command)
    eleg_ortho_rows=cursor.fetchall()            
    cursor.execute(command2)
    eleg_ortho_rows=eleg_ortho_rows+cursor.fetchall()    
    l=float(len(eleg_ortho_rows))
    for row in eleg_ortho_rows:
        Rcount+=1   
        foundList=['%s' %speciesList[0]]                
        progress(Rcount,l)        
        global eleg_ortho_id
        eleg_ortho_id=row[2]
        eleg_gid=row[3]
        species1=speciesList[0]
        species1_ids=[eleg_ortho_id]
        global species2_ids,species3_ids,species4_ids,species5_ids
        offset=int(row[0])
        for species in speciesList[1:]:                        
            species_ids=get_eleg_gid_list(eleg_gid,species,region_overlap,row[1],row[0])                        
            if len(species_ids) > 0:
                foundList.append(species)
                if len(foundList)==2:species2,species2_ids=species,species_ids
                elif len(foundList)==3:species3,species3_ids=species,species_ids
                elif len(foundList)==4:species4,species4_ids=species,species_ids
                elif len(foundList)==5:species5,species5_ids=species,species_ids
        global conservation
        conservation = len(foundList)
        if conservation == 5: all_species_ids_list=[species1_ids,species2_ids,species3_ids,species4_ids,species5_ids]
        elif conservation == 4:all_species_ids_list=[species1_ids,species2_ids,species3_ids,species4_ids]
        elif conservation == 3:all_species_ids_list=[species1_ids,species2_ids,species3_ids]
        elif conservation == 2:all_species_ids_list=[species1_ids,species2_ids]
        elif conservation ==1: all_species_ids_list=[species1_ids]
        get_all_gene_combos(all_species_ids_list,foundList)
    conn.close()
main()

