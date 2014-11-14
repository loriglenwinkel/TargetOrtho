#!/usr/bin/python
#Author: Lori Glenwinkel
# CREATE DATE: 2014
# DESCRIPTION: creates All_conserved_hits tables from 1-5 conserved hits tables (made by conserved_ortho.py). Updates binding site features, calculates cumulative site score, rank orders all binding sites. Creats top_ranked_hits tables.

import sys
import numpy as np
import  os
import matplotlib
matplotlib.use('Agg')#disables interactive mode
import matplotlib.pyplot as plt
import sys,scipy.stats
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-m','--PSSMnum',dest='PSSMnum',help='--PSSMnum 3')
parser.add_option( '-v', '--max_offset_variance', dest='max_offset_variance', help='-v 1')
parser.add_option('-i', '--output_dir', dest='output_dir',help='-i /data/newTargetOrtho/run/output/123232')
parser.add_option('-l', '--speciesStr', dest='speciesStr',help='-l \'eleg brig bren rema japo\'')
parser.add_option('-n', '--PSSMname',dest='PSSMname',help='-f unc-3')
parser.add_option('-j','--jobID',dest='jobID',help='-j 2012420')
parser.add_option('-t','--top_ranked',dest='top_ranked_per_gene',help='-t 1 ..number of toped ranked sites to show per gene in top ranked hits table')
parser.add_option('-a','--scale_A',dest='scale_A',help='-a 1 ..factor to scale the avg species site score before ranking final results')
parser.add_option('-b','--scale_B',dest='scale_B',help='-b 1 ..factor to scale the avg species region score before ranking final results')
parser.add_option('-c','--scale_C',dest='scale_C',help='-c 1 ..factor to scale the avg species region site count before ranking final results')
parser.add_option('-d','--scale_D',dest='scale_D',help='-d 1 ..factor to scale the avg species total gene score before ranking final results')
parser.add_option('-e','--scale_E',dest='scale_E',help='-e 1 ..factor to scale the avg species gene site counts  before ranking final results')
parser.add_option('-f','--scale_F',dest='scale_F',help='-f 1 ..factor to scale the conservation score before ranking final results')
parser.add_option('-g','--scale_G',dest='scale_G',help='-g 1 ..factor to scale the avg species offset variance  before ranking final results')
parser.add_option( '-k', '--TargetOrtho_path', dest='TargetOrtho_path', help='-k /data/newTargetOrtho' )
(options, args) = parser.parse_args()

TargetOrtho_path=options.TargetOrtho_path
scale_A=options.scale_A
scale_B=options.scale_B
scale_C=options.scale_C
scale_D=options.scale_D
scale_E=options.scale_E
scale_F=options.scale_F
scale_G=options.scale_G
print scale_G,'scale_G'
jobID=options.jobID
PSSMnum=int(options.PSSMnum)
max_offset_variance=options.max_offset_variance
output_dir=options.output_dir
PSSMname=options.PSSMname
top_ranked_per_gene=options.top_ranked_per_gene
speciesStr=options.speciesStr.split('-')
speciesList=[]
for n in speciesStr[:-1]:
    speciesList.append('%s' %(n))

import sqlite3 as lite
conn = lite.connect('%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,PSSMnum),isolation_level=None)
cursor = conn.cursor()
for species in speciesList:
    cursor.execute("attach database '%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' as %s_%s%s_db" %(TargetOrtho_path,jobID,species,PSSMnum,jobID,species,PSSMnum))

def prepare_sql_tables():
    def addSpecies(cons):
        if cons==1:return "species%s varchar(100),rank int(10),cum_score double(10,4),associated_gene%s varchar(100),site_id%s int(10),site_position%s varchar(100),site_strand%s varchar(100),site_sequence%s varchar(100),site_region%s varchar(100),site_offset%s int(10), conservation int(10),site_score%s double(10,4),avg_region_score%s double(10,4),avg_gene_score%s double(10,4),region_site_count%s double(10,4),gene_site_count%s double(10,4)"%(cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons)
        else:return "species%s varchar(100),associated_gene%s varchar(100),site_id%s int(10),site_position%s varchar(100),site_strand%s varchar(100),site_sequence%s varchar(100),site_region%s varchar(100),site_offset%s int(10),site_score%s double(10,4),avg_region_score%s double(10,4),avg_gene_score%s double(10,4),region_site_count%s double(10,4),gene_site_count%s double(10,4)"%(cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons,cons)
    cursor.execute("drop table if exists %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum))    
    cursor.execute("create table %s_All_conserved_hits%s_ranked (ID INTEGER PRIMARY KEY,%s,%s,%s,%s,%s, avg_species_site_score varchar(100),avg_species_region_score varchar(100),avg_species_gene_score varchar(100),avg_species_region_site_count varchar(100),avg_species_gene_site_count varchar(100),offset_variance varchar(100))" %(jobID,PSSMnum,addSpecies(1),addSpecies(2),addSpecies(3),addSpecies(4),addSpecies(5)))                            
    #print "create table %s_All_conserved_hits%s_ranked (ID INTEGER PRIMARY KEY,%s,%s,%s,%s,%s, avg_species_site_score varchar(100),avg_species_region_score varchar(100),avg_species_gene_score varchar(100),avg_species_region_site_count varchar(100),avg_species_gene_site_count varchar(100),offset_variance varchar(100))" %(jobID,PSSMnum,addSpecies(1),addSpecies(2),addSpecies(3),addSpecies(4),addSpecies(5))
    cursor.execute("drop table if exists %s_top_ranked_hits%s" %(jobID,PSSMnum))
    cursor.execute("create table  %s_top_ranked_hits%s as select * from %s_All_conserved_hits%s_ranked where 0" %(jobID,PSSMnum,jobID,PSSMnum))
    cursor.execute("drop table if exists %s_inputSummary" %(jobID))
    cursor.execute("create table %s_inputSummary(parameter text,usr_input text)" %(jobID))
    cursor.execute("drop table if exists %s_ResultsSummary%s" %(jobID,PSSMnum))
    cursor.execute("create table %s_ResultsSummary%s(result text,count text)" %(jobID,PSSMnum))

def autolabel(rects):
    for rect in rects:
        height = rect.get_height()
        try:plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),ha='center', va='bottom',fontsize='8')
        except:print 'not using autolabel, error occured'

def mkPlot1_2(PSSMnum,uniqueHitsList,uniqueGenesList,SpeciesCountList):
    plt.figure(1+(PSSMnum-1)*5)
    ind=np.arange(5)
    width = .35
    plt.subplot(111)
    rects1 = plt.bar(ind, uniqueHitsList, width,color='b')
    rects2 = plt.bar(ind+width, uniqueGenesList, width,color='orange')
    plt.ylabel('frequency')
    plt.title('%s UniqueGeneCounts' %(PSSMname))
    plt.xticks(ind + .5*width, ('c_eleg ','2 genomes','3 genomes','4 genomes','5 genomes'))
    plt.legend((rects1[0],rects2[0]), ('All sites','Unique Genes'))
    autolabel(rects1)
    autolabel(rects2)
    plt.savefig("%s/%s/Results_plots/1_UniqueGenesCounts_%s.png" %(output_dir,PSSMname.split()[0],PSSMname),format = 'png')    
    plt.figure(2+(PSSMnum-1)*5)
    plt.subplot(111)
    SpeciesCountPlot=plt.bar(ind,SpeciesCountList, width, color = 'g')
    plt.ylabel('percent species coverage')
    plt.title('Species representation among conserved Hits for %s' %(PSSMname))
    width=.35
    plt.xticks(ind + .5*width, speciesList)
    autolabel(SpeciesCountPlot)
    plt.savefig("%s/%s/Results_plots/2_SpeciesRepresentation%s.png" %(output_dir,PSSMname.split()[0],PSSMname),format = 'png')
    
    def writeResults1():
        cursor=conn.cursor()
        cursor.execute("delete from %s_ResultsSummary%s" %(jobID,PSSMnum))
        xlabelList=['reference species ','2 species','3 species','4 species','5 species']
        n=5
        for i in range(len(xlabelList)):
            n-=1
            if n!=0:cursor.execute("insert into %s_ResultsSummary%s (result,count)  values( 'total motif matches conserved in %s, (number of genes with at least one conserved motif match across %s)','%s, (%s)')" %(jobID,PSSMnum,xlabelList[n],xlabelList[n],uniqueHitsList[n],uniqueGenesList[n]))
            else:cursor.execute("insert into %s_ResultsSummary%s (result,count)  values( 'total motif matches found in the %s only, (number of genes with at least one non-conserved  motif match)','%s, (%s)')" %(jobID,PSSMnum,xlabelList[n],uniqueHitsList[n],uniqueGenesList[n]))
        for n in range(len(speciesList)):                
            cursor.execute("insert into %s_ResultsSummary%s (result,count)  values('percentage of motif matches with at least one %s orthologous motif match','%.2f')" %(jobID,PSSMnum,speciesList[n],SpeciesCountList[n]))
    writeResults1()

def mkPlot3(regionsList,regionCountList):
    plt.figure(3+(PSSMnum-1)*5)
    plt.title('motif (%s) site regions' %(PSSMname))
    ind=np.arange(len(regionsList))
    width=.35*2
    consPlot=plt.bar(ind,regionCountList, width, color = 'b')
    plt.ylabel('site count')
    plt.xticks(ind + .5*width, regionsList,fontsize='8',rotation=-90)
    autolabel(consPlot)
    plt.savefig("%s/%s/Results_plots/3_Site_Regions_%s.png" %(output_dir,PSSMname.split()[0],PSSMname),format = 'png')

def mkPlot4(PSSMnum):
    cursor=conn.cursor()
    def getRegionCount(region,speciesCount):
        cursor.execute("select region_site_count%s from %s_All_conserved_hits%s_ranked where site_region%s ='%s' and region_site_count%s>0" %(speciesCount,jobID,PSSMnum,speciesCount,region,speciesCount))
        count=cursor.fetchall()        
        count=[int(i[0]) for i in count]
        return count
    up=getRegionCount('upstream',1)
    down=getRegionCount('downstream',1)
    exon=getRegionCount('EXON',1)
    intron=getRegionCount('INTRON',1)
    print 'making HitsPerGene plot'
    plt.figure(4+(PSSMnum-1)*5)
    plt.subplot(111)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    axisCoords=plt.axis()
    try:bins=max(up+down+intron+exon)
    except:bins=10
    plt.title('motif (%s) elegans_HitsPerGene%s'%(PSSMname,PSSMnum))
    try:plt.hist(exon,bins, label='exon',histtype='stepfilled',facecolor='red',alpha=1)
    except:None
    try:plt.hist(down,bins,label='downstream',histtype='stepfilled',facecolor='green',alpha=1)    
    except:None
    try:
        plt.hist(up,bins, label='upstream',histtype='stepfilled',facecolor='blue',alpha=1)
    except:None
    try:plt.hist(intron,bins,label='intron',histtype='stepfilled',facecolor='purple',alpha=1)    
    except:None    
    plt.ylabel('gene frequency')
    plt.xlabel('motif frequency')
    plt.legend()
    width=.35
    plt.xticks(np.arange(bins))
    plt.savefig("%s/%s/Results_plots/4_elegans_HitsPerGene%s_%s.png" %(output_dir,PSSMname.split()[0],PSSMnum,PSSMname),format = 'png')

    
def mkPlot5(PSSMnum,allScores,elegScores,brigScores,brenScores,remaScores,japoScores):
    plt.figure(5+(PSSMnum-1)*5)
    plt.subplot(321)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    if len(set(allScores))<13:bins=len(set(allScores))
    else:bins=12
    plt.title('%s Score Dist by Species'%(PSSMname))
    if len(allScores)>0:his1=plt.hist(allScores,bins,facecolor='blue',label='all_species')
    plt.ylabel('score frequency')
    plt.legend(loc=2)
    axisCoords=plt.axis()

    plt.subplot(322)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    plt.axis(axisCoords)
    if len(elegScores)>0:his2=plt.hist(elegScores,bins, facecolor='pink',label='eleg')
    plt.legend(loc=2)

    plt.subplot(323)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    plt.axis(axisCoords,size='xx-small')
    plt.ylabel('score frequency')
    if len(brigScores)>0:his3=plt.hist(brigScores,bins, facecolor='orange',label='brig')
    plt.legend(loc=2)

    plt.subplot(324)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    plt.axis(axisCoords)
    if len(remaScores)>0:his4=plt.hist(remaScores,bins, facecolor='green',label='rema')
    plt.legend(loc=2)

    plt.subplot(325)
    plt.axis(axisCoords)
    if len(japoScores) > 0: his5=plt.hist(japoScores,bins, facecolor='purple', label='japo')
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    plt.ylabel('score frequency')
    plt.xlabel('score')
    plt.legend(loc=2)

    plt.subplot(326)
    plt.yticks(fontsize='8')
    plt.xticks(fontsize = '8')
    plt.axis(axisCoords)
    if len(brenScores)>0:his6=plt.hist(brenScores,bins, facecolor='m',label='bren')
    plt.legend(loc=2)
    plt.xlabel('score')
    plt.savefig("%s/%s/Results_plots/5_ScoreDistBySpecies_%s.png" %(output_dir,PSSMname.split()[0],PSSMname),format = 'png')

def mkPlot6():
    def getList(conservation):
        cursor.execute("select distinct associated_gene1, site_offset1,site_id1 from %s_All_conserved_hits%s_ranked where conservation=%s" %(jobID,PSSMnum,conservation))
        offsets=cursor.fetchall()
        offset1List=[i[1] for i in offsets]
        return offset1List

    plt.figure(6+(PSSMnum-1)*5)
    for i in range(0,5):
        offset1List=getList(i+1)
        plt.subplot(1,1,1)
        plt.yticks(fontsize='8')
        plt.xticks(fontsize = '8')        
        if len(offset1List)>0:
            events,edges,patches=plt.hist(offset1List,histtype='stepfilled',label="%s" %(i+1))
            plt.legend()
    plt.ylabel("site frequency per gene")
    plt.xlabel("distance")
    plt.title("site frequency and distance by conservation")
    plt.savefig("%s/%s/Results_plots/6_Distance_%s.png" %(output_dir,PSSMname.split()[0],PSSMname),format = 'png')

def rankPercentile(score,allscores):
    rank=100-(scipy.stats.percentileofscore(allscores,score*(-1),kind='strict'))
    return rank

def mk_site_region_string(site_region):
    s=site_region.split('\'')
    site_region_string=''
    for n in s:
        if len(n) > 2:site_region_string=site_region_string +' ' + n 
    return site_region_string

def getTopRanked(PSSMnum):
    print 'getting top ranked per gene'
    cursor.execute("delete from %s_top_ranked_hits%s" %(jobID,PSSMnum))
    cursor.execute("select distinct associated_gene1 from %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum))
    eleg_gids=cursor.fetchall()            
    l=float(len(eleg_gids))
    count=0
    for eleg_gid in eleg_gids:        
        count+=1
        if PSSMnum=='1':sys.stdout.write("progress_%s: %d%%   \r" % (PSSMnum,count/l*100) )
        else:sys.stdout.write("%sprogress_%s: %d%%   \r" % ('\t\t\t'*(int(PSSMnum)-1),PSSMnum,count/l*100) )
        sys.stdout.flush()
        #cursor.execute("select * from %s_All_conserved_hits%s_ranked where associated_gene1 = '%s' order by cum_score limit 0,%s" %(jobID,PSSMnum,eleg_gid[0],top_ranked_per_gene))
        #row=cursor.fetchall()
        #print row
        #id=row[0][0]
        #print
        cursor.execute("insert into %s_top_ranked_hits%s select * from %s_All_conserved_hits%s_ranked where associated_gene1='%s' order by rank limit 0, %s" %(jobID,PSSMnum,jobID,PSSMnum,eleg_gid[0],int(top_ranked_per_gene)))
        #cursor.execute("select * from %s_top_ranked_hits%s where associated_gene1 = '%s' " %(jobID,PSSMnum,eleg_gid[0]))
        #row=cursor.fetchall()
        #print row,'top ranked row'
        #print
        #print
        
def mkAll_conserved_hits_ranked(PSSMnum):
    cursor.execute("select * from %s_5_conserved_hits%s" %(jobID,PSSMnum))
    rows=cursor.fetchall()
    for row in rows:
        conservation=5
        coord="%s:%s..%s"%(row[6],row[5],int(row[5])+len(row[3])-1)
        coord2="%s:%s..%s"%(row[16],row[15],int(row[15])+len(row[13])-1)
        coord3="%s:%s..%s"%(row[26],row[25],int(row[25])+len(row[23])-1)
        coord4="%s:%s..%s"%(row[36],row[35],int(row[35])+len(row[33])-1)
        coord5="%s:%s..%s"%(row[46],row[45],int(row[45])+len(row[43])-1)
        avg_species_site_score="%0.4f" %((row[4]+row[14]+row[24]+row[34]+row[44])/conservation)       
        offset_variance ="%0.4f" %( abs(scipy.stats.variation([row[8],row[18],row[28],row[38],row[48]])))

        cursor.execute("insert into %s_All_conserved_hits%s_ranked (ID,conservation,avg_species_site_score,offset_variance,species1,site_id1,site_position1,site_strand1,site_sequence1,site_region1,site_offset1,site_score1,associated_gene1,species2,site_id2,site_position2,site_strand2,site_sequence2,site_region2,site_offset2,site_score2,associated_gene2,species3,site_id3,site_position3,site_strand3,site_sequence3,site_region3,site_offset3,site_score3,associated_gene3,species4,site_id4,site_position4,site_strand4,site_sequence4,site_region4,site_offset4,site_score4,associated_gene4,species5,site_id5,site_position5,site_strand5,site_sequence5,site_region5,site_offset5,site_score5,associated_gene5) values(Null,\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,PSSMnum,conservation,avg_species_site_score,offset_variance, row[0],row[2],coord,row[7],row[3],row[9],row[8],row[4], row[1],row[10],row[12],coord2,row[17],row[13],row[19],row[18], row[14],row[11],row[20],row[22],coord3,row[27],row[23],row[29],row[28], row[24],row[21],row[30],row[32],coord4,row[37],row[33],row[39],row[38], row[34],row[31],row[40],row[42],coord5,row[47],row[43],row[49],row[48], row[44],row[41]))

    cursor.execute("select * from %s_4_conserved_hits%s" %(jobID,PSSMnum))
    rows=cursor.fetchall()
    for row in rows:
        conservation=4
        coord="%s:%s..%s"%(row[6],row[5],int(row[5])+len(row[3])-1)
        coord2="%s:%s..%s"%(row[16],row[15],int(row[15])+len(row[13])-1)
        coord3="%s:%s..%s"%(row[26],row[25],int(row[25])+len(row[23])-1)
        coord4="%s:%s..%s"%(row[36],row[35],int(row[35])+len(row[33])-1)
        avg_species_site_score="%0.4f" %((row[4]+row[14]+row[24]+row[34])/conservation)
        offset_variance ="%0.4f" %( abs(scipy.stats.variation([row[8],row[18],row[28],row[38]])))
        cursor.execute("insert into %s_All_conserved_hits%s_ranked (ID,conservation,avg_species_site_score,offset_variance,species1,site_id1,site_position1,site_strand1,site_sequence1,site_region1,site_offset1,site_score1,associated_gene1,species2,site_id2,site_position2,site_strand2,site_sequence2,site_region2,site_offset2,site_score2,associated_gene2,species3,site_id3,site_position3,site_strand3,site_sequence3,site_region3,site_offset3,site_score3,associated_gene3,species4,site_id4,site_position4,site_strand4,site_sequence4,site_region4,site_offset4,site_score4,associated_gene4) values(Null,\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,PSSMnum,conservation,avg_species_site_score,offset_variance, row[0],row[2],coord,row[7],row[3],row[9],row[8],row[4],row[1],row[10],row[12],coord2,row[17],row[13],row[19],row[18], row[14],row[11],row[20],row[22],coord3,row[27],row[23],row[29],row[28], row[24],row[21],row[30],row[32],coord4,row[37],row[33],row[39],row[38], row[34],row[31]))

    cursor.execute("select * from %s_3_conserved_hits%s" %(jobID,PSSMnum))
    rows=cursor.fetchall()
    for row in rows:
        conservation=3        
        coord="%s:%s..%s"%(row[6],row[5],int(row[5])+len(row[3])-1)
        coord2="%s:%s..%s"%(row[16],row[15],int(row[15])+len(row[13])-1)
        coord3="%s:%s..%s"%(row[26],row[25],int(row[25])+len(row[23])-1)
        avg_species_site_score="%0.4f" %((row[4]+row[14]+row[24])/conservation)
        offset_variance ="%0.4f" %( abs(scipy.stats.variation([row[8],row[18],row[28]])))
        cursor.execute("insert into %s_All_conserved_hits%s_ranked (ID,conservation,avg_species_site_score,offset_variance,species1,site_id1,site_position1,site_strand1,site_sequence1,site_region1,site_offset1,site_score1,associated_gene1,species2,site_id2,site_position2,site_strand2,site_sequence2,site_region2,site_offset2,site_score2,associated_gene2,species3,site_id3,site_position3,site_strand3,site_sequence3,site_region3,site_offset3,site_score3,associated_gene3) values(Null,\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,PSSMnum,conservation,avg_species_site_score,offset_variance, row[0],row[2],coord,row[7],row[3],row[9],row[8],row[4],row[1],row[10],row[12],coord2,row[17],row[13],row[19],row[18], row[14],row[11],row[20],row[22],coord3,row[27],row[23],row[29],row[28], row[24],row[21]))

    cursor.execute("select * from %s_2_conserved_hits%s" %(jobID,PSSMnum))
    rows=cursor.fetchall()
    for row in rows:
        conservation=2
        coord="%s:%s..%s"%(row[6],row[5],int(row[5])+len(row[3])-1)
        coord2="%s:%s..%s"%(row[16],row[15],int(row[15])+len(row[13])-1)
        avg_species_site_score="%0.4f" %((row[4]+row[14])/conservation)
        offset_variance ="%0.4f" %( abs(scipy.stats.variation([row[8],row[18]])))
        cursor.execute("insert into %s_All_conserved_hits%s_ranked (ID,conservation,avg_species_site_score,offset_variance,species1,site_id1,site_position1,site_strand1,site_sequence1,site_region1,site_offset1,site_score1,associated_gene1,species2,site_id2,site_position2,site_strand2,site_sequence2,site_region2,site_offset2,site_score2,associated_gene2) values(Null,\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,PSSMnum,conservation,avg_species_site_score,offset_variance,row[0],row[2],coord,row[7],row[3],row[9],row[8],row[4], row[1],row[10],row[12],coord2,row[17],row[13],row[19],row[18], row[14],row[11]))

    cursor.execute("select * from %s_1_conserved_hits%s" %(jobID,PSSMnum))
    rows=cursor.fetchall()
    for row in rows:
        conservation=1
        coord="%s:%s..%s"%(row[6],row[5],int(row[5])+len(row[3])-1)
        cursor.execute("insert into %s_All_conserved_hits%s_ranked (ID,conservation,species1,site_id1,site_position1,site_strand1,site_sequence1,site_region1,site_offset1,site_score1,associated_gene1) values(Null,\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\',\'%s\')" %(jobID,PSSMnum,conservation,row[0],row[2],coord,row[7],row[3],row[9],row[8],row[4],row[1]))

def rankTable():    
    def scale(list,order):
        if order=='high':
            k=0
            j=1
        else:
            k=1
            j=-1
        try:
            scaled=[((k+j*float(i)-min(list))/(max(list)-min(list))*100) for i in list]            
        except:
            trimmed=[]
            for i in list:                    
                if i !=None and i !='':trimmed.append(k+j*float(i))
                else:trimmed.append(0)
            minL=min(trimmed)
            maxL=max(trimmed)
            if maxL-minL==0:scaled=[100 for i in list]
            else:
                scaled=[]
                for i in list:
                    if i==None or i=='':scaled.append(0)
                    else:scaled.append((k+j*float(i)-minL)/(maxL-minL)*100)
        return scaled

    def rankData(data):         
        uni=list(set(data))            
        uni.sort()
        uni.reverse()
        uniDic={}
        for i in range(len(uni)):
            uniDic[uni[i]]=i+1
        rankedList=[uniDic[n] for n in data]
        return rankedList
            
    cursor=conn.cursor()        
    cursor.execute("select count(*) from %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum))#conservation
    count=cursor.fetchone()
    d=np.zeros((12,int(count[0])))
    
    def eleg_reg_criteriaABC():
        if float(scale_A)!=0:
            print 'eleg site score'
            print 'ranking 1 of 12 lists'
            cursor.execute("select site_score1 from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#site score
            AvgSp=cursor.fetchall()
            cursor.execute("select site_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#site score
            eleg=cursor.fetchall()
            l=AvgSp+eleg
            d[1]=scale([i[0] for i in AvgSp+eleg],'high')            
        if float(scale_B)!=0:
            print 'eleg avg_region_score'
            print 'ranking 2 of 12 lists'
            cursor.execute("select avg_region_score1 from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg region site score
            AvgSp=cursor.fetchall()
            cursor.execute("select avg_region_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg region site score
            eleg=cursor.fetchall()
            d[2]=scale([i[0] for i in AvgSp+eleg],'high')    
        if float(scale_C)!=0:
            print 'eleg region site count'
            print 'ranking 3 of 12 lists'
            cursor.execute("select region_site_count1 from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#region site count
            AvgSp=cursor.fetchall()
            cursor.execute("select region_site_count1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#region site count
            eleg=cursor.fetchall()
            d[4]=scale([i[0] for i in AvgSp+eleg],'high')
        
    def eleg_gene_criteriaDE():
        if float(scale_D)!=0:
            print 'eleg gene site count'
            print 'ranking 4 of 12 lists'
            cursor.execute("select gene_site_count1 from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#gene_site_count
            AvgSp=cursor.fetchall()
            cursor.execute("select gene_site_count1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#gene_site_count
            eleg=cursor.fetchall()
            d[5]=scale([i[0] for i in AvgSp+eleg],'high')    
        if float(scale_E)!=0:
            print 'eleg avg_gene_score'
            print 'ranking 5 of 12 lists'
            cursor.execute("select avg_gene_score1 from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg gene score          
            AvgSp=cursor.fetchall()
            cursor.execute("select avg_gene_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg gene score          
            eleg=cursor.fetchall()
            d[3]=scale([i[0] for i in AvgSp+eleg],'high')

    def avg_species_reg_criteriaFG():
        
        if float(scale_F)!=0:
            print 'conservation'
            print 'ranking 6 of 7 lists'
            cursor.execute("select conservation from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#conservation
            AvgSp=cursor.fetchall()
            cursor.execute("select conservation from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#conservation
            eleg=cursor.fetchall()
            print eleg[:10],'eleg cons'
            result=scale([int(i[0]) for i in AvgSp+eleg],'high')        
            d[0]=[float(scale_F)*r for r in result]
            
        if float(scale_G)!=0:
            print 'ranking 10 of 12'
            cursor.execute("select offset_variance from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#offset_variance
            AvgSp=cursor.fetchall()
            cursor.execute("select offset_variance from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#offset_variance
            eleg=cursor.fetchall()
            l=AvgSp+eleg
            maxS=max(l)
            eleg=[maxS for i in eleg]
            result=scale([i[0] for i in AvgSp+eleg],'low')
            d[11]=[float(scale_G)*r for r in result]

    def avg_species_gene_criteriaDE():
        if float(scale_D)!=0:
            print 'ranking 4 of 7'
            cursor.execute("select avg_species_gene_site_count from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg species gene site count        
            AvgSp=cursor.fetchall()
            cursor.execute("select gene_site_count1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg species gene site count        
            eleg=cursor.fetchall()
            l=AvgSp+eleg#set non conserved site scores to min
            minS=min(l)#
            eleg=[minS for i in eleg]#
            result=scale([i[0] for i in AvgSp+eleg],'high')
            d[10]=[float(scale_D)*r for r in result]
        if float(scale_E)!=0:
            print 'ranking 5 of 7'
            cursor.execute("select avg_species_gene_score from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg species gene score
            AvgSp=cursor.fetchall()
            cursor.execute("select avg_gene_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg species gene score
            eleg=cursor.fetchall()
            l=AvgSp+eleg#set non conserved site scores to min
            minS=min(l)#
            eleg=[minS for i in eleg]#
            result=scale([i[0] for i in AvgSp+eleg],'high')
            d[8]=[float(scale_E)*r for r in result]
        
    def avg_species_reg_criteriaABC():
        if float(scale_A)!=0:
            print 'ranking 1 of 7 lists'
            cursor.execute("select avg_species_site_score from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg species region score        
            AvgSp=cursor.fetchall()
            cursor.execute("select site_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg species region score        
            eleg=cursor.fetchall()
            l=AvgSp+eleg#set non conserved site scores to min
            minS=min(l)#
            eleg=[minS for i in eleg]#
            result=scale([i[0] for i in AvgSp+eleg],'high')
            d[6]=[float(scale_A)*r for r in result]
        if float(scale_B)!=0:
            print 'ranking 2 of 7'
            cursor.execute("select avg_species_region_score from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg species region score
            AvgSp=cursor.fetchall()        
            cursor.execute("select avg_region_score1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg species region score
            eleg=cursor.fetchall()
            l=AvgSp+eleg#set non conserved site scores to min
            minS=min(l)#
            eleg=[minS for i in eleg]#
            result=scale([i[0] for i in AvgSp+eleg],'high')        
            d[7]=[float(scale_B)*r for r in result]
        if float(scale_C)!=0:
            print 'ranking 3 of 7'
            cursor.execute("select avg_species_region_site_count from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))#avg species region site count
            AvgSp=cursor.fetchall()
            cursor.execute("select region_site_count1 from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))#avg species region site count
            eleg=cursor.fetchall()
            l=AvgSp+eleg#set non conserved site scores to min
            minS=min(l)#
            eleg=[minS for i in eleg]#
            result=scale([i[0] for i in AvgSp+eleg],'high')
            d[9]=[float(scale_C)*r for r in result]


    #scale_A=1
    #scale_B=0
    #scale_C=0
    #scale_D=0
    #scale_F=0
    #scale_G=0
    
    avg_species_reg_criteriaFG()
    avg_species_gene_criteriaDE()
    avg_species_reg_criteriaABC()
    #eleg_reg_criteriaABC()
    #eleg_gene_criteriaDE()        

    div=7
    div1=6    
    
    print 'summing ranked lists'
    ranksList=[sum(i) for i in d.T]
    
    cursor.execute("select * from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum))
    AvgSp_rows=cursor.fetchall()
    cursor.execute("select * from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum))
    eleg_rows=cursor.fetchall()
    rows=AvgSp_rows+eleg_rows
    print 'averaging normalized sums'
    count=0
    ranksAdjusted=[]
    s=1#scale factor to weight unconserverd site cum scores.    
    for row in rows:
        if int(row[11])>1:#conservation
            ranksAdjusted.append(ranksList[count]/div)
        else:ranksAdjusted.append(ranksList[count]/div1/s)
        count+=1
    print 'assigning rank order to normalized averages'
    ranksList=rankData(ranksAdjusted)
    cum_scoreList=ranksAdjusted
    command="select ID from %s_All_conserved_hits%s_ranked where conservation>1" %(jobID,PSSMnum)
    cursor.execute(command)
    AvgSp_IDs=cursor.fetchall()
    command="select ID from %s_All_conserved_hits%s_ranked where conservation=1" %(jobID,PSSMnum)
    cursor.execute(command)
    eleg_IDs=cursor.fetchall()
    IDs=AvgSp_IDs+eleg_IDs
    count=0
    for ID in IDs:
        cursor.execute("update %s_All_conserved_hits%s_ranked set rank = %s where ID = %s" %(jobID,PSSMnum,ranksList[count],ID[0]))
        cursor.execute("update %s_All_conserved_hits%s_ranked set cum_score = %s where ID= %s" %(jobID,PSSMnum,cum_scoreList[count],ID[0]))
        count+=1
    
def update():
    def update_site_counts_avg_score(region,speciesCount):
        if region=='All':
            s=''
            s2=''
        else: 
            s='and site_region%s =\'%s\' ' %(speciesCount,region)
            s2=' and site_region=\'%s\' ' %(region)
        gids=[gid for gid in cursor.execute("select distinct associated_gene%s,species%s from %s_All_conserved_hits%s_ranked where species%s not null %s" %(speciesCount,speciesCount,jobID,PSSMnum,speciesCount,s))]
        for gid in gids:            
            gene_id=gid[0]
            species=gid[1]            
            id_scores=[id for id in cursor.execute("select distinct site_id,score from %s_%s_hit_ortho%s where gene_name = '%s' %s" %(jobID,species,PSSMnum,gene_id,s2))]
            count=len(id_scores)
            avg_gene_score=sum([float(id[1]) for id in id_scores])/count
            if region=='All':
                cursor.execute("update %s_All_conserved_hits%s_ranked set avg_gene_score%s=%s where associated_gene%s = '%s'" %(jobID,PSSMnum,speciesCount,avg_gene_score,speciesCount,gid[0]))
                cursor.execute("update %s_All_conserved_hits%s_ranked set gene_site_count%s=%s where associated_gene%s = '%s'" %(jobID,PSSMnum,speciesCount,count,speciesCount,gid[0]))
            else:
                cursor.execute("update %s_All_conserved_hits%s_ranked set avg_region_score%s=%s where associated_gene%s = '%s' %s"%(jobID,PSSMnum,speciesCount,avg_gene_score,speciesCount,gid[0],s))
                cursor.execute("update %s_All_conserved_hits%s_ranked set region_site_count%s=%s where associated_gene%s = '%s' %s" %(jobID,PSSMnum,speciesCount,count,speciesCount,gid[0],s))
           
                
    def updateAvgSpecies(k,v):
        cons2=" + %s2" %(v) 
        cons3=" %s +%s3" %(cons2,v)
        cons4=" %s + %s4" %(cons3,v)
        cons5=" %s + %s5" %(cons4,v)
        for cons in [cons2,cons3,cons4,cons5]:
            command="update %s_All_conserved_hits%s_ranked set %s = (%s1 %s)/conservation where conservation=%s" %(jobID,PSSMnum,k,v,cons,cons[-1])
            cursor.execute(command)

    regions=['All','upstream','downstream','INTRON','EXON']
    for region in regions:
        [update_site_counts_avg_score(region,speciesCount+1) for speciesCount in range(5)]
        
    for k,v in {'avg_species_region_score':'avg_region_score','avg_species_gene_site_count':'gene_site_count','avg_species_gene_score':'avg_gene_score','avg_species_region_site_count':'region_site_count'}.iteritems():
        updateAvgSpecies(k,v)
    
def getRegions():
    print 'gathering conserved hits regions'
    cursor=conn.cursor()
    regionList,regionsList,regionCountList=['upstream','downstream','INTRON','EXON'],[],[]
    for region in regionList:
        regionCount=0
        cursor.execute("select count(distinct site_id1) as count from  %s_All_conserved_hits%s_ranked  where site_region1='%s'" %(jobID,PSSMnum,region))
        numrows=cursor.fetchone()[0]
        regionCount=regionCount+numrows        
        if regionCount>0:
            regionsList.append(region)
            regionCountList.append(regionCount)
    return regionsList,regionCountList

def getSpeciesCounts():
    cursor=conn.cursor()
    cursor.execute("select distinct associated_gene1 from %s_All_conserved_hits%s_ranked" %(jobID,PSSMnum))
    gids=cursor.fetchall()
    brigCount=brenCount=remaCount=japoCount=elegCount=0
    for gid in gids:
        cursor.execute("select distinct species1,species2,species3,species4,species5 from %s_All_conserved_hits%s_ranked where associated_gene1=\'%s\'" %(jobID,PSSMnum,gid[0]))
        results=cursor.fetchall()
        eleg=brig=bren=japo=rema=0
        for i in range(len(results)):
            for j in range(5):
                if results[i][j]!=None:
                    if str(speciesList[1]) in results[i][j]:brig=1.0
                    if str(speciesList[2]) in results[i][j]:bren=1.0
                    if str(speciesList[3]) in results[i][j]:rema=1.0
                    if str(speciesList[4]) in results[i][j]:japo=1.0
                    if str(speciesList[0]) in results[i][j]:eleg=1.0
        brigCount+=brig
        brenCount+=bren
        remaCount+=rema
        japoCount+=japo
        elegCount+=eleg
    total=elegCount
    if total==0:SpeciesCountList=[0,0,0,0,0]
    else:SpeciesCountList=[100,brigCount/total*100,brenCount/total*100,remaCount/total*100,japoCount/total*100]
    return SpeciesCountList

def countCons():
    def getCount(conservation,attribute):
        rows=[row for row in cursor.execute("select distinct conservation, %s from %s_All_conserved_hits%s_ranked where conservation=%s" %(attribute,jobID,PSSMnum,conservation))]
        return len(rows)
    uniqueGenesList=[getCount(n+1,'associated_gene1') for n in range(5)]
    uniqueHitsList=[getCount(n+1,'site_id1') for n in range(5)]
    return uniqueGenesList,uniqueHitsList
def getScores():
    def getSpeciesScores(species):
        cursor.execute("select score from %s_%s_hits%s" %(jobID,species,PSSMnum))
        scores=cursor.fetchall()
        return [score[0] for score in scores]
    print 'getting eleg scores'
    elegScores=getSpeciesScores(speciesList[0])
    print 'getting brig scores'
    brigScores=getSpeciesScores(speciesList[1])
    print 'getting bren scores'
    brenScores=getSpeciesScores(speciesList[2])
    print 'getting japo scores'
    japoScores=getSpeciesScores(speciesList[4])
    print 'getting rema scores'
    remaScores=getSpeciesScores(speciesList[3])
    print 'combining all scores'
    allScores=elegScores+brigScores+brenScores+japoScores+remaScores
    return allScores,elegScores,brigScores,brenScores,remaScores,japoScores

def main():    
    
    prepare_sql_tables()                  
    
    print 'making All conserved hits table'
    mkAll_conserved_hits_ranked(PSSMnum)    
    print 'updating site counts and scores'
    
    update()
    print 'ranking conserved hits'
    
    rankTable()
    print 'making top ranked hits table'
    
    getTopRanked(PSSMnum)    
    
    allScores,elegScores,brigScores,brenScores,remaScores,japoScores=getScores()
    print 'getting species counts'
    SpeciesCountList=getSpeciesCounts()
    print 'getting conservered hits counts'
    uniqueGenesList,uniqueHitsList=countCons()        
    print 'making plots 1,2, and 3'
    mkPlot1_2(PSSMnum,uniqueHitsList,uniqueGenesList,SpeciesCountList)       
    print 'making score distributions plot'
    mkPlot5(PSSMnum,allScores,elegScores,brigScores,brenScores,remaScores,japoScores)
    print 'getting regions lists'
    
    regionsList,regionCountList=getRegions()
    print regionsList,regionCountList,'regionsList,regionCountList'
    print 'making plot 3'
    mkPlot3(regionsList,regionCountList)
    print 'making HitsPerGene plot'
    mkPlot4(PSSMnum)
    print 'making plot6'
    mkPlot6()
    
main()
