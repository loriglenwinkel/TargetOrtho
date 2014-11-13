#!/Usr/bin/python
import sys,os
from optparse import OptionParser
import sqlite3 as lite
parser = OptionParser()
parser.add_option('-i', '--output_dir', dest='output_dir',help='-d /data/newTargetOrtho/run/TargetOrtho_output/123232')
parser.add_option('-j', '--jobID',dest='jobID',help='-g 2012419')
parser.add_option('-g', '--jobTag',dest='jobTag',help='-j job_name')
parser.add_option( '-c', '--TargetOrtho_path', dest='TargetOrtho_path', help='-c /data/newTargetOrtho/run' )
parser.add_option('-m','--matrix_count',dest='matrix_count',help='matrix_count: -m 3')
parser.add_option('-n','--matrixNames',dest='matrixNamesList',help='-n name1*name2*name3')
parser.add_option('-q','--file',dest='queryFile_path',help='-f queryFile.txt')
parser.add_option('-w','--QueryOnly',dest='QueryOnly', help='=w True')
parser.add_option('-S','--speciesList',dest='speciesList',help='-S c_eleg-c_brig-c_bren-')
(options, args) = parser.parse_args()

speciesList=options.speciesList
speciesList=speciesList.split('-')[:-1]
QueryOnly=options.QueryOnly
queryFile_path=options.queryFile_path
matrix_count=int(options.matrix_count)
TargetOrtho_path=options.TargetOrtho_path
output_dir=options.output_dir
jobID=options.jobID
jobTag=options.jobTag
matrixNamesList=options.matrixNamesList.split('*')[:-1]
conn = lite.connect('%s/run/sqlite_tmp_files/%s_TargetOrtho.db' %(TargetOrtho_path,jobID),isolation_level=None)
cursor = conn.cursor()    
NameDic={'c_eleg':('C.&nbsp;elegans','c_elegans'),
         'c_bren':('C.&nbsp;brenneri','c_brenneri'),
         'c_brig':('C.&nbsp;briggae','c_briggsae'),
         'c_rema':('C.&nbsp;remanei','c_remanei'),
         'c_japo':('C.&nbsp;japonica','c_japonica'),
         'd_mel':('D.&nbsp;melanogaster','d_melanogaster'),
         'd_sec':('D.&nbsp;sechellia','d_sechellia'),
         'd_sim':('D.&nbsp;simulans','d_simulans'),
         'd_yak':('D.&nbsp;yakuba','d_yakuba'),
         'd_ere':('D.&nbsp;erecta','erecta')}

print speciesList,'speciesList'
def color(seq):
    colorDic={'G':'orange\">G','A':'red\">A','C':'blue\">C','T':'green\">T','N':'black\">N','X':'black\">X'}
    final_str=''
    for j in seq:
        fontStr='<span style=\"color:%s</span>' %colorDic[j]
        final_str=final_str+fontStr
    return final_str

def mvTablesFromSqltoTEXT(table,m,db):
    if db!='':cursor.execute("attach '%s' as a" %(db))
    if m!='':dir="%s/%s/Results_text" %(output_dir,matrixNamesList[m-1])
    
    else:dir=output_dir
    def exportTable(table):
        outfile=file('%s/%s.txt' %(dir,table),'w')
        [outfile.write('%s\t' %str(i).strip()) for i in [col[1] for col in cursor.execute("PRAGMA table_info(%s)" %table)]]
        outfile.write('\n')
        if ('top_ranked_hits' in table or 'All_conserved' in table or 'QueryList' in table) and ('CRM' not in table):
            [([outfile.write('%s\t' %str(i).strip()) for i in row],outfile.write('\n')) for row in cursor.execute("select * from %s order by rank ASC " %(table))]
        else:
            [([outfile.write('%s\t' %str(i).strip()) for i in row],outfile.write('\n')) for row in cursor.execute("select *  from %s" %(table))]
    exportTable(table)
    if db!='':cursor.execute("detach a")

def galaxy_results_file(matrix_count):
    outfile=file('%s/index.html' %(output_dir),'w')
    outfile.write("<html>\n<head>\n<title>TargetOrthoV1.0 Results</title>\n<h1>TargetOrtho v1.0 Results</h1>\n<p>\n\
                  job name: %s <br>\n\
                  jobID: %s <br>\n\
                  motif names: %s<br>\n</p>\n</head>\n<body>\n<p>\n" %(jobTag,jobID,matrixNamesList))
    outfile.write("<a href=\"Download_%s.tar.gz\"><img src=\"download_icon.png\" alt=\"Download full results directory\" width=\"100\" height=\"100\"></a><br>" %(jobID))
    outfile.write("<a href=\"Download_%s.tar.gz\">Download full results directory (tar.gz) </a><br><br>\n" %(jobID))
    pref=''
    if matrix_count>1:
        ref1=''
        ref2='<a href=\"%s_CRM_results.html\">CRM_results.html (whole genome)</a>/<a href=\"%s_CRM_results.txt\">(text)</a>:<br> The subset of motif matches and their associated genes in which at least one motif match for each input motif is found in the same gene region. Each line in this table shows a unique combination of motif matches per gene. Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_OACisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria></a> for an explanation of column headers. <br><br>\n' %(jobID,jobID)


        ref3='<a href=\"%s_top_ranked_per_gene_CRM.html\">CRM_results_top_ranked_per_gene.html</a>/<a href=\"%s_top_ranked_per_gene_CRM.txt\">(text)</a>:<br> Data for the top ranked motif match per gene. For the complete set of motif match data, see the CRM_results table or the QueryListCRM table if the Query Only option was True. Each row  shows data for the best ranked motif match per gene.  Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers. <br><br><br><br>\n' %(jobID,jobID)
        if queryFile_path != 'None' and queryFile_path!=None:            
            ref1='<a href=\"%s_QueryListResults_CRM.html\">CRM_QueryList_results.html</a>/<a href=\"%s_QueryListResults_CRM.txt\">(text)</a>:<br> The subset of result CRM data in which each associated gene is present in the user input query list file. The subset of motif matches and their associated genes in which at least one motif match for each input motif is found in the same gene region. Each line in this table shows a unique combination of motif matches per gene. Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers. <br><br>\n' %(jobID,jobID)
            if QueryOnly=='True':ref2=''
        outfile.write('<h2>Cis Regulatory Module (CRM) Results</h2>\n')
        outfile.write('%s\n%s\n%s\n' %(ref1,ref2,ref3))
        outfile.write('<h3>Individual motif results</h3>\n')        
        pref='.'
    for n in range(matrix_count):
        outfile.write('<B>Results for motif %s</B>: %s<br>\n' %(n+1,matrixNamesList[n]))
        ref1=''
        ref2='<a href=\"%s%s_All_conserved_hits%s_ranked.html\">All_conserved_hits_ranked.html (whole genome)</a>/<a href=\"%s%s_All_conserved_hits%s_ranked.txt\">(text)</a>:<br> Each row shows data associated with one motif match so that data for any  gene associated with multiple motif matches will be shown in separate rows. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers. <br><br>\n' %(pref,jobID,n+1,pref,jobID,n+1)
        ref3='<a href=\"%s%s_top_ranked_hits%s.html\">Top_ranked_hits_per_gene.html</a>/<a href=\"%s%s_top_ranked_hits%s.txt\">(text)</a>:<br>Data for the top ranked motif match per gene. For the complete set of motif match data, see the All conserved hits ranked table (or the QueryList Results table if the query only option was used). Each row  shows data for the best ranked motif match per gene.  See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers. <br><br>\n' %(pref,jobID,n+1,pref,jobID,n+1)


        if queryFile_path != 'None' and queryFile_path!=None:
            ref1='<a href=\"%s%s_QueryListResults%s.html\">Query List Results.html</a>/<a href=\"%s%s_QueryListResults%s.txt\">(text)</a>:<br>The subset of result data in which each associated gene is present in the user input query list file.Each row shows data associated with one motif match so that data for any  gene associated with multiple motif matches will be shown in separate rows. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers. <br><br>\n' %(pref,jobID,n+1,pref,jobID,n+1)


            if QueryOnly=='True':ref2=''
        outfile.write("%s%s%s\n" %(ref1,ref2,ref3))
        outfile.write("<a href=\"%s%s_%s_%s.bed\">GenomeBrowserFile.bed</a>:<br> genome browser track file. Each site is assigned a color score correlated with the log-odds score of the motif match so that stronger predicted binding sites are colored darker than weaker binding sites. Each motif match score is normalized between 0 and 1000 where 1000 represents the strongest possible binding site score.  <br><br><a href=\"%s%s_ResultsSummary%s.html\">ResultsSummary.html</a>:<br> Summary of TargetOrtho results<br><br>\n" %(pref,jobID,speciesList[0],matrixNamesList[n],pref,jobID,n+1))
    outfile.write("<a href=\"%s_inputSummary.txt\">inputSummary.txt</a>:<br> Summary of TargetOrtho input parameters<br><br>\n" %(jobID))
    outfile.write("<a href=\"%sTargetOrthoStdout.out\">job_output.txt</a>: Standard output from job execution<br>\n" %(jobID))
    outfile.write("If you use TargetOrtho, please cite:\n Lori Glenwinkel, Di Wu, Gregory Minevich and Oliver Hobert. TargetOrtho: A Phylogenetic Footprinting Tool to Identify Transcription Factor Targets. Genetics May 1, 2014 vol. 197 no. 1 61-76\n")
def mkHTMLtables(table,m,db):
    if db !='':cursor.execute("attach '%s' as c" %(db))
    if m !='':dir="%s/%s/Results_html" %(output_dir,matrixNamesList[m-1])
    else:dir=output_dir
    outfile=file('%s/%s.html' %(dir,table),'w')
    if ('CRM' in table or 'per_gene_CRM' in table or 'QueryListResults_CRM' in table) and (matrix_count>1):
        if matrix_count >2:rank3=' +motif3_site_rank'
        else:rank3=''
        if matrix_count >3:rank4=' +motif4_site_rank'
        else:rank4=''
        if matrix_count >4:rank5=' +motif5_site_rank'
        else:rank5=''
        
        
        command="""select * from %s order by (motif1_site_rank+motif2_site_rank%s%s%s) ASC"""%(table,rank3,rank4,rank5)        
        rows=list([row for row in cursor.execute(command)])
    else:
        if ('top_ranked_hits' in table or 'QueryList' in table or 'All_conserved' in table) and ('CRM' not in table):
            command="""select * from %s order by rank ASC""" %table
            rows=list([row[1:] for row in cursor.execute(command)])
        else:
            command="""select * from %s""" %table
            rows=list([row for row in cursor.execute(command)])
        
    outfile.write("<html>\n<head>\n<title>TargetOrthoV1.0 Results</title>\n<script src=\".sorttable.js\"></script>\n")
    outfile.write("<script type=\"text/javascript\" src=\".jquery-1.8.2.min.js\"></script>\n")
    num=table[-1]
    try:matrixName="matrix_name: %s" %(matrixNamesList[int(num)-1])
    except:matrixName=''    
    t=table
    #assign table descriptions for html headers
    if 'CRM' in t:
        desc=''
        if 'if top' in t:desc='Data for the top ranked motif match per gene. For the complete set of motif match data, see the CRM_results table or the QueryListCRM table if the Query Only option was True. Each row  shows data for the best ranked motif match per gene.  Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'        
        elif 'Query' in t:desc='The subset of result CRM data in which each associated gene is present in the user input query list file. The subset of motif matches and their associated genes in which at least one motif match for each input motif is found in the same gene region. Each line in this table shows a unique combination of motif matches per gene. Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'
        else:
            desc='The subset of motif matches and their associated genes in which at least one motif match for each input motif is found in the same gene region. Each line in this table shows a unique combination of motif matches per gene. Each input motif match rank is shown where the ranks is based on the cumulative site score for the motif match in the CRM. Results are ordered by average motif rank. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'
    else:    
        if 'ResultsSummary' in t:desc='Summary of TargetOrtho results'
        if 'top_ranked' in t:desc='Data for the top ranked motif match per gene. For the complete set of motif match data, see the All conserved hits ranked table.. Each row  shows data for the best ranked motif match per gene.  See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'        
        if 'Query' in t:desc='The subset of result data in which each associated gene is present in the user input query list file.Each row shows data associated with one motif match so that data for any  gene associated with multiple motif matches will be shown in separate rows. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers'
        if 'hit_ortho' in t:desc=''
        if 'All' in t:desc='Each row shows data associated with one motif match so that data for any  gene associated with multiple motif matches will be shown in separate rows. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'
        else:
            if 'conserved_hits' in t:desc='Data species for motif matches conserved in 1,2,3,4 or 5 species. Each row shows data associated with one motif match so that data for any  gene associated with multiple motif matches will be shown in separate rows. See <a href="http://hobertlab.org/TargetOrtho/5_CisOrtho_RankingCriteria.pdf" target="_blank">TargetOrtho Ranking Criteria</a> for an explanation of column headers.'

    outfile.write("\
 <script type=\"text/javascript\">\n\
  $(document).ready(function(){\n\
   //hide all the rows that have class=item\n\
   $(\".ITEM\").hide();\n\
   \n\
   //Slide up and down on click\n\
   $(\"#data a.toggle\").click(function(){\n\
    var toToggle = \".\"+$(this).attr(\"href\");\n\
    if($(this).closest(\"tr\").hasClass(\"rowhighlight\")){\n\
     $(this).closest(\"tr\").removeClass(\"rowhighlight\");\n\
    }\n\
    else{\n\
     $(this).closest(\"tr\").addClass(\"rowhighlight\");\n\
    }\n\
    $(toToggle).slideToggle(\"slow\");\n\
    return false;\n\
   });\n\
  });\n\
 </script>\n\
 <style type=\"text/css\">\n\
  .hide { display: none; }\n\
  .show { display: block; }\n\
  .row0 { background-color:#ddd; }\n\
  .row1 { background-color:#eee; }\n\
  .rowhighlight { background-color:#ddd;}\n\
  .rowhighligt, .rowhighlight a { color: blue; }\n\
  tr:hover, tr.atr:hover { background-color :yellow; }\n\
 </style>\n\
 <h1>TargetOrtho v1.0 Results</h1>\n\
 <p>\n\
                  <B>table name</B>: %s<br>\n\
                  <B>description</B>: %s<br>\n\
                  <B>%s</B><br>\n\
                  <B>job name</B>: %s<br>\n\
                  <B>jobID: %s</B><br>\n\
 </p>\n\
</head>\n\
<body>\n" %(t[len(jobID)+1:],desc,matrixName,jobTag,jobID))    
    if '%s_CRM_results' %(jobID) in table or 'QueryListResults_CRM' in table or 'top_ranked_per_gene_CRM' in table:                
        for n in range(matrix_count):
            outfile.write('<th><br><B>motif%s</B>:%s </br></th>\n' %(n+1,matrixNamesList[n]))    
    outfile.write('<table class="sortable" id="data" border="0" style="width: 100%">\n')
    outfile.write('<tr>\n')    
    command="""pragma table_info(%s)""" %table
    fields=[col[1] for col in cursor.execute(command)]
    n=len(fields)
    headers=fields
    if '%s_All_conserved_hits' %(jobID) in table or '%s_top_ranked_hits' %(jobID) in table  or '%s_QueryListResults' %(jobID) in table:
        if 'CRM' not in table:
            outfile.write('<th>%s</th>\n' %(''))
            n=16
            headers=fields[1:n+1]+fields[-6:]                               
    for field in headers:
        field=field.strip()
        if ('%s_All_conserved_hits' %(jobID) in table or '%s_top_ranked_hits' %(jobID) in table  or '%s_QueryListResults' %(jobID) in table) and ('CRM' not in table):
            if field.endswith('1'):outfile.write('<th>%s</th>\n' %(field[:-1]))
            else:outfile.write('<th>%s</th>\n' %(field))
        else:outfile.write('<th>%s</th>\n' %(field))
    outfile.write('</tr>\n')    
    parentID=0
    def writeCommon(row):
        #write html rows that will have blank exandable rows beneath        
        for j in row[-6:]:                        
            try:outfile.write('<td>%.2f</td>'%(float(j)))                
            except:outfile.write('<td>%s</td>'%(j))                
    def write_child_rows(c):
        #write expandable html rows for each additional species
        NameDic={'c_bren':('C.&nbsp;brenneri','c_brenneri'),
                 'c_brig':('C.&nbsp;briggae','c_briggsae'),
                 'c_rema':('C.&nbsp;remanei','c_remanei'),
                 'c_japo':('C.&nbsp;japonica','c_japonica'),
                 'd_mel':('D.&nbsp;melanogaster','d_melanogaster'),
                 'd_sec':('D.&nbsp;sechellia','d_sechellia'),
                 'd_sim':('D.&nbsp;simulans','d_simulans'),
                 'd_yak':('D.&nbsp;yakuba','d_yakuba'),
                 'd_ere':('D.&nbsp;erecta','d_erecta')}

        if i==c:
            outfile.write('\n<tr class="ITEM row1 WOW%s">\n' %(parentID))
            outfile.write('<td></td><td><i>%s</i></td><td></td></td><td></td>' %(NameDic[row[i].strip()][0]))#species
        if i>c and i<c+13:
            if i==c+1:
                if speciesList[0] =='d_mel':outfile.write('<td><a href="http://flybase.org/reports/%s.html" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))#associated_gene
                else:outfile.write('<td><a href="http://www.wormbase.org/db/get?name=%s;class=Gene" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))#associated_gene
            elif i==c+3:
                if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(row[c][0],row[c][2:],row[i],row[i]))#site_position
                else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[row[c].strip()][1],row[i],row[i]))#site_position
            elif i==c+5:
                outfile.write('<td>%s</td>' %(color(row[i])))#site sequence
            elif i==c+8:
                outfile.write('<td></td><td>%s</td>' %(row[i]))#conservation
            elif i==c+12:
                outfile.write("<td>%s</td>%s" %(row[i],'<td></td>'*6))#gene_site_count
            else:outfile.write('<td>%s</td>' %(row[i]))
    def mkTable2(row,i):
        if i==0:
            outfile.write('<tr class="PARENT row0">\n')
            if speciesList[0] =='d_mel':outfile.write('<td><a href="http://flybase.org/reports/%s.html" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))#associated_gene
            else:outfile.write('<td><a href="http://www.wormbase.org/db/get?name=%s;class=Gene" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))                    
        else:
            if matrix_count==2:
                if i==7 or i==8:
                    outfile.write('<td >%s</td>' %(color(row[i])))
                elif i==9 or i==10:
                    if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(speciesList[0][0],speciesList[0][2:],row[i],row[i]))#site_position
                    else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[speciesList[0]][1],row[i],row[i]))
                else:outfile.write('<td>%s</td>' %(row[i]))                
            elif matrix_count==3:
                if i==10 or i==11 or i==12:
                    outfile.write('<td >%s</td>' %(color(row[i])))
                elif i==13 or i==14 or i==15:
                    if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(speciesList[0][0],speciesList[0][2:],row[i],row[i]))#site_position
                    else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[speciesList[0]][1],row[i],row[i]))
                else:outfile.write('<td>%s</td>' %(row[i]))                
            elif matrix_count==4:
                if i==14 or i==15 or i==16 or i==17:
                    outfile.write('<td >%s</td>' %(color(row[i])))
                elif i==18 or i==19 or i==20 or i==21:
                    if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(speciesList[0][0],speciesList[0][2:],row[i],row[i]))#site_position
                    else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[speciesList[0]][1],row[i],row[i]))
                else:outfile.write('<td>%s</td>' %(row[i]))                
            elif matrix_count==5:
                if i==16 or i==17 or i==18 or i==19 or i==20:
                    outfile.write('<td >%s</td>' %(color(row[i])))
                elif i==21 or i==22 or i==23 or i==24 or i==25:
                    if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(speciesList[0][0],speciesList[0][2:],row[i],row[i]))#site_position
                    else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[speciesList[0]][1],row[i],row[i]))
                else:outfile.write('<td>%s</td>' %(row[i]))                
            if i==len(fields):outfile.write('</tr>\n')
    def mkTable(row,i):    
        if i==0:
            outfile.write('<tr class="PARENT row0">\n')
            outfile.write('<td><a href="WOW%s" STYLE="TEXT-DECORATION: NONE" alink="#FF99CC" class="toggle">+</a></td><td><i>%s</i></td>' %(parentID,NameDic[speciesList[0]][0])) 
        elif i>0 and  i <n:     
            if i==3:
                if speciesList[0] =='d_mel':outfile.write('<td><a href="http://flybase.org/reports/%s.html" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))#associated_gene
                else:outfile.write('<td><a href="http://www.wormbase.org/db/get?name=%s;class=Gene" target="_blank"><i>%s</i></a></td>' %(row[i],row[i]))
            elif i==5:
                if 'd_mel' in speciesList[0]:outfile.write('<td><a href="http://flybase.org/cgi-bin/gbrowse/%s%s/?name=%s" target="_blank">%s</a></td>' %(speciesList[0][0],speciesList[0][2:],row[i],row[i]))#site_position
                else:outfile.write('<td><a href="http://www.wormbase.org/tools/genome/gbrowse/%s/?name=%s" target="_blank">%s</a></td>' %(NameDic[speciesList[0]][1],row[i],row[i]))
            elif i==7:outfile.write('<td >%s</td>' %(color(row[i])))
            elif i in [1,4,9,10,11]:outfile.write('<td>%s</td>'%(row[i]))                            
            else:
                try:outfile.write('<td>%.2f</td>'%float((row[i])))
                except:outfile.write('<td>%s</td>'%(row[i]))                
        if i==n:writeCommon(row)
        if int(row[10])>1:write_child_rows(n)
        if int(row[10])>2:write_child_rows(2*n-3)
        if int(row[10])>3:write_child_rows(3*n-6)
        if int(row[10])>4:write_child_rows(4*n-9)
    for row in rows:
        parentID+=1
        for i in range(len(row)):            
            if ('%s_All_conserved_hits' %(jobID) in table or '%s_top_ranked_hits' %(jobID) in table  or '%s_QueryListResults' %(jobID) in table) and ('CRM' not in table):                                
                mkTable(row,i)                                            
            elif ('%s_CRM_results' %(jobID) in table or 'QueryListResults_CRM' in table or 'top_ranked_per_gene_CRM' in table) and (matrix_count>1):
                mkTable2(row,i)                 
            else:
                if i==0:outfile.write('<tr class=\"PARENT row1\">')
                outfile.write('<td>%s</td>' %(row[i]))                        
        outfile.write('</tr>\n')
    outfile.write('</table>\n</body></html>')
    if db!='':cursor.execute("detach c")        

def mkViewAllPlots(matrix_dir):
    outfile=file('%s/%s/Results_plots/View_All_plots.html' %(output_dir,matrix_dir),'w')
    plotFiles=os.listdir('%s/%s/Results_plots/' %(output_dir,matrix_dir))
    plotFiles.sort()
    for plot in plotFiles:
        outfile.write('<img src=\"%s\" />' %(plot))

def mktableList():
    
    for species in speciesList:
        for m in range(matrix_count):
            db='%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,species,m+1)
            table='%s_%s_hit_ortho%s'%(jobID,species,m+1)
            mkHTMLtables(table,m+1,db)
            mvTablesFromSqltoTEXT(table,m+1,db)
    for n in range(5):
        for m in range(matrix_count):
            db='%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,m+1)
            table='%s_%s_conserved_hits%s' %(jobID,n+1,m+1)
            mkHTMLtables(table,m+1,db)
            mvTablesFromSqltoTEXT(table,m+1,db)
    
    for n in range(matrix_count):
        db='%s/run/sqlite_tmp_files/%s%s_TargetOrtho.db' %(TargetOrtho_path,jobID,n+1)
        
        table='%s_ResultsSummary%s' %(jobID,n+1)
        mkHTMLtables(table,n+1,db)
        mvTablesFromSqltoTEXT(table,'',db)
        table='%s_top_ranked_hits%s' %(jobID,n+1)
        mkHTMLtables(table,n+1,db)
        mvTablesFromSqltoTEXT(table,n+1,db)
        
        table='%s_All_conserved_hits%s_ranked' %(jobID,n+1)
        mkHTMLtables(table,n+1,db)
        
        mvTablesFromSqltoTEXT(table,n+1,db)
       
        if queryFile_path != 'None' and queryFile_path!=None:
            table='%s_QueryListResults%s' %(jobID,n+1)
            mkHTMLtables(table,n+1,db)
            mvTablesFromSqltoTEXT(table,n+1,db)
            if n==1:
                table='%s_QueryListResults_CRM' %(jobID)
                mkHTMLtables(table,'',db)
                mvTablesFromSqltoTEXT(table,'',db)
        if n==1:
            table='%s_CRM_results' %(jobID)
            mkHTMLtables(table,'',db)
            mvTablesFromSqltoTEXT(table,'',db)
            table='%s_top_ranked_per_gene_CRM' %(jobID)
            mkHTMLtables(table,'',db)
            mvTablesFromSqltoTEXT(table,'',db)
       
def arrange_hier():
    os.system("cp %s/run/scripts/jquery-1.8.2.min.js %s/.jquery-1.8.2.min.js" %(TargetOrtho_path,output_dir))
    os.system("cp %s/run/scripts/sorttable.js %s/.sorttable.js" %(TargetOrtho_path,output_dir))
    os.system("cp %s/run/images/download_icon.png %s/" %(TargetOrtho_path,output_dir))
    for n in range(matrix_count):           
        os.system("cp %s/run/scripts/sorttable.js %s/%s/Results_html/.sorttable.js" %(TargetOrtho_path,output_dir,matrixNamesList[n]))
        os.system("cp %s/run/scripts/jquery-1.8.2.min.js %s/%s/Results_html/.jquery-1.8.2.min.js" %(TargetOrtho_path,output_dir,matrixNamesList[n]))
        if matrix_count==1:pref=''
        else:pref='.'
        os.system("cp %s/%s/Results_html/%s_top_ranked_hits%s.html %s/%s%s_top_ranked_hits%s.html" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
        os.system("cp %s/%s/Results_text/%s_top_ranked_hits%s.txt %s/%s%s_top_ranked_hits%s.txt" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
        os.system("cp %s/%s/Results_html/%s_ResultsSummary%s.html %s/%s%s_ResultsSummary%s.html" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
        if queryFile_path !='None' or queryFile_path!=None:
            os.system("cp %s/%s/Results_html/%s_QueryListResults%s.html %s/%s%s_QueryListResults%s.html" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
            os.system("cp %s/%s/Results_text/%s_QueryListResults%s.txt %s/%s%s_QueryListResults%s.txt" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
            if QueryOnly!='True':
                os.system("cp %s/%s/Results_html/%s_All_conserved_hits%s_ranked.html %s/%s%s_All_conserved_hits%s_ranked.html" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
                os.system("cp %s/%s/Results_text/%s_All_conserved_hits%s_ranked.txt %s/%s%s_All_conserved_hits%s_ranked.txt" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
        else:
            os.system("cp %s/%s/Results_html/%s_All_conserved_hits%s_ranked.html %s/%s%s_All_conserved_hits%s_ranked.html" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))
            os.system("cp %s/%s/Results_text/%s_All_conserved_hits%s_ranked.txt %s/%s%s_All_conserved_hits%s_ranked.txt" %(output_dir,matrixNamesList[n],jobID,n+1,output_dir,pref,jobID,n+1))

        os.system("cp %s/%s/Results_text/%s_%s_%s.bed %s/%s%s_%s_%s.bed" %(output_dir,matrixNamesList[n],jobID,speciesList[0],matrixNamesList[n],output_dir,pref,jobID,speciesList[0],matrixNamesList[n]))    
    if matrix_count>1:
        try:os.makedirs("%s/All_tables"%(output_dir))
        except:'already exists'
        names=os.listdir(output_dir)
        for name in names:
            if name in matrixNamesList:
                os.system("mv %s/%s %s/All_tables/" %(output_dir,name,output_dir))

def mk_genome_browse_track_file(motifnum,species):
    outfile=file('%s/%s/Results_text/%s_%s_%s.bed' %(output_dir,matrixNamesList[motifnum-1],jobID,species,matrixNamesList[motifnum-1]),'w')
    outfile.write("track name= \'%s\' description=\"TargetOrtho Results\" useScore=1\n" %(matrixNamesList[motifnum-1]))
    cursor.execute("attach '%s/run/sqlite_tmp_files/%s_%s%s_TargetOrtho.db' as %s_%s%s_db" %(TargetOrtho_path,jobID,species,motifnum,jobID,species,motifnum))
    cursor.execute("select score from %s_%s_hit_ortho%s order by score desc limit 1" %(jobID,species,motifnum))
    try:max_score=float(cursor.fetchone()[0])
    except:max_score=10
    
    cursor.execute("select score from %s_%s_hit_ortho%s order by score limit 1" %(jobID,species,motifnum))
    try:min_score=float(cursor.fetchone()[0])
    except:min_score=-8
    
    rows1=[row for row in cursor.execute("select distinct dna, start, end, '%s_binding_site',((score)-(%s))/((%s)-(%s))*(1000-500)+500, '+' from %s_%s_hit_ortho%s where strand = '+'" %(matrixNamesList[motifnum-1],min_score,max_score,min_score,jobID,species,motifnum))]    
    rows2=[row for row in cursor.execute("select distinct dna, start, end, '%s_binding_site', ((score)-(%s))/((%s)-(%s))*(1000-500)+500, '-' from %s_%s_hit_ortho%s where strand = '-'" %(matrixNamesList[motifnum-1],min_score,max_score,min_score,jobID,species,motifnum))]    
    [([outfile.write('%s\t' %str(i).strip()) for i in row],outfile.write('\n')) for row in (rows1+rows2)]
    cursor.execute("detach '%s_%s%s_db'"%(jobID,species,motifnum)) 

def main():    
    tableList=mktableList()
    
    for n in range(matrix_count):
        motifnum = n+1
        mkViewAllPlots(matrixNamesList[n])
        for species in speciesList:
            mk_genome_browse_track_file(motifnum,species)    
    arrange_hier()
    galaxy_results_file(matrix_count)
    os.system("tar -zcvf %s/run/output/Download_%s.tar.gz %s/run/output/%s" %(TargetOrtho_path,jobID,TargetOrtho_path,jobID))
    os.system("mv %s/run/output/Download_%s.tar.gz %s/" %(TargetOrtho_path,jobID,output_dir))
    
    conn.close()
main()
