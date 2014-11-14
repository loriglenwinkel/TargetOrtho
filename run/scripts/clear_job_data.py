#!/usr/bin/python
#Author: Lori Glenwinkel
from optparse import OptionParser
import os
import sqlite3 as lite
parser = OptionParser()
parser.add_option('-j','--jobID',dest='jobID',help='-j 123345')
parser.add_option('-c','--TargetOrtho_path',dest='TargetOrtho_path',help='-c /data/newTargetOrtho')
( options, args ) = parser.parse_args()

TargetOrtho_path=options.TargetOrtho_path
jobID=options.jobID

def clear_tables(jobID):
    cursor.execute("SELECT name FROM sqlite_master WHERE type = 'table'")
    tables=cursor.fetchall()
    for table in tables:
        if table[0].startswith("%s" %(jobID)):cursor.execute("drop table %s" %(table[0]))

def clear_input_files(TargetOrtho_path,jobID):
    inputFiles=os.listdir("%s/run/input" %(TargetOrtho_path))
    for name in inputFiles:
        if name.startswith("%s" %(jobID)):
            os.system("rm -f %s/run/input/%s" %(TargetOrtho_path,name))

def clear_output_dir(TargetOrtho_path,jobID):
    try:os.system("rm -r %s/run/output/%s" %(TargetOrtho_path,jobID))
    except:print 'no output dir of that jobID'

def clear_fimo_out(TargetOrtho_path,jobID):
    inputFiles=os.listdir("%s/run/fimo_out" %(TargetOrtho_path))
    for name in inputFiles:
        if name.startswith("%s" %(jobID)):
            os.system("rm -f %s/run/fimo_out/%s" %(TargetOrtho_path,name))

def clear_ag_out(TargetOrtho_path,jobID):
    inputFiles=os.listdir("%s/run/associate_genes_out" %(TargetOrtho_path))
    for name in inputFiles:
        if name.startswith("%s" %(jobID)):
            os.system("rm -f %s/run/associate_genes_out/%s" %(TargetOrtho_path,name))
            
def clear_mysqld_db_dir():
    files=os.listdir("%s/run/sqlite_tmp_files" %TargetOrtho_path)
    for name in files:
        if name.startswith('%s' %(jobID)):
            os.system("rm %s/run/sqlite_tmp_files/%s" %(TargetOrtho_path,name))

def main():
    print 'clearing data for jobID: ', jobID
    clear_ag_out(TargetOrtho_path,jobID)
    clear_input_files(TargetOrtho_path,jobID)
    #clear_output_dir(TargetOrtho_path,jobID)
    clear_fimo_out(TargetOrtho_path,jobID)
    clear_mysqld_db_dir()
main()
