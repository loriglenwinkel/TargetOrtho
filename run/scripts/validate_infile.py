#!/usr/bin/python
#Author: Lori Glenwinkel
import os
import sys

def validate(infile,TargetOrtho_path,jobID):
    print 'validating motif input'
    try:
        lines=file(infile,'r').readlines()
        lineCount=0
        motif_names=[]
        motifCount=0
        foundCount=0
        nameCount=0
        foundCount2=0
        def write_meme_header(outfile,motifCount,bg):
            outfile.write('MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies (from uniform background):\n%s\n' %(bg))
            print 'writing outfile'
        for n in range(len(lines)):        
            if 'Background' in lines[n]:
                bg= lines[n+1]
            if lines[n].startswith('MOTIF'):
                name=lines[n].split()[1]
                motif_names.append(name)
                motifCount+=1
                nameCount=1
                outfile=file('%s/run/input/%s_meme%s.txt' %(TargetOrtho_path,jobID,motifCount),'w')
                lineCount=1
                write_meme_header(outfile,motifCount,bg)
            if lines[0].startswith('***'):
                if lines[n].startswith('log-odds'):
                #outfile.write(lines[n])
                    foundCount=1                        
                    foundCount2=1
                    if nameCount==1:outfile.write('\nMOTIF %s\n\n' %(name))                
                    else:outfile.write('\n')
                    nameCount+=1
                if lines[n].startswith('---'):                
                    foundCount=0
                if foundCount==1:
                    outfile.write(lines[n])        
                    
                    print lines[n]
            elif lineCount>0:
                outfile.write(lines[n])    
                print lines[n]
                if lines[n].startswith('log-odds'):
                    foundCount2=1
        if foundCount2!=1:
            print 'no log-odds matrix found'
            raise
    except:
        print "Unexpected error: problem with input file. Make sure to use PSSM format motif file (log-odds motif file) or use entire MEME text file as input)", sys.exc_info()[0]
        raise
        
    print motifCount,motif_names
    return motifCount,motif_names
                               
                              
	     

def main():
    validate(options.infile,options.TargetOrtho_path,options.jobID)
#main()
 
