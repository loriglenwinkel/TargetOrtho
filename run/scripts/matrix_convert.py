#Lori Glenwinkel May 31 2012
#convert letter-probability matrix to log-odds matrix. Prints to stdout
#takes a bg file order 0

import math
from optparse import OptionParser
import decimal

parser = OptionParser()
parser.add_option( '-b', '--bgfile', dest='bgfile', help='background file order 0 only' )
parser.add_option( '-f', '--file', dest='infile', help='lett-probability matrix input file')
parser.add_option( '-m', '--output_matrixType', dest='output_matrixType',help='-m log-odds, -m pfm')
parser.add_option('-i','--input_matrixType', dest='input_matrixType', help= '-i letter-prob')
( options, args ) = parser.parse_args()
psuedo_count= .000000001 #adjust later to use psuedo_count if 0 division error encountered or log2(0) encountered.

if options.bgfile!=None:
    bgfile=file('%s'%(options.bgfile),'r')
    bg=bgfile.readlines()[1:]
    bfs=[]
    for f in bg:
        bfs.append(float(f[2:-1]))
else:bfs=[.25,.25,.25,.25]

A=bfs[0]
C=bfs[1]
G=bfs[2]
T=bfs[3]

print 'background letter frequencies [A,C,G,T]: ', bfs

def lpmConvert(output_matrixType,e):
    if output_matrixType=='log-odds':
        print "%s\t%s\t%s\t%s" %(int(math.log(float(e[0]+ psuedo_count)/A,2)*100), int(math.log(float(e[1]+psuedo_count)/C,2)*100), int(math.log(float(e[2]+ psuedo_count)/G,2)*100), int(math.log(float(e[3]+ psuedo_count)/T,2)*100))
    elif output_matrixType=='pfm':
        m=abs(min([decimal.Decimal(str(e[i])).as_tuple().exponent for i in range(4)]))
        print "%s\t%s\t%s\t%s" %(int(e[0]*(10**m)),int(e[1]*(10**m)),int(e[2]*(10**m)),int(e[3]*(10**m)),)

def pfmConvert(output_matrixType,e):
    print 'converting'

def read(infile):
    f=file(infile,'r')
    matrix=f.readlines()
    for line in matrix:
        if 'name' in line:
            print line.strip()
    print options.output_matrixType, 'matrix: alength=4 width = n' #log-odds matrix: alength= 4'
    for line in matrix:
        try:
            e=[float(x) for x in line.split()]
            if options.input_matrixType=='letter-prob':lpmConvert(options.output_matrixType,e)
            elif options.input_matrixType=='pfm':pfmConvert(options.output_matrixType,e)
        except:None
def main():
    read(options.infile)
main()
