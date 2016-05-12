#######################################################################################################################
#                                                                                                                     #
# HERMES is a straightforward index which tries to summarize the mitochondrial evolution pace using a single number.  #
# Several mitogenomic features are evaluated in a factor analysis framework; namely, in the current version, they are #
# %URs, AMIGA, SU skew, root-to-tip distance and pairwise ML distance.                                                #
#                                                                                                                     #
# Copyright (C) 2016 Guglielmo Puccio, Federico Plazzi                                                                #
#                                                                                                                     #
# This program is free software: you can redistribute it and/or modify                                                #
# it under the terms of the GNU General Public License as published by                                                #
# the Free Software Foundation, either version 3 of the License, or                                                   #
# (at your option) any later version.                                                                                 #
#                                                                                                                     #
# This program is distributed in the hope that it will be useful,                                                     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                       #
# GNU General Public License for more details.                                                                        #
#                                                                                                                     #
# You should have received a copy of the GNU General Public License                                                   #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                               #
#                                                                                                                     #
#######################################################################################################################

import subprocess
from Bio import SeqIO
from ete2 import Tree
import os, argparse

#arguments parsing
parser=argparse.ArgumentParser()
#change the optional title
parser._optionals.title = "Arguments"
#input file
parser.add_argument('-I',dest='gb_file',required=True,help='Input File containing the records in Genebank format')
#list of tree names and NCBI ids
parser.add_argument('-D',dest='entry_names',required=True,help='File listing all entry names with the corresponding NCBI ids')
#tree in Newick format
parser.add_argument('-L',dest='tree',required=True,help='Tree file in Newick format')
#outgroup in the tree
parser.add_argument('-O',dest='outgroup',required=True,help='Name of the outgroup')
#arguments for RAxML
parser.add_argument('-s',dest='alignment',required=True,help='Alignment file')
parser.add_argument('-m',dest='model',required=True,help='Molecular evolution model (for RAxML)')
parser.add_argument('-q',dest='partition',help='Partition file (for RAxML)')
parser.add_argument('-t',dest='threads',help='Number of CPU cores to be used (for RAxML multithreading)',default="1")
#alpha for RMSEA
parser.add_argument('-a',dest='alpha',type=float,help='Confidence level for RMSEA CI (for factor analysis)',default="0.05")
args = parser.parse_args()


output="HERMES"
version="1.0"
try:
    os.mkdir("./Results")
except:
    pass
file_for_R=open("./Results/HERMES_variables.txt","w")

try:
	os.remove("RAxML_distances."+output)
except:
	pass
try:
	os.remove("RAxML_info."+output)
except:
	pass
try:
	os.remove("RAxML_parsimonyTree."+output)
except:
	pass

def append_gene(name,feature):  #it takes the NCBI name and the feature from genbank
    Q = feature.qualifiers
    global f,D
    if 'gene' in Q.keys():
        if feature.strand == 1:
            f[name].append(["+",D[Q['gene'][0].lower()]])
        else:
            f[name].append(["-",D[Q['gene'][0].lower()]])
    else:
        try:
            if feature.strand == 1:
                f[name].append(["+",D[Q['product'][0].lower()]])
            else:
                f[name].append(["-",D[Q['product'][0].lower()]])
        except:     #if the feature isn't annotated as gene or product
            if feature.strand == 1:
                f[name].append(["+","orfan"])
            else:
                f[name].append(["-","orfan"])
    return


def change_sign(gene):
    if gene[0]== '-':
        gene[0] = '+'
    else:
        gene[0] = '-'
    return

#Here we create the general mitochondrial Dictionary from the file
DD = open("./Genes.dict").readlines()
D={}
for i in range(len(DD)):
    DD[i]=DD[i].strip().split(",")
    D[DD[i][0]]=DD[i][1].upper()

#Here we create the dictionary to translate the genbank ids with the names used in the tree (usually species)
NN=open(args.entry_names).readlines()
N={}
M={}
O={}
if len(NN[0].split(",")) > 2:
    for i in range(len(NN)):
        NN[i]=NN[i].strip().split(",")
        N[NN[i][0]]=NN[i][1]    #Dictionary with species as keys and NCBI IDs as values
        M[NN[i][1]]=NN[i][0]    #Dictionary with NCBI IDs as keys and species as values
        O[NN[i][1]]=int(NN[i][2])    #Dictionary with NCBI IDs as keys and color number as value
else:
    for i in range(len(NN)):
        NN[i]=NN[i].strip().split(",")
        N[NN[i][0]]=NN[i][1]    #Dictionary with species as keys and NCBI IDs as values
        M[NN[i][1]]=NN[i][0]    #Dictionary with NCBI IDs as keys and species as values

#-----main----#

f = {}
g = {}
h={}
for gb_record in SeqIO.parse(open(args.gb_file, "r"), "genbank"):
    F = gb_record.features
    LL=len(gb_record.seq)
    L = [1 for i in range(len(gb_record.seq))]
    f[(gb_record.id).split(".")[0]]=[]
    h[(gb_record.id).split(".")[0]]=[[],[]]
    for i in range(1,len(F)):
        if F[i].type in ["CDS","rRNA","tRNA"]:
            #for SU skew
            h[(gb_record.id).split(".")[0]][1].append(F[i].strand)
            #for URs
            for x in F[i].location.parts:
                for y in x:
                    if L[y]==1:
                        L[y]=0
            #for GeneOrder
        if F[i].type == "CDS":
            append_gene((gb_record.id).split(".")[0],F[i])
    L=sum(L)
    g[(gb_record.id).split(".")[0]]=(float(L)/LL)*100



first=[]

#to obtain the correct gene order
for key in f:
    for j in range(len(f[key])):
        if f[key][j][1]== 'CO1' and f[key][j][0]=='+':
            if [y[1] for y in f[key]].index('CO1') == 0:  #if cox1 is the first gene and it's on the +
                c=0
            else:   #if not, we move all the genes before cox1 to the end
                for k in range(0,[y[1] for y in f[key]].index('CO1')):
                    f[key].insert(len(f[key]), f[key].pop(0))
        elif f[key][j][1]== 'CO1' and f[key][j][0]=='-':
            if [y[1] for y in f[key]].index('CO1') == len(f[key]):  #if cox1 is the last gene on the - strand
                f[key].reverse()    #we just need to reverse the order
                for k in f[key]:
                    change_sign(k)  #and change the signs
            else:
                for k in reversed(range([y[1] for y in f[key]].index('CO1')+1,len(f[key]))):
                    f[key].insert(0, f[key].pop(len(f[key])-1))
                f[key].reverse()
                for k in f[key]:
                    change_sign(k)
    first.append("_".join([x[1] for x in f[key]]))


#to evaluate AMIGA score    
second=set(first)
third={}
fourth={}
for i in second:
    third[i]=0
for key in f:
    j="_".join([x[1] for x in f[key]])
    third[j] += 1
for key in f:
    j="_".join([x[1] for x in f[key]])
    fourth[key] = (third[j]-1)/float(len(f)-1)


#to evaluate the absolute value of SU skew
for key in h:
    m,p = 0,0
    for i in h[key][1]:
        if i == 1:
            p=p+1
        else:
            m=m+1
    h[key][0]= abs(float(p-m)/(p+m))


#to find Root-to-tip distance
t=Tree(args.tree,format=1)
T={}
for el in t.get_leaves():
    T[N[el.name]] = round(el.get_distance(t.name),3)    #by using N (the dictionary with the tree names) we obtain the NCBI id 


#ML_distance
W={}
W[N[args.outgroup]]="NA"   #the distance between the outgroup and itself is set here.
optionals=[args.threads,]
if args.partition:
    pa_bool=True
else:
    pa_bool=False
subprocess.call(["./raxml","-f","x","-s",args.alignment,"-n",output,"-m",args.model,"-q"*pa_bool,str(args.partition)*pa_bool,"-T",args.threads,"-p","123456"], stdout=open(os.devnull, "wb"))
distances=open("RAxML_distances."+output).readlines()
for i in range(len(distances)):
    distances[i]=distances[i].split()
    b=[distances[i][0],distances[i][1]]
    if args.outgroup in b:
        b.remove(args.outgroup)
        W[N[b[0]]] = distances[i][2]


#save the 5 fields in one file
if len(O) > 0:  #if len(O) > 0 means there is a third column in the species file, and we have to add the Shade column to HERMES_variables.txt file.
    file_for_R.write("ID\tURs\tAMIGA\tSUskew\tRtoTdist\tMLdist\tShades\n")
    for key in sorted(O , key=O.get, reverse=True):
        if key == N[args.outgroup]:
            continue
        else:
            if key in T.keys():
                file_for_R.write(M[key]+"\t"+str(g[key])+"\t"+str(fourth[key])+"\t"+str(h[key][0])+"\t"+str(T[key])+"\t"+str(W[key])+"\t"+str(O[key])+"\n")
            else:   #if the species isn't in the tree
                file_for_R.write(M[key]+"\t"+str(g[key])+"\t"+str(fourth[key])+"\t"+str(h[key][0])+"\t"+"NA"+"\t"+str(W[key])+"\t"+str(O[key])+"\n")
else:   #if O is empty, means that there isn't a third column in the species file and that there is no column O[key] to be added (shades).
    file_for_R.write("ID\tURs\tAMIGA\tSUskew\tRtoTdist\tMLdist\n")
    for key in f:
        if key == N[args.outgroup]:
            continue
        else:
            if key in T.keys():
                file_for_R.write(M[key]+"\t"+str(g[key])+"\t"+str(fourth[key])+"\t"+str(h[key][0])+"\t"+str(T[key])+"\t"+str(W[key])+"\n")
            else:
                file_for_R.write(M[key]+"\t"+str(g[key])+"\t"+str(fourth[key])+"\t"+str(h[key][0])+"\t"+"NA"+"\t"+str(W[key])+"\n")
file_for_R.close()

if type(args.alpha)==float:     #if alpha isn't specified in the command line, it is equal to 'None' and so it isn't a float.
    subprocess.call(["Rscript","--vanilla","HERMES-v"+version+".R",str(args.alpha)])
else:
    subprocess.call(["Rscript","--vanilla","HERMES-v"+version+".R"])

os.rename("./Rplots.pdf","./Results/HERMES.pdf")
os.remove("RAxML_distances."+output)
os.remove("RAxML_info."+output)
os.remove("RAxML_parsimonyTree."+output)
