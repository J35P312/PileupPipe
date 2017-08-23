import sys
import os
#accepts the path to a vcf as first command, prints each line containing the word given as second argument
inputFile=sys.argv[1];
folder=inputFile.split("/");
folder=folder[:len(folder)-1]
name=inputFile.split("/");
name =name[-1];
name = name.split(".");
name = name[:len(name)-1]
name=".".join(name);

path="/".join(folder)


outputFile=os.path.join(path,name+".FILTERED.pathogenic.csv")
bins=[]
events=[];
Nevents=0;
first =1
variants=sys.argv[2:]
genes=[]
for line in open(sys.argv[2]):
    genes.append("\|"+line.strip() + "\|" );

gene_list="|".join(genes)
os.system("grep -E \"#|{}\" {}".format(gene_list, inputFile))
