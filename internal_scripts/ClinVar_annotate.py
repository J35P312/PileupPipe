import sys
import argparse
import subprocess
import os
import fnmatch
import sqlite3

severity_2_text={"0":"Uncertain significance","1":"not provided,","2":"Benign","3":"Likely benign","4":"Likely pathogenic","5":"Pathogenic","6":"drug response","7":"histocompatibility","255":"other"}
def main(args):

    conn = sqlite3.connect(':memory:')
    c = conn.cursor()
    A="CREATE TABLE SVDB (chr TEXT, pos INT,alt TEXT, severity TEXT)"
    c.execute(A)

    exac_db=[]
    for line in open(args.db):
        if line[0] == "#":
            continue
        content=line.strip().split()
        severity=line.strip().split("\t")[7].split(";CLNSIG=")[-1].split(";")[0]
        if "|" in severity and not "," in severity:
            severity= max(severity.split("|"))

        if severity in severity_2_text:
            severity=severity_2_text[severity]
        else:
            severity="none"


        if not "," in content[4]:
            exac_db.append([ content[0], content[1] , content[4] ,severity ])
        else:
            alt_list=content[3].split(",")
            severity=line.strip().split("\t")[7].split(";CLNSIG=")[-1].split(";")[0].split(",")
            for i in range(0,len(alt_list)):

                if "|" in severity[i]:
                    severity[i]= max(severity[i].split("|"))
                if severity[i] in severity_2_text:
                    severity[i]=severity_2_text[severity[i]]
                else:
                    severity[i]="none"
                exac_db.append([ content[0], content[1] , alt_list[i] , severity[i] ])
                
        if len(exac_db) > 100000:
            c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',exac_db)          
            exac_db=[]
    if exac_db:
        c.executemany('INSERT INTO SVDB VALUES (?,?,?,?)',exac_db)
                    
    del exac_db
    
    A="CREATE INDEX SNP ON SVDB (chr, pos, alt)"
    c.execute(A)
    conn.commit()
    
    for line in open(args.vcf):
        if line[0] == "#" and line[1] == "#":
            print line.strip()
        elif line[0] == "#":
            print "##INFO=<ID={},Number=1,Type=String,Description=\"Clinical signifiance(CLINVAR).\">".format(args.tag)
            print line.strip()
        else:
            FRQ="none"
            annotation=";{}={}"
            content=line.strip().split("\t")

            chromosome=content[0]
            position=content[1]
            end=position
            alt=content[4]
            #cadd annotation
            
            A='SELECT severity FROM SVDB WHERE chr == \'{}\' AND pos == {} AND alt == \'{}\' '.format(chromosome,position,alt)
            
            hits=[]
            for hit in c.execute(A):
                FRQ = hit[0]
            content[7] += annotation.format(args.tag,FRQ)
            print("\t".join(content).strip())
	   		#popfreq
        
    conn.close()
parser = argparse.ArgumentParser("""this scripts annotates a SNP using CADD and popfreq, the popfreq and CADD file must be tabix indexed and tabbix must be installed""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--folder',type=str,help="used instead of vcf to annotate each vcf in a folder")
parser.add_argument('--db',type=str,help="the path to the exac DB")
parser.add_argument('--tag',type=str,default="CLINVAR",help="the vcf frequency tag(default = CLINVAR)")
args, unknown = parser.parse_known_args()

if args.vcf:
    main(args)
elif args.folder:
    for root, dirnames, filenames in os.walk(args.folder):
            for filename in fnmatch.filter(filenames, '*.vcf'):
                bam_file=os.path.join(root, filename)
                args.vcf=bam_file
                main(args)
else:
    print("|>L3453 5|>3(1/=Y --vcf || --folder")
