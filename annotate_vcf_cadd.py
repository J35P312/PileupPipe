import sys
import argparse
import subprocess
import os
import fnmatch

def main(args):
    cadd_found=False
    for line in open(args.vcf):
            if line[0] == "#" and line[1] == "#":
                print line.strip()
                if "##INFO=<ID=CADD," in line:
                    cadd_found=True
            elif line[0] == "#":
                if not cadd_found:
                    print "##INFO=<ID=CADD,Number=1,Type=Integer,Description=\"The CADD relative score for this alternative.\">"
                print line.strip()
            else:
                CADD=False
                annotation=";CADD={}"
                content=line.strip().split("\t")

                chromosome=content[0]
                position=int(content[1])
                end=position+1
                alt=content[4]
                #cadd annotation
                if args.cadd:
                    command=["tabix {} {}:{}-{}".format(args.cadd,chromosome,position,end)]
                    tmp=subprocess.check_output(command, shell = True);
                    output=tmp.split("\n")
                    for entry in output:
                        db_content=entry.split("\t")
                        if len(db_content) > 3:
                            if db_content[3] == alt and str(position) == db_content[1]:
                                CADD=db_content[5]
                                break
                if CADD:
                    content[7]+=annotation.format(CADD)
                print("\t".join(content).strip())

parser = argparse.ArgumentParser("""this scripts annotates a SNP using CADD, the CADD file must be tabix indexed and tabbix must be installed""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--folder',type=str,help="used instead of vcf to annotate each vcf in a folder")
parser.add_argument('--cadd',type=str,help="The path to the CADD DB")
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
