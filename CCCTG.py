import sys
import argparse
from operator import itemgetter, attrgetter
import xlwt

def EFF(effects,pos):
    snp_dictionary={}
    if pos == "27763405":
        print effects
    for effect in effects:
        effect=effect.split(";")[-1]
        if "HIGH" in effect or "MODERATE" in effect:
            gene=effect.split("|")[3]
            feature=effect.split("|")[15]
            cdna=effect.split("|")[12]
            consequence=effect.split("|")[1]
            sift=effect.split("|")[-2]
            phen=effect.split("|")[-1]
            #for each snp, each unique effect of each gene shoudl only be reported once, ie not multiple intron variant GENEX for snp z
            if not gene in snp_dictionary:
                snp_dictionary[gene]={}
            snp_dictionary[gene][consequence]=[feature,cdna,sift,phen]
        #if sequence_feature is not the only entry of a gene, 
    return(snp_dictionary)

    
parser = argparse.ArgumentParser("""turns a snpeff vcf into a csv file, output is printed to the stdout""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
args, unknown = parser.parse_known_args()
variant_list=[]

for line in open(args.vcf):
    if not "#" == line[0]:
        content=line.strip().split("\t")
        zygosity = "";
        #check the zygozity of the variant, using the GT field
        if len(content) >= 10:
            format_=content[8].split(":")
            if "GT" in format_:
                sample=content[9].split(":")
                pos=format_.index("GT")
                if( sample[pos] == "1/0" or sample[pos] == "0/1"):
                    zygosity="het"
                elif(sample[pos] == "1/1"):
                    zygosity="hom"
        #grab all information from the different fields
        chrom=content[0]
        pos=content[1]
        id_=content[2]
        if pos == "27763405":
            print line
        if id_ ==0:
            id= "."
        ref=content[3]
        alt=content[4]
        cadd=""
        txt=content[7].split(";CADD=")
        tolerated_cadd=0
        if len(txt) == 2:
            txt=txt[-1]
            cadd=txt.split(";")[0]
            if int(CADD) < 10:
                tolerated_cadd=1

        popfreq=""
        txt=content[7].split(";1000GAF=")
        if len(txt) == 2:
            txt=txt[-1]
            popfreq=txt.split(";")[0]
        sweaf=""
        txt=content[7].split(";SWEAF=")
        if len(txt) == 2:
            txt=txt[-1]
            sweaf=txt.split(";")[0]
        exac=""
        txt=content[7].split(";EXACAF=")
        if len(txt) == 2:
            txt=txt[-1]
            exac=txt.split(";")[0]
        clinvar=""
        txt=content[7].split(";CLINVAR=")
        if len(txt) == 2:
            txt=txt[-1]
            clinvar=txt.split(";")[0]        

        #have a look in the snpeff field
        snp_dictionary=EFF(line.strip().split("\t")[7].split("CSQ=")[-1].split(";")[0].split(","),pos )
        if pos == "27763405":
            print snp_dictionary
            print clinvar
        for gene in snp_dictionary:
            for variant in snp_dictionary[gene]:
                tolerated=0+tolerated_cadd
                feature=snp_dictionary[gene][variant][0]
                cdna=snp_dictionary[gene][variant][1]
                poly=snp_dictionary[gene][variant][-1]
                if pos == "27763405":
                    print tolerated
                if "benign" in poly:
                    tolerated+=1
                sift=snp_dictionary[gene][variant][-2]
                if "tolerated" in sift:
                    tolerated+=1
                if tolerated >= 2:
                    continue
                
                variant_list.append([chrom,pos,id_,ref, alt,feature,cdna,variant,gene,zygosity,clinvar,poly,sift,cadd,popfreq,sweaf,exac])

filename=args.vcf.replace(".vcf",".xls")

wb =  xlwt.Workbook()
ws0 = wb.add_sheet("variants",cell_overwrite_ok=True)
i=0;
header=["Chromosome","Position","ID","Ref","Alt","AA-change","cDNA","Consequencec","Gene","zygosity","ClinVar","PolyPhen","SIFT","CADD","thousand_genome","sweFreq","exac"]
j=0
for item in header:
    ws0.write(i, j, item)
    j+=1
i=1
for entry in sorted(variant_list, key=itemgetter(13), reverse = True): 
    j=0;
    for item in entry:
        ws0.write(i, j, item)
        j+=1
    i+=1

wb.save(filename)
