import sys
import argparse
from operator import itemgetter, attrgetter
import xlwt

def EFF(effects,pos):
    snp_dictionary={}
    for effect in effects:
        effect=effect.split(";")[-1]
        if "HIGH" in effect or "MODERATE" in effect:
            #get all information from the csq field of the present transcript
            gene=effect.split("|")[3]
            feature=effect.split("|")[11].split(":")[-1]
            cdna=effect.split("|")[10].split(":")[-1]
            consequence=effect.split("|")[1]
            sift=effect.split("|")[-3]
            phen=effect.split("|")[-2]
            high=0
            if "|HIGH|" in effect:
                high=1

            if not gene in snp_dictionary:
                snp_dictionary[gene]={}
            
            snp_dictionary[gene][consequence]=[feature,cdna,sift,phen,high]
        #if sequence_feature is not the only entry of a gene, 
    return(snp_dictionary)

def compute_rankscore(variant_dictionary):
    rank_score=0
    #add cadd score
    if variant_dictionary["cadd"] > 30:
        rank_score+=2
    elif variant_dictionary["cadd"] > 20:
        rank_score +=1

    #add points if the impact of the variant is high
    if variant_dictionary["high"]:
        rank_score +=2

    #add points for compound and homozygous variants
    if variant_dictionary["zygosity"] == "hom":
        rank_score += 2
    elif not variant_dictionary["zygosity"] == "het":
        rank_score += 1

    #increase the score of known pathogenic or likely pathogenic variants, decrease the score of benign variants
    if variant_dictionary["clinvar"] == "Pathogenic":
        rank_score += 5
    elif variant_dictionary["clinvar"] == "Likely pathogenic":
        rank_score += 5
    elif variant_dictionary["clinvar"] == "benign":    
        rank_score += -4

    #variants not located in in a blacklisted region gets a higher score 
    if not variant_dictionary["black_listed"]:
        rank_score +=1

    #increase the score if the maximum population frequency is below 0.01
    if max([ variant_dictionary["1kgaf"],variant_dictionary["sweaf"],variant_dictionary["exac"],variant_dictionary["gnomad"] ]) < 0.01:
        rank_score += 3

    return(rank_score)

    
parser = argparse.ArgumentParser("""turns a snpeff vcf into a csv file, output is printed to the stdout""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--frequency',type=float,default=0.01,help="frequency cut-off")
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
        qual=content[5]
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
            try:
                if float(cadd) < 10:
                   tolerated_cadd=1
                   cadd=float(cadd)
            except:
               cadd=""

        popfreq="0.0"
        txt=content[7].split(";1000GAF=")
        if len(txt) == 2:
            txt=txt[-1]
            popfreq=txt.split(";")[0]
        sweaf="0.0"
        txt=content[7].split(";SWEAF=")
        if len(txt) == 2:
            txt=txt[-1]
            sweaf=txt.split(";")[0]
        exac="0.0"
        txt=content[7].split(";EXACAF=")
        if len(txt) == 2:
            txt=txt[-1]
            exac=txt.split(";")[0]
        clinvar=""
        txt=content[7].split(";CLINVAR=")

        gnomad="0.0"
        txt=content[7].split(";GNOMADAF=")
        if len(txt) == 2:
            txt=txt[-1]
            gnomad=txt.split(";")[0]

        blacklist=0
        txt=content[7].split(";BLACKLIST=")
        if len(txt) == 2:
            txt=txt[-1]
            blacklist=int(txt.split(";")[0])

        if len(txt) == 2:
            txt=txt[-1]
            clinvar=txt.split(";")[0]        

        #have a look in the snpeff field
        snp_dictionary=EFF(line.strip().split("\t")[7].split("CSQ=")[-1].split(";")[0].split(","),pos )
        for gene in snp_dictionary:
            for variant in snp_dictionary[gene]:
                tolerated=0+tolerated_cadd
                feature=snp_dictionary[gene][variant][0]
                cdna=snp_dictionary[gene][variant][1]
                high=snp_dictionary[gene][variant][-1]
                poly=snp_dictionary[gene][variant][-2]

                sift=snp_dictionary[gene][variant][-3]

                if float(popfreq) > args.frequency or float(sweaf) > args.frequency or float(exac) > args.frequency or float(gnomad) > args.frequency:
                     continue
                rankscore=compute_rankscore({"zygosity":zygosity,"clinvar":clinvar,"poly":poly,"sift":sift,"cadd":cadd,"1kgaf":popfreq,"sweaf":sweaf,"exac":exac,"gnomad":gnomad,"black_listed":blacklist,"high":high})
                variant_list.append([chrom,pos,id_,ref, alt,feature,cdna,variant,gene,zygosity,clinvar,rankscore,blacklist,poly,sift,cadd,popfreq,sweaf,exac,gnomad])

filename=args.vcf.replace(".vcf",".xls")

wb =  xlwt.Workbook()
ws0 = wb.add_sheet("variants",cell_overwrite_ok=True)
i=0;
header=["Chromosome","Position","ID","Ref","Alt","AA-change","CDNA","Consequence","Gene","zygosity","ClinVar","RankScore","blackListed","PolyPhen","SIFT","CADD","thousand_genome","sweFreq","exac","gnomad"]
j=0
for item in header:
    ws0.write(i, j, item)
    j+=1
i=1
for entry in sorted(variant_list, key=itemgetter(11), reverse = True): 
    j=0;
    for item in entry:
        ws0.write(i, j, item)
        j+=1
    i+=1

wb.save(filename)
