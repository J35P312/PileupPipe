import sys
import argparse
from operator import itemgetter, attrgetter
import xlwt


def EFF(effects,order):
    snp_dictionary={}
    for effect in effects:
        effect=effect.split(";")[-1]
        effects=effect.split("|")
        gene=effects[order["SYMBOL"]]
        #get all information from the csq field of the present transcript
        feature=effects[order["HGVSp"]]
        cdna=effects[order["HGVSc"]]
        impact=effects[order["IMPACT"]]
        sift=effects[order["SIFT"]]
        phen=effects[order["PolyPhen"]]
        high=0
        if "HIGH" == effects[order["IMPACT"]]:
            high=1
        cadd=""
        if "cadd_CADD" in order:
            cadd=effects[order["cadd_CADD"]].split("&")[0]
            if cadd != "":
                cadd = float(cadd)

        if not gene in snp_dictionary:
            snp_dictionary[gene]={}
        gnomAD_AF=effects[order["gnomAD_AF"]]
        max_AF=effects[order["MAX_AF"]]
        kg_AF=effects[order["AF"]]
	clin_sig=effects[order["CLIN_SIG"]]
        snp_dictionary[gene][impact]={"CADD":cadd,"CLIN_SIG":clin_sig,"feature":feature,"cdna":cdna,"sift":sift,"poly":phen,"high":high,"gnomAD_AF":gnomAD_AF,"max_af":max_AF,"kg_AF":kg_AF,"consequence":effects[order["Consequence"]],"rsid":effects[order["Existing_variation"]]}

        
        #if sequence_feature is not the only entry of a gene, 
    return(snp_dictionary)

def compute_rankscore(variant_dictionary):
    rank_score=0
    #add cadd score
    if variant_dictionary["cadd"] != "":
        if float(variant_dictionary["cadd"]) > 30:
            rank_score+=2
        elif float(variant_dictionary["cadd"]) > 20:
            rank_score +=1
        elif float(variant_dictionary["cadd"]) < 10:
            rank_score -2
  
    #add points if the impact of the variant is high
    if variant_dictionary["high"]:
        rank_score +=2

    #add points for compound and homozygous variants
    if variant_dictionary["zygosity"] == "hom":
        rank_score += 2
    elif not variant_dictionary["zygosity"] == "het":
        rank_score += 1

    #increase the score of known pathogenic or likely pathogenic variants, decrease the score of benign variants
    if "pathogenic" in variant_dictionary["clinvar"]:
        rank_score += 5
    elif variant_dictionary["clinvar"] == "benign":    
        rank_score += -4

    #variants not located in in a blacklisted region gets a higher score 
    if not variant_dictionary["black_listed"]:
        rank_score +=1
s
    if variant_dictionary["impact"] == "LOW" or "MODIFIER" == variant_dictionary["impact"]:
        rank_score+= -5

    #increase the score if the maximum population frequency is below 0.01
    if max([ float(variant_dictionary["1kgaf"]),float(variant_dictionary["max_af"]),float(variant_dictionary["gnomad"]) ]) < 0.01:
        rank_score += 3

    return(rank_score)

    
parser = argparse.ArgumentParser("""turns a snpeff vcf into a csv file, output is printed to the stdout""")
parser.add_argument('--vcf',type=str,required=True,help="the path to the vcf file")
parser.add_argument('--frequency',type=float,default=0.01,help="frequency cut-off")
args, unknown = parser.parse_known_args()
variant_list=[]

order={}
for line in open(args.vcf):
    if line[0] == "#":
        if "##INFO=<ID=CSQ,Number" in line:
           tmp=line.strip().split("Format: ")[-1].split("\">")[0].split("|")
           i=0
           for entry in tmp:
               order[entry]=i
               i+=1 
    elif not "#" == line[0]:
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
        qual=content[5]
        snp_dictionary=EFF(line.strip().split("\t")[7].split("CSQ=")[-1].split(";")[0].split(","),order )
        ref=content[3]
        alt=content[4]

        blacklist=0
        txt=content[7].split(";BLACKLIST=")
        if len(txt) == 2:
            txt=txt[-1]
            blacklist=int(txt.split(";")[0])    

        #have a look in the snpeff field
        for gene in snp_dictionary:
            for variant in snp_dictionary[gene]:
                if snp_dictionary[gene][variant]["consequence"] == "downstream_gene_variant" or "upstream_gene_variant" == snp_dictionary[gene][variant]["consequence"] or "intron_variant" in snp_dictionary[gene][variant]["consequence"] or "intergenic" in snp_dictionary[gene][variant]["consequence"]:
                    continue

                if "HIGH" in snp_dictionary[gene] and variant != "HIGH":
                    continue
                elif "MODERATE" in snp_dictionary[gene] and ( variant == "LOW" or variant == "MODIFIER" ):
                    continue

                gnomad=snp_dictionary[gene][variant]["gnomAD_AF"]
                if gnomad == "":
                    gnomad="0.0"
	

                popfreq=snp_dictionary[gene][variant]["kg_AF"]
                if popfreq == "":
                    popfreq="0.0"

                max_af=snp_dictionary[gene][variant]["max_af"]
                if max_af == "":
                    max_af="0.0"
                cadd=""
                if "CADD" in snp_dictionary[gene][variant]:
                    cadd=snp_dictionary[gene][variant]["CADD"]

                feature=snp_dictionary[gene][variant]["feature"]
                cdna=snp_dictionary[gene][variant]["cdna"]
                high=snp_dictionary[gene][variant]["high"]
                poly=snp_dictionary[gene][variant]["poly"]
                id_=snp_dictionary[gene][variant]["rsid"]
                clin_sig=snp_dictionary[gene][variant]["CLIN_SIG"]

                sift=snp_dictionary[gene][variant]["sift"]

                if float(max_af) > args.frequency:
                     continue
                rankscore=compute_rankscore({"zygosity":zygosity,"clinvar":clin_sig,"poly":poly,"sift":sift,"cadd":cadd,"1kgaf":popfreq,"gnomad":gnomad,"max_af":max_af,"black_listed":blacklist,"high":high,"impact":variant,"cadd":cadd})
                variant_list.append([chrom,pos,id_,ref, alt,feature,cdna,snp_dictionary[gene][variant]["consequence"],variant,gene,zygosity,clin_sig,rankscore,cadd,blacklist,poly,sift,popfreq,gnomad,max_af])

filename=args.vcf.replace(".vcf",".xls")

wb =  xlwt.Workbook()
ws0 = wb.add_sheet("variants",cell_overwrite_ok=True)
i=0;
header=["Chromosome","Position","ID","Ref","Alt","AA-change","CDNA","Consequence","severity","Gene","zygosity","ClinVar","RankScore","CADD","blackListed","PolyPhen","SIFT","thousand_genome","gnomad","max_af"]
j=0
for item in header:
    ws0.write(i, j, item)
    j+=1
i=1
for entry in sorted(variant_list, key=itemgetter(12), reverse = True): 
    j=0;
    for item in entry:
        ws0.write(i, j, item)
        j+=1
    i+=1

wb.save(filename)
