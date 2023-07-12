#!/usr/bin/env python

import sys
import os.path
import argparse
import pandas as pd
import re

def main(infilename, gene_db_file, cancer_gene_file, out_dir):
    

    #print( "\n")

    #if buffer_size < 0:
    #    buffer_size = 0
    #gene_db_file_format = os.path.splitext(gene_db_file)[1].lower()
    if not out_dir:
        out_dir = os.path.dirname(infilename)

    filename_short = os.path.basename(infilename)
    basename, suffix = os.path.splitext(filename_short)
    outfilename = os.path.join(out_dir, basename + '_annotated' + suffix)
    outbedpename = os.path.join(out_dir, basename + '_proximity' + '.bedpe')
    #print('{}, {}'.format(outfilename, outbedpename))
    #sys.exit(0)

# Input gene database BED file:
#chr1    11869   14409   DDX11L1
#chr1    14404   29570   WASH7P
#chr1    17369   17436   MIR6859-1
#chr1    29554   31109   MIR1302-2HG
#chr1    30366   30503   MIR1302-2
#chr1    34554   36081   FAM138A

# gtf format
#chr1	ENSEMBL	gene	17369	17436	.	-	.	gene_id "ENSG00000278267.1"; gene_type "miRNA"; gene_name "MIR6859-1"; level 3;
#chr1	HAVANA	gene	29554	31109	.	+	.	gene_id "ENSG00000243485.5"; gene_type "lincRNA"; gene_name "MIR1302-2HG"; level 2; tag "ncRNA_host"; havana_gene "OTTHUMG00000000959.2";
#chr1	ENSEMBL	gene	30366	30503	.	+	.	gene_id "ENSG00000284332.1"; gene_type "miRNA"; gene_name "MIR1302-2"; level 3;

# gff3 format
#chr1	HAVANA	gene	14404	29570	.	-	.	ID=ENSG00000227232.5;gene_id=ENSG00000227232.5;gene_type=unprocessed_pseudogene;gene_name=WASH7P;level=2;havana_gene=OTTHUMG00000000958.1
#chr1	ENSEMBL	gene	17369	17436	.	-	.	ID=ENSG00000278267.1;gene_id=ENSG00000278267.1;gene_type=miRNA;gene_name=MIR6859-1;level=3
#chr1	HAVANA	gene	29554	31109	.	+	.	ID=ENSG00000243485.5;gene_id=ENSG00000243485.5;gene_type=lincRNA;gene_name=MIR1302-2HG;level=2;tag=ncRNA_host;havana_gene=OTTHUMG00000000959.2

    dict_gene_db = {}
    with open(gene_db_file, 'r') as GENE_DB:
        for line in GENE_DB:
            if line.startswith('#'):
                continue
            x = line.split('\t')

            # Start and End positions of the gene, with sequence numbering starting at 1
            chr, start, end, gene_name = (x[0], int(x[1]), int(x[2]), x[3].strip())
            if not chr.lower().startswith('chr'):
                chr = "chr" + chr

            if chr in dict_gene_db:
                dict_gene_db[chr].add((start, end, gene_name))
            else:
                dict_gene_db[chr] = {(start, end, gene_name)}

            #print(dict_gene_db[chr])
            #print("{}\t{}\t{}\t{}".format(chr, start, end, gene_name))
            #sys.exit(0)

    dict_cancer_genes = {}
    with open(cancer_gene_file, 'r') as CANCER_GENES:
        for line in CANCER_GENES:
            if line.startswith('#'):
                continue
            x = line.split('\t')

            # Start and End positions of the gene, with sequence numbering starting at 1
            chr, start, end, gene_name = (x[0], int(x[1]), int(x[2]), x[3].strip())
            if not chr.lower().startswith('chr'):
                chr = "chr" + chr

            if chr in dict_cancer_genes:
                dict_cancer_genes[chr].add((start, end, gene_name))
            else:
                dict_cancer_genes[chr] = {(start, end, gene_name)}

            #print(dict_cancer_genes[chr])
            #print("{}\t{}\t{}\t{}".format(chr, start, end, gene_name))
            #sys.exit(0)

#chr1	x1	x2	chr2	y1	y2	strand1	strand2	resolution	=-logP	Sample
# chr19	36224000	36249000	chr1	27667000	27680000	-	-	1kb	175.962	RT4_rep1_S1
# chrX	5120000 	5176000	    chrX	130137000	130213000	-	+	1kb	299.344	WSU-DLCL2_rep2_S4
# chr13	46957000	47006000	chr13	50520000	50536000	-	+	1kb	77.8997	WSU-DLCL2_rep2_S4

    buffer = 1000000

    with open(infilename, 'r') as INFILE, open(outfilename, 'w') as OUTFILE, open(outbedpename, 'w') as OUTBEDPE:
        for line in INFILE:
            line = line.strip()
            if line.startswith('#'):
                OUTBEDPE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("#chr1", "Prox_start1", "Prox_end1", "chr2", "Prox_start2", "Prox_end2", "strand1", "strand2"))

                OUTFILE.write(line + "\n")
                continue
            x = line.split('\t')


# The strand predictions are mean to be estimates of whether the rearrangement breakpoint is at a given edge of the sub-matrix.
# A "+" value means that we predict the "end" coordinate to be closests to the actually breakpoint.
# A "-" value indicates we believe the "start" coordinate is closest to the breakpoint.

            chr_1, start_1, end_1, strand_1 = (x[0], int(x[1]), int(x[2]), x[6])
            chr_2, start_2, end_2, strand_2 = (x[3], int(x[4]), int(x[5]), x[7])

            if not chr_1.lower().startswith('chr'):
                chr_1 += "chr"
            if not chr_2.lower().startswith('chr'):
                chr_2 += "chr"

            # Define the "exact" position and the proximity region of each breakpoint
            # TODO: take care of the 0-based vs 1-based situation!
            if strand_1 == "-":
                proximity_start_1 = start_1 + 1
                proximity_end_1 = start_1 + buffer
                breakpoint_1 = start_1
            else:
                proximity_start_1 = max(0, end_1 - buffer)
                proximity_end_1 = max(0, end_1 - 1)
                breakpoint_1 = end_1

            if strand_2 == "-":
                proximity_start_2 = start_2 + 1
                proximity_end_2 = start_2 + buffer
                breakpoint_2 = start_2
            else:
                proximity_start_2 = max(0, end_2 - buffer)
                proximity_end_2 = max(0, end_2 - 1)
                breakpoint_2 = end_2

            breakpoint1_genes = set()
            breakpoint2_genes = set()
            breakpoint1_cancer_genes = set()
            breakpoint2_cancer_genes = set()
            breakpoint1_proximity_cancer_genes = set()
            breakpoint2_proximity_cancer_genes = set()

            # Check breakpoint1 against the cancer gene list:
            if chr_1 in dict_cancer_genes:
                for gene in dict_cancer_genes[chr_1]:
                    #print(dict_cancer_genes[chr_1])
                    #sys.exit(0)
                    if overlap(gene[0], gene[1], breakpoint_1, breakpoint_1):
                        breakpoint1_cancer_genes.add(gene[2])

            # Check breakpoint2 against the cancer gene list:
            if chr_2 in dict_cancer_genes:
                for gene in dict_cancer_genes[chr_2]:
                    #print(dict_cancer_genes[chr_2])
                    #sys.exit(0)
                    if overlap(gene[0], gene[1], breakpoint_2, breakpoint_2):
                        breakpoint2_cancer_genes.add(gene[2])

            # Check breakpoint1 against the entire gene database:
            if chr_1 in dict_gene_db:
                for gene in dict_gene_db[chr_1]:
                    #print(dict_gene_db[chr_1])
                    #sys.exit(0)
                    if overlap(gene[0], gene[1], breakpoint_1, breakpoint_1):
                        breakpoint1_genes.add(gene[2])

            # Check breakpoint2 against the entire gene database:
            if chr_2 in dict_gene_db:
                for gene in dict_gene_db[chr_2]:
                    #print(dict_gene_db[chr_2])
                    #sys.exit(0)
                    if overlap(gene[0], gene[1], breakpoint_2, breakpoint_2):
                        breakpoint2_genes.add(gene[2])

            # Calculate proximity distances
            breakpoint1_proximity_distance = []
            breakpoint2_proximity_distance = []

            if chr_1 in dict_cancer_genes:
                for gene in dict_cancer_genes[chr_1]:
                    if overlap(gene[0], gene[1], proximity_start_1, proximity_end_1):
                        breakpoint1_proximity_cancer_genes.add(gene[2])
                        if strand_1 == "-":
                            breakpoint1_proximity_distance.append(gene[0] - start_1)
                        else:
                            breakpoint1_proximity_distance.append(end_1 - gene[1])

            if chr_2 in dict_cancer_genes:
                for gene in dict_cancer_genes[chr_2]:
                    if overlap(gene[0], gene[1], proximity_start_2, proximity_end_2):
                        breakpoint2_proximity_cancer_genes.add(gene[2])
                        if strand_2 == "-":
                            breakpoint2_proximity_distance.append(gene[0] - start_2)
                        else:
                            breakpoint2_proximity_distance.append(end_2 - gene[1])

            breakpoint1_genes_str = ",".join(sorted(breakpoint1_genes))
            breakpoint2_genes_str = ",".join(sorted(breakpoint2_genes))
            breakpoint1_cancer_genes_str = ",".join(sorted(breakpoint1_cancer_genes))
            breakpoint2_cancer_genes_str = ",".join(sorted(breakpoint2_cancer_genes))
            breakpoint1_proximity_cancer_genes_str = ",".join(sorted(breakpoint1_proximity_cancer_genes))
            breakpoint2_proximity_cancer_genes_str = ",".join(sorted(breakpoint2_proximity_cancer_genes))

            # Assign gene fusion type

            # Class 1: Gene To Gene Fusion (Type 1)
            # Class 1A: B1: cancer gene, B2: cancer gene
            # Class 1B: B1: cancer gene, B2: human gene (B1 and B2 can be reversed)
            # Class 1C: B1: human gene, B2: human gene, B1P and B2P: make sure no cancer genes are proximal

            if len(breakpoint1_cancer_genes_str) > 0 and len(breakpoint2_cancer_genes_str) > 0:
                fusion_type = "1A"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, breakpoint1_cancer_genes_str, breakpoint2_cancer_genes_str))
            elif ( len(breakpoint1_cancer_genes_str) > 0 and len(breakpoint2_genes_str) > 0 or
                   len(breakpoint1_genes_str) > 0 and len(breakpoint2_cancer_genes_str) > 0 ):
                fusion_type = "1B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, breakpoint1_cancer_genes_str + "[" + breakpoint1_genes_str + "]", breakpoint2_cancer_genes_str + "[" + breakpoint2_genes_str + "]"))
            elif ( len(breakpoint1_genes_str) > 0 and len(breakpoint2_genes_str) > 0 and
                   len(breakpoint1_proximity_cancer_genes_str) == 0 and len(breakpoint2_proximity_cancer_genes_str) == 0 ):
                fusion_type = "1C"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, "[" + breakpoint1_genes_str + "]", "[" + breakpoint2_genes_str + "]"))

            # Class 2: Gene Fusion (Type 2)
            # Class 2: Break in cancer gene, partner intergenic OR Break in intergenic, partner genic break
            # Class 2A: B1: cancer gene, B2: intergenic (B1 and B2 can be reversed)
            # Class 2B: B1: human gene, B2: intergenic, B1P and B2P: make sure no cancer genes are proximal (B1 and B2 can be reversed)

            elif ( len(breakpoint1_cancer_genes_str) > 0 and len(breakpoint2_genes_str) == 0 or
                   len(breakpoint1_genes_str) == 0 and len(breakpoint2_cancer_genes_str) > 0 ):
                fusion_type = "2A"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, breakpoint1_cancer_genes_str, breakpoint2_cancer_genes_str))

            elif (( len(breakpoint1_genes_str) > 0 and len(breakpoint2_genes_str) == 0 or
                    len(breakpoint1_genes_str) == 0 and len(breakpoint2_genes_str) > 0 ) and
                    len(breakpoint1_proximity_cancer_genes_str) == 0 and len(breakpoint2_proximity_cancer_genes_str) == 0 ):
                fusion_type = "2B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, "[" + breakpoint1_genes_str + "]", "[" + breakpoint2_genes_str + "]"))

            # Proximity Fusions
            # Class 3A: Break in gene adjacent to cancer gene, partner genic break
            # B1: human gene, B2: human gene, B1P: cancer gene, B2P: cancer gene or not? Outout gene list if it exists
            # B2: human gene, B1: human gene, B2P: cancer gene, B1P: cancer gene or not? Outout gene list if it exists

            # Class 3B: Break in gene adjacent to cancer gene, partner intergenic break
            # B1: human gene, B2: intergenic, B1P: cancer gene, B2P: cancer gene or not? Outout gene list if it exists
            # B2: human gene, B1: intergenic, B2P: cancer gene, B1P: cancer gene or not? Outout gene list if it exists

            # Class 3B: Break in intergenic region adjacent to cancer gene, partner genic break
            # B1: intergenic, B2: human gene, B1P: cancer gene, B2P: cancer gene or not? Outout gene list if it exists
            # B2: intergenic, B1: human gene, B2P: cancer gene, B1P: cancer gene or not? Outout gene list if it exists

            # Class 3C: Break in intergenic region adjacent to cancer gene, partner intergenic break
            # B1: intergenic, B2: intergenic, B1P: cancer gene, B2P: cancer gene or not? Outout gene list if it exists
            # B2: intergenic, B1: intergenic, B2P: cancer gene, B1P: cancer gene or not? Outout gene list if it exists

            #elif len(breakpoint1_cancer_genes_str) == 0 and len(breakpoint2_cancer_genes_str) == 0:
            elif ( len(breakpoint1_genes_str) > 0 and len(breakpoint2_genes_str) > 0 and
                 ( len(breakpoint1_proximity_cancer_genes_str) > 0 or len(breakpoint2_proximity_cancer_genes_str) > 0 ) ):
                fusion_type = "3A"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str + "[" + breakpoint1_genes_str + "]", breakpoint2_proximity_cancer_genes_str + "[" + breakpoint2_genes_str + "]"))

            elif len(breakpoint1_genes_str) > 0 and len(breakpoint2_genes_str) == 0 and len(breakpoint1_proximity_cancer_genes_str) > 0:
                fusion_type = "3B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str + "[" + breakpoint1_genes_str + "]", breakpoint2_proximity_cancer_genes_str))
            elif len(breakpoint1_genes_str) == 0 and len(breakpoint2_genes_str) > 0 and len(breakpoint2_proximity_cancer_genes_str) > 0:
                fusion_type = "3B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str, breakpoint2_proximity_cancer_genes_str + "[" + breakpoint2_genes_str + "]"))

            elif len(breakpoint1_genes_str) == 0 and len(breakpoint2_genes_str) > 0 and len(breakpoint1_proximity_cancer_genes_str) > 0:
                fusion_type = "3B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str + "[" + breakpoint1_genes_str + "]", breakpoint2_proximity_cancer_genes_str + "[" + breakpoint2_genes_str + "]"))
            elif len(breakpoint1_genes_str) > 0 and len(breakpoint2_genes_str) == 0 and len(breakpoint2_proximity_cancer_genes_str) > 0:
                fusion_type = "3B"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str + "[" + breakpoint1_genes_str + "]", breakpoint2_proximity_cancer_genes_str + "[" + breakpoint2_genes_str + "]"))

            elif ( len(breakpoint1_genes_str) == 0 and len(breakpoint2_genes_str) == 0 and
                 ( len(breakpoint1_proximity_cancer_genes_str) > 0 or len(breakpoint2_proximity_cancer_genes_str) > 0 ) ):
                fusion_type = "3C"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, min(breakpoint1_proximity_distance, default=0), min(breakpoint2_proximity_distance, default=0), breakpoint1_proximity_cancer_genes_str, breakpoint2_proximity_cancer_genes_str))

            else:
                fusion_type = "0"
                OUTFILE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(line, fusion_type, proximity_start_1, proximity_end_1, proximity_start_2, proximity_end_2, 0, 0, "", ""))

            #else:
            #    print("ERROR: Wrong gene fusion type!\n" + line)
            #    sys.exit(1)

            # Output the proximity region to a BEDPE file
            OUTBEDPE.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(chr_1, proximity_start_1, proximity_end_1, chr_2, proximity_start_2, proximity_end_2, strand_1, strand_2))
    return outfilename


def overlap(left_1, right_1, left_2, right_2):
    if right_1 < left_2 or left_1 > right_2:
        return False
    else:
        return True

def get_cancer_gene_list(outfilename, path):
    '''a function that takes gene annotation bedpe file as input 
    and returns the cancer gene list in a text file'''
    df = pd.read_csv(f'{path}/{outfilename}', 
                     sep='\t',
                     skiprows=1, ### skipping first row because the headers are incomplete ####
                     names=['#chr1',	'x1',	'x2',	'chr2',	'y1',	'y2',	
                            'strand1',	'strand2',	'Fusion_type',	'Prox_start1',	
                            'Prox_end1',	'Prox_start2',	'Prox_end2',	'Dist1',	
                            'Dist2',	'Proximity_genes1',	'Proximity_genes2'])
    # remove the rows that do not have any proximity data
    df = df[df.Fusion_type != '0']
    # split the fusion type into fusion_class (1,2,3) and fusion_type (A,B,C)
    df['fusion_class'] = df['Fusion_type'].apply(lambda x: x[0])
    df['fusion_type']  = df['Fusion_type'].apply(lambda x: x[1])
    # only get rows with class 3
    df1 = df[df.fusion_class == '3']
    # replace "nan" with ''
    df1 = df1.fillna('')
    # extract proximity data from the dataframe as two lists
    prox_list1=list(df1.Proximity_genes1)
    prox_list2=list(df1.Proximity_genes2)
    ar1=list({i.strip() for i in set(prox_list1) if i != ''})
    ar2=list({i.strip() for i in set(prox_list2) if i != ''})
    # finally generate the gene list
    gene_list = []
    # parse data based on presence of "[" and report genes that are not in
    # square brackets
    for i in ar1:
        # if the data has "["
        if '[' in i: 
            parse_ar1 = i.split('[')
            if parse_ar1[0]:
                gene_list.append(parse_ar1[0])
                print (parse_ar1)
        # if the data doesn't have "["
        elif '[' not in i:
            gene_list.append(i)

    for i in ar2:
        if '[' in i: 
            parse_ar2 = i.split('[')
            if parse_ar2[0]:
                gene_list.append(parse_ar2[0])
        elif '[' not in i:
            gene_list.append(i)
    # remove duplicates by using set
    final_set = list(set(gene_list))

    # write the file
    with open(f'{path}/proximity_cancer_genelist.txt', 'w') as file:
        for i in final_set:
            file.write(f'{i}\n')

    

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Tool for assigning gene fusion type using bedpe file containing SVs.")

    parser.add_argument('-i', metavar='INFILE', type=str, dest='infilename', help='Input HiC-breakfinder BEDPE file including strand information', required=True)
    parser.add_argument('-g', metavar='GENE_DB_FILE', type=str, dest='gene_db_file', help='Gene database BED file', required=True)
    parser.add_argument('-c', metavar='CANCER_GENE_FILE', type=str, dest='cancer_gene_file', help='Cancer gene BED file', required=True)
    parser.add_argument('-o', metavar='OUT_DIR', type=str, dest='out_dir', help='Output directory [default: the same as the input file]')
    # parser.add_argument('-b', metavar='BUFFER_SIZE', type=int, dest='buffer_size', default=0, help='Buffer size for overlapping genes (bp) [default: 0]')

    args= parser.parse_args()
    
    infilename = args.infilename
    gene_db_file = args.gene_db_file
    cancer_gene_file = args.cancer_gene_file
    out_dir = args.out_dir
    
    outfilename = main(infilename, gene_db_file, cancer_gene_file, out_dir)
    path = os.getcwd()
    
    get_cancer_gene_list(outfilename,path)
    
