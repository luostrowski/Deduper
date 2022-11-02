#!/usr/bin/env python

'''
This is a script to eliminate PCR duplicates from a SAM file. Keep in mind that the input file must be sorted by
chromosome. This is a very conservative script, so we are probably removing more reads than expected.

'''

import re
import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program to eliminate PCR duplicates from a sorted SAM file.")
    parser.add_argument("-f", "--file", help="Absolute file path to sorted sam file", type=str, required=True)
    parser.add_argument("-o", "--output", help="Absolute file path to sorted output sam file", type=str, required=True)
    parser.add_argument("-u", "--umi", help="File containing the list of UMIs", type=str, required=True)
    return parser.parse_args()

args = get_args()
f = args.file
o = args.output
u = args.umi

##use a sorted file!
#test_file = "/projects/bgmp/luizah/bioinfo/Bi624/Deduper/test_sorted.sam"
#umi_file = "/projects/bgmp/luizah/bioinfo/Bi624/Deduper/STL96.txt"

def get_umi(umi:str) -> str:
    '''This function finds the UMI in the header of the SAM file
    and check if the UMI is in the known UMIs file.'''
    #umi = re.findall('[A,T,C,G]{8}', qname)
    umi = qname.split(":")[-1]
    return umi

def is_minus_strand(flag:int) -> bool:
    '''This function check the bitwise flag and returns if the read is
    in minus strand.'''
    return (flag & 16) == 16

def parse_cigar(cigar:str) -> list:
    '''
    This function parses the CIGAR string and returns a list with tuples of the 
    form (digit,letter).
    '''
    parsed = re.findall('(\d+)([A-Z]{1})', cigar)
    return parsed

def adjust_pos(old_pos:int, is_minus_strand:bool, parsed_cigar:list) -> int:
    i=0
    shift=0 
    if is_minus_strand:
        shift = -1 ##variable for how much we need to adjust in the position
    for length, op in parsed_cigar: ##separating tuple returned by parsed_cigar function in lenght and operation
        if is_minus_strand: 
            if op == "N" or op == "D" or op == "M":
                shift = shift + int(length) ##sum the length 
            if op == "S" and i > 0: ##soft-clipping not at the begining of the CIGAR
                shift = shift + int(length)
        else:
            if op == "S" and i == 0: ##soft-clipping at the beginning
                shift = shift - int(length)
        i+=1
    return old_pos + shift


umi_list=[]
read_dict={} # saving the non dup reads
curr_chrom = None

with open(u, "r") as fq: ##saving known UMIs in a list 
    for line in fq:
        umi_list.append(line.strip())

pcr_duplicates = open("PCR_duplicates.sam", "w") ##opening files to store PCR duplicates
wrong_umis = open("wrong_UMI.sam", "w") ##opening files to store reads with wrong UMIs

with open(f, "r") as fh, open(o, "w") as of:
    for line in fh:
        if line.startswith("@"): ##print headers directly into the output file
            print(line.strip("\n"), file=of)
        else: ##if line is not a header
            line = line.strip("\n").split("\t")
            qname = line[0]
            flag = int(line[1])
            chr = line[2]
            pos = int(line[3])
            cigar = line[5]
            umi = get_umi(qname)
            is_minus = is_minus_strand(flag)
            if is_minus:
                strand = "-"
            else:
                strand = "+"
            parsed = parse_cigar(cigar)
            #print(len(read_dict), chr)
            if chr != curr_chrom:
                if curr_chrom is not None:
                    # TODO: print all lines in the dictionary
                    for k, v in read_dict.items():
                        print(v, file=of)
                    read_dict.clear() ##reset dictionary for the first read
                curr_chrom = chr
            if umi in umi_list:
                adj_pos = adjust_pos(pos, is_minus, parsed)
                key = (umi, chr, strand, adj_pos)
                if key not in read_dict:
                    read_dict[key]=line
                else:
                    print(line, file=pcr_duplicates)
            else:
                print(line, file=wrong_umis)
    for k, v in read_dict.items():
        print(v, file=of)
    read_dict.clear() ##reset dictionary for the first read

            #print(read_dict)
of.close()
pcr_duplicates.close()
wrong_umis.close()


