# Deduper
## Part 1
- Defining the problem:

PCR duplicate is when you have identical molecules produced by PCR. It is necessary to remove them before the downstream analysis to prevent bias. It is easier to remove after the alignment because the .sam file is more manageable in comparison to a FASTQ file. 

- Some considerations: 
    - To be a PCR duplicate, the reads must be at the same alignment position on the genome. It is necessary to compare the chromosome number (RNAME), the position (POS), and the strand (FLAG). 
    - Soft-clipping: it is necessary to parse the CIGAR string and adjust the position before checking the duplicates, otherwise we can fail to identify PCR duplicates.
    - Check for unique molecular index (UMIs): look at the UMIs used in the analysis and discard reads with incorrect UMIs.
    - Adjust position before checking if it is a duplicate. Do not change the read, we need to save it to the file as the original read.


- Develop your algorithm using pseudocode:

    - First sort data using samtools.
    - Adjust position of all the minus strand and the soft-clipped reads.
    - Discard reads with UMIs that are not on the list.
    - Open data structure to store reads. An idea is a tuple with chromosome, start position, strand and UMI. 
    - Iterate over the first read. 
        - Save chromosome, position, strand and UMI in the data structure.
    - Iterate over the second read. 
        - check to see if the data is equal to the data saved in the set. 
            - if is different:
                - Save previous read in the output file and save the new read in the data structure (replacing the old value, so you are not saving all the reads)
            - if it is equal:
                - Discard read and move to the next line.

- Determine high level functions: 

    - Soft-clipping position adjustment:
```
def soft_clipping(pos:str) -> str:
    '''This function parse the CIGAR string and adjusts the position of reads that were soft-clipped.''''
    pass
    return new_pos
```
Test example:

```
Input: 
Reference: CTTTATTA
Read: AGTTATTA
Pos: 113

Ouput:
Reference: TTATTA
Read: TTATTA
Pos: 111
```
