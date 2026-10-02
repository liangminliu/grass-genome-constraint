import sys
import os.path

def rates_add_coords(inrates,outbed):
    end_pos_bed = sum(1 for line in open(inrates))
    end_pos = end_pos_bed - 1
    chrom = os.path.basename(inrates)[:-10]
    startcount = 0
    endcount = 1
    with open(outbed,'a') as out:
        with open(inrates,'r') as gerp:
            for line in gerp:
                line = line.strip()
                line = line.split("\t")
                ExpS = line[0]
                RejS = line[1]
                print(*[chrom,startcount,endcount,ExpS,RejS],sep="\t",file=out)
                startcount += 1
                endcount += 1
rates_add_coords(snakemake.input[0], snakemake.output[0])



