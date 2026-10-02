import sys
def sitefold2bed(infold,outbed):
    with open(outbed,'a') as out:
        with open(infold,'r') as sites:
            for line in sites:
                if not line.startswith("#"):
                    line = line.strip()
                    line = line.split(" ")
                    chrom = line[0]
                    pos = line[1]
                    fold = int(line[2])
                    # conver 1-indexed to 0-indexed BED position
                    bedstart = int(line[1]) - 1
                    bedstop = int(line[1])
                    if fold == 4:
                        print(*[chrom,bedstart,bedstop],sep="\t",file=out)
sitefold2bed(snakemake.input[0], snakemake.output[0])



