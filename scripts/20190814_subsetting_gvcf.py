#%%
import vcf
import multiprocessing as mp
import pandas as pd, numpy as np
#%%
# orig = vcf.Reader(filename='./vcfgz/Tiger_SNP_6th_N32_final.vcf.gz')
new = vcf.Reader(filename=snakemake.input[1])
origpostsv=pd.read_table(snakemake.input[0], header=None)
#%%
# def extractsitesold(chrompos):
#     chrom=chrompos[0]
#     pos=chrompos[1]
#     record=next(orig.fetch(chrom=chrom, start=pos-1))
#     genolist=[record.REF, str(record.ALT[0])]
#     print(genolist)
#     # return(genolist)
def extractsitesnew(chrompos):
    chrom=chrompos[0]
    pos=chrompos[1]
    try:
        record=next(new.fetch(chrom=chrom, start=pos-1,end=pos))
        genolist=[record.REF, str(record.ALT[0])]
        # print(genolist)
        genotype=list(record.genotype(snakemake.wildcards["name"])['GT'])
        genotype2=genolist[int(genotype[0])]+"\t"+genolist[int(genotype[2])]
        return(genotype2)
    except StopIteration:
        return("0\t0")
    except ValueError:
        return("0\t0")
    

#%%
def getpos(rownum):
   return((origpostsv.iloc[rownum]))

#%%
### get all chromosome and position from orig data
# origchrompos=list()
# for record in orig:
#     origchrompos.append([record.CHROM, record.POS])
#%%
origchromposit=map(getpos, range(origpostsv.shape[0]))

#%%
origchromposlist=list(origchromposit)

#%%
### from chrompos get genolist
# if __name__=='__main__':
#     with mp.Pool(processes=24) as pool:
pedlist=map(extractsitesnew, origchromposlist)
pedstring="\t".join(pedlist)
#%%
# genolist = list(genomap)
finalstring=snakemake.wildcards["name"]+"\t"+snakemake.wildcards["name"]+"\t0\t0\t0\t0\t"+pedstring
#%%
with open(snakemake.output[0],'w') as f:
    print(finalstring, file = f)

#%%
