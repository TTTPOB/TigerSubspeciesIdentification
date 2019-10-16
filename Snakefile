import hashlib
from snakemake.utils import validate
# import sys
REF="genome/P_tigris.scaffold.fa"
THR=46
# NAME=sys.argv[1]
configfile: "workfiles/config.json"
ADMRUN=list(range(1,21))
SAMPLEINFO="workfiles/sampleinfo"
SAMPLENAME="workfiles/samplename"
NAME=config["name"]
scaffolds_to_use=config['scaffolds']
scaffold_file_md5=hashlib.md5(scaffolds_to_use.encode('utf8')).hexdigest()
# with open("workfiles/scaffolds_to_use.txt",'r')as f:
#     scaffold_file_md5=hashlib.md5(f.readline().encode('utf8')).hexdigest()
# validate(config, "config.schema.json")
rule all:
    input:
        # scaffoldtouse="workfiles/scaffolds_to_use.txt",
        # fa="genome/subsetted/"+scaffold_file_md5+".fa",
        # subsetvcfs="vcf/subsetted/"+scaffold_file_md5+".vcf.gz",
        # totsv="vcf/subsetted/"+scaffold_file_md5+".vcf.gz",
        # toped="vcf/subsetted/"+scaffold_file_md5+".vcf.gz",
        # ln_s="workfiles/plink/"+scaffold_file_md5+".vcf.gz",
        # plinkped=scaffold_file_md5+".ped",
        # plinkmapfile=scaffold_file_md5+".map",
        # bwa="workfiles/bwa/{name}.bam",
        # markdup="workfiles/markdup/{name}.bam",
        # gvcf="workfiles/haplotypecaller/{name}.vcf.gz",
        # vcf="workfiles/genotypegvcf/{name}.vcf.gz",
        # filteredvcf="workfiles/filtered/{name}.vcf.gz",
        # filteredvcftbi="workfiles/filetered/{name}.vcf.gz.tbi",
        # catbed="workfiles/admixture/{name}.bed",
        # catbim="workfiles/admixture/{name}.bim",
        # catfam="workfiles/admixture/{name}.fam",
        # admqfile=expand("workfiles/admixture/{name}_{admrun}.5.Q", admrun=ADMRUN, name=NAME),
        rplots=expand("workfiles/plot/"+config["name"]+"_{admrun}.pdf", admrun=ADMRUN)
rule subsetgenome:
    input:
        ref=REF,
        # scaffolds=scaffolds_to_use,
    params:
        scaffolds=scaffolds_to_use
    output:
        fa="genome/subsetted/"+scaffold_file_md5+".fa",
    shell:
        r"""
        if [ {params.scaffolds}==all ]; then
            cp {input.ref} {output.fa}
        else
            samtools faidx {input.ref} `sed 's/,/ /g' {params.scaffolds}` -o {output.fa}
        fi
        """
rule bwaindex:
    input:
        ref=rules.subsetgenome.output.fa
    output:
        idx="genome/subsetted"+scaffold_file_md5+".fa.bwt"
    shell:
        r"""
        bwa index {input.ref}
        """
rule bwa:
    input:
        ref=rules.subsetgenome.output.fa,
        fq1="workfiles/fastq/{name}_L001_R1_001.fastq.gz",
        fq2="workfiles/fastq/{name}_L001_R2_001.fastq.gz",
        idx=rules.bwaindex.output.idx
    output:
        bwa="workfiles/bwa/{name}.bam"
    params:
        rg="\"@RG\\tID:{name}\\tSM:{name}\\tPL:illumina\""
    threads: THR
    shell:
        r"""
        bwa mem -t{threads} -M -R {params.rg} {input.ref} {input.fq1} {input.fq2}|samtools view -bh -f3 -F3852 -t{threads} > {output.bwa}
        """ 
rule sortandmarkdup:
    input:
        rules.bwa.output.bwa
    output:
        markdup="workfiles/markdup/{name}.bam"
    threads: THR
    shell:
        r"""
        samtools sort {input} -@{threads} -n -o workfiles/markdup/{wildcards.name}_namesort.bam
        samtools fixmate workfiles/markdup/{wildcards.name}_namesort.bam -@{threads} -m workfiles/markdup/{wildcards.name}_fixmate.bam
        samtools sort workfiles/markdup/{wildcards.name}_fixmate.bam -@{threads} -o workfiles/markdup/{wildcards.name}_positionsort.bam
        samtools markdup workfiles/markdup/{wildcards.name}_positionsort.bam -@{threads} -r -s {output.markdup}
        rm workfiles/markdup/*_*
        samtools index {output.markdup}
        """
rule haplotypecaller:
    input:
        markdup=rules.sortandmarkdup.output,
        ref=rules.subsetgenome.output.fa
    output:
        gvcf="workfiles/haplotypecaller/{name}.vcf.gz",
    threads: THR
    shell:
        r"""
        if [ ! -f 'genome/subsetted/{scaffold_file_md5}.dict' ]; then
            picard    CreateSequenceDictionary.jar \
                    R= {input.ref} \
                    O= genome/subsetted/{scaffold_file_md5}.dict
        fi
        if [ ! -f 'genome/subsetted/{scaffold_file_md5}.fa.fai' ]; then
            samtools    faidx \
                        genome/subsetted/{scaffold_file_md5}.fa
        fi
        gatk3   -Xmx70G \
                -T HaplotypeCaller \
                -R {input.ref} \
                -I {input.markdup} \
                --emitRefConfidence GVCF \
                -o {output.gvcf} \
                -nct {threads} \
                -variant_index_type LINEAR -variant_index_parameter 128000
        """
rule genotypegvcfs:
    input:
        gvcf=rules.haplotypecaller.output.gvcf,
        ref=rules.subsetgenome.output.fa
    output:
        vcf="workfiles/genotypegvcfs/{name}.vcf.gz",
        vcfidx="workfiles/genotypegvcfs/{name}.vcf.gz.tbi"
    threads: THR
    shell:
        r"""
        gatk3   -Xmx70G \
                -T GenotypeGVCFs \
                -R {input.ref} \
                --variant {input.gvcf} \
                -nt {threads} \
                -o {output.vcf} 
        tabix -p vcf {output.vcf} -f
        """

rule filtering:
    input:
        vcf=rules.genotypegvcfs.output.vcf,
        idx=rules.genotypegvcfs.output.vcfidx
    output:
        filteredvcf="workfiles/filteredvcf/{name}.vcf.gz",
        tabixfile="workfiles/filteredvcf/{name}.vcf.gz.tbi",
    threads: THR
    shell:
        r"""
#        bcftools filter --threads {threads} -i {input.vcf} | \
#        bcftools view -m2 -M2 -v snps -Oz -o {output.filteredvcf} --threads {threads} 
        bcftools view {input.vcf} -m2 -M2 -v snps -Oz -o {output.filteredvcf} --threads {threads}
        tabix -p vcf {output.filteredvcf}
        """

rule subsetvcfs:
    input:
        origvcfs="vcf/Tiger_SNP_6th_N32_final.vcf.gz"
    params:
        scaffolds=rules.subsetgenome.params.scaffolds
    output:
        subsettedvcf="vcf/subsetted/"+scaffold_file_md5+".vcf.gz",
        index="vcf/subsetted/"+scaffold_file_md5+".vcf.gz.tbi"
    shell:
        r"""
        if [ {params.scaffolds}==all ]; then
            cp {input.origvcfs} {output.subsettedvcf}
        else
            bcftools view {input.origvcfs} -r `echo {params.scaffolds} |\
            sed 's/ /,/g'` -Oz -o {output.subsettedvcf}
        fi
        tabix -p vcf {output.subsettedvcf}
        """

rule totsv:
    input:
        rules.subsetvcfs.output.subsettedvcf
    output:
        "vcf/subsetted/"+scaffold_file_md5+".tsv"
    shell:
        r"""
        bcftools query {input} -f '%CHROM\t%POS\n' > {output}
        """

rule toped:
    input:
        chrompostsv=rules.totsv.output,
        filteredvcf=rules.filtering.output.filteredvcf,
        filteredvcftbi=rules.filtering.output.tabixfile,
    output:
        "workfiles/{name}.ped"
    script:
        "scripts/20190814_subsetting_gvcf.py"

rule ln_s:
    input:
        "vcf/subsetted/"+scaffold_file_md5+".vcf.gz"
    output:
        "workfiles/plink/"+scaffold_file_md5+".vcf.gz"
    shell:
        r"""
        mkdir -p workfiles/plink
        cd workfiles/plink
        ln -s ../../vcf/subsetted/{scaffold_file_md5}.vcf.gz ./
        """

rule plink:
    input:
        origvcfsubset=rules.ln_s.output
    output:
        ped="workfiles/plink/"+scaffold_file_md5+".ped",
        mapfile="workfiles/plink/"+scaffold_file_md5+".map",
    shell:
        r"""
        plink --vcf {input.origvcfsubset} --recode --aec --out workfiles/plink/{scaffold_file_md5}
        """

rule cat:
    input:
        origped=rules.plink.output.ped,
        newped=rules.toped.output,
        mapfile=rules.plink.output.mapfile,
    output:
        bed="workfiles/admixture/{name}.bed",
        bim="workfiles/admixture/{name}.bim",
        fam="workfiles/admixture/{name}.fam",
    shell:
        r"""
        cat {input.origped} {input.newped} > workfiles/admixture/{wildcards.name}.ped
        cp {input.mapfile} workfiles/admixture/{wildcards.name}.map
        plink --file workfiles/admixture/{wildcards.name} --make-bed --aec --out workfiles/admixture/{wildcards.name}
        cp {output.bim} {output.bim}.bak
        awk '{{$1=0; print $0}}' {output.bim}.bak > {output.bim}
        """

rule admixture:
    input:
        bed=rules.cat.output.bed,
        bim=rules.cat.output.bim,
        fam=rules.cat.output.fam,
    output:
        qfile="workfiles/admixture/{name}_{admrun}.5.Q"
    params:
        seed="{admrun}"
    threads: THR
    shell:
        r"""
        mkdir -p workfiles/admixture/{wildcards.admrun}
        cp {input.bed} {input.bim} {input.fam} workfiles/admixture/{wildcards.admrun}
        cd workfiles/admixture/{wildcards.admrun}
        admixture -j{threads} --cv {wildcards.name}.bed 5 --seed={params.seed}
        mv {wildcards.name}.5.Q ../../../{output.qfile}
        """

rule rplot:
    input:
        qfile=rules.admixture.output.qfile,
        sampleinfo=SAMPLEINFO,
        samplename=SAMPLENAME,
    output:
        plots="workfiles/plot/{name}_{admrun}.pdf"
    script:
        "scripts/plot.R"
