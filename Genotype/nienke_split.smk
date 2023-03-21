import yaml
import os
import json
import sys
import re
import string
import random
from glob import glob
from snakemake.logging import logger
import traceback
import logging


indir = config["i"] # path to input folder
outdir = config["o"] # path to output folder
check_sex = config["check_sex"] # check sex or not
chunks = 100 # how many files dosage should be split into

ciim_tools = "resources/tools"
G1000 = "1000G"
HRC = "HRC"
biotools = "/vol/biotools/bin"
rscript = "/vol/biotools/alma8/R-4.2.0/bin/Rscript"


tools = base + "/tools"
logdir = outdir+"/logs"


ext = ["bed", "bim", "fam"]     # plink
ext2 = ["pgen", "psam", "pvar"] # plink2


onstart:
    print("##### Snakemake starting #####\n")
    print("\t Creating logdir...\n")
    shell("mkdir -p {logdir}")
    #Path(logdir).mkdir(parents=True, exist_ok=True)

    #TODO: check validity api token before start
    #{tools}/imputationbot/bot jobs

### START
rule all:
	input: expand(outdir+"/PostImp/merge_filtered_rsid.{ext}", ext=ext2)

rule split_by_chromosome:
    input:  expand(indir+sample+".{ext}", ext=ext)
    output: outdir+"/QC/"+sample+".harmonized-chr{chr}.vcf"
    log:    logdir+"/split_by_chromosome/chr{chr}.log"
    shell:
        "{ciim_tools}/plink --bfile {indir}/{sample} \
        --chr {wildcards.chr} --recode vcf-iid \
        --out {outdir}/QC/{sample}.harmonized-chr{wildcards.chr}"

rule sort:
    input:  rules.split_by_chromosome.output
    output: outdir+"/QC/"+sample+".harmonized-chr{chr}_sorted.vcf.gz"
    log:    logdir+"/sort/chr{chr}.log"
    shell:
        "{ciim_tools}/bcftools sort {input} | \
        {biotools}/bgzip -c > {output}"

def generate_password(n):
    return ''.join(random.choice(string.ascii_letters + string.digits) for i in range(n))

rule imputation: # add your own api key to  bot for this to work
    input:  expand(rules.sort.output, chr=range(1,23))
    output: #zip = temp(expand(outdir+"/imputation/local/chr_{chr}.zip", chr=range(1,23))),
            vcf = expand(outdir+"/imputation/local/chr{chr}.dose.vcf.gz", chr=range(1,23)),
            info = expand(outdir+"/imputation/local/chr{chr}.info.gz", chr=range(1,23)),
    params: password = generate_password(8),
    log:    logdir+"/imputation.log"
    resources: runtime="48:00:00"  # can take a really long time due to being queued sometimes
    run:
        try:
            # output is created in working dir, therefore 'cd' first
            cmd = "cd %s && %s/imputationbot/imputationbot impute --files %s \
            --refpanel hrc-r1.1 --population eur --name %s --autoDownload \
            --password %s 2>&1 | tee -a %s" % (outdir, tools, input, sample, params.password, log)

            p = subprocess.run(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            out = p.stdout.decode("ascii", errors="ignore")
            print(out)

            search = re.search("job(-\d+){3}", out)
            if search is not None:
                jobid = search.group(0)
            else:
                print("Jobid not found in output. Please rerun imputation job.")
                sys.exit()

            print("Imputation job either completed or failed. Checking job status...")
            # p.returncode is always 0 I think, so checking output instead
            search = re.search("State: (.+)", out)
            if search is not None:
                state = search.group(1)
            else:
                print("Job status not found in output. Please rerun imputation job.")
                sys.exit()
                # could be : Error: 400 - {"success":false,"message":"Only 3 jobs per user can be executed simultaneously."}

            if "Canceled" in state: # checking like this instead of == because it has ansi color formatting
                print("Job was cancelled. Please rerun imputation job.")
                sys.exit()
            elif "Failed" in state:
                print("Job has failed. Please rerun imputation job.")
                sys.exit()
            elif "Success" in state:  # if successful, change folder name and gunzip .info.gz files
                search = re.search("Error:", out)
                if search is not None:
                    cmd = "%s/imputationbot/imputationbot download %s --password %s" % (tools, jobid, params.password)

                    print("Error after job was successful. Attempting to redownload results as failsafe.")
                    print("CMD:", cmd)
                    shell(cmd)
                else:
                    print("Job was successful. Moving files to imputation folder.")
                    shell("mv %s/%s-%s/local/* %s/imputation/local" % (outdir, jobid, sample, outdir))
                    shell("touch {output}")  # to update the modification time
            else:
                cmd = "%s/imputationbot/imputationbot download %s --password %s" % (tools, jobid, params.password)

                print("Unknown job status. Attempting to download results as failsafe.")
                print("CMD:", cmd)
                shell(cmd)
        except Exception as e:
            logging.error(traceback.format_exc())

rule imputation_info:
    input:  outdir+"/imputation/local/chr{chr}.info.gz"
    output: outdir+"/imputation/local/chr{chr}.info"
    log:    logdir+"/imputation_info_chr{chr}.log"
    shell:  "gunzip {input} -fk | tee -a {log}"

rule pfile:
    input:  outdir+"/imputation/local/chr{chr}.dose.vcf.gz"
    output: expand(outdir+"/PostImp/chr{{chr}}.{ext}", ext=ext2)
    log:    logdir+"/pfile/{chr}.log"
    resources: mem_mb=20000
    shell:
        "plink2 --vcf {input} dosage=HDS --make-pgen \
        --out {outdir}/PostImp/chr{wildcards.chr}" #--double-id

rule filter:
    input:  vcf = rules.pfile.output,
            info = rules.imputation_info.output #outdir+"/imputation/local/chr{chr}.info"
    output: expand(outdir+"/PostImp/filtered_chr{{chr}}.{ext}", ext=ext2)
    log:    logdir+"/filter/chr{chr}.log"
    resources: mem_mb=10000
    shell:
        "plink2 --pfile {outdir}/PostImp/chr{wildcards.chr} \
         --qual-scores {input.info} 7 1 1 --qual-threshold 0.5 --make-pgen \
         --out {outdir}/PostImp/filtered_chr{wildcards.chr}"

rule merge:
    input:  expand(rules.filter.output, chr=range(1,23))
    output: expand(outdir+"/PostImp/merge.{ext}", ext=ext2)
    log:    logdir+"/merge.log"
    shell:
        """
        ls {outdir}/PostImp/filtered_chr*.pgen | sed 's/.pgen//g' > \
        {outdir}/PostImp/merge.list

        plink2 --pmerge-list {outdir}/PostImp/merge.list \
        --out {outdir}/PostImp/merge
        """

rule filter_merged:
    input:  rules.merge.output
    output: expand(outdir+"/PostImp/merge_filtered.{ext}", ext=ext2)
    params: maf = config['filter_maf'],
            hwe = config['filter_hwe']
    log:    logdir+"/filter_merged.log"
    shell:
        "plink2 --pfile {outdir}/PostImp/merge --maf {params.maf} \
        --hwe {params.hwe} --make-pgen --out {outdir}/PostImp/merge_filtered"

rule rsid:
    input:  rules.filter_merged.output
    output: expand(outdir+"/PostImp/merge_filtered_rsid.{ext}", ext=ext2)
    log:    logdir+"/rsid.log"
    shell:
        """
        zcat {HRC}/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz | \
        awk \'{{print $1\":\"$2\":\"$4\":\"$5\"\t\"$3}}\' > \
        {outdir}/PostImp/rsID_chrpos.txt

        awk \'{{print $1\" \"$1\"%\"$2}}\' {outdir}/PostImp/rsID_chrpos.txt > \
        {outdir}/PostImp/rsID_chrpos_pasted.txt

        plink2 --pfile {outdir}/PostImp/merge_filtered \
        --update-name {outdir}/PostImp/rsID_chrpos_pasted.txt --make-pgen \
        --out {outdir}/PostImp/merge_filtered_rsid
	"""
