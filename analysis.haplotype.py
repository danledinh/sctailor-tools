# Haplotyping
import pysam
import pandas as pd
import numpy as np
import glob, os, itertools,statistics
                        
def make_codebank(bam_l):
    """
    For a single gene, for each overlapping read, encode SNV profile.
    bam_l = list, paths of BAM files named using SNV and allele for a single gene
    """
    
    ### get query names from list of BAM files
    qn_l = []
    snv_l = []
    col_l = []
    for i, bamf in enumerate(bam_l):
        snv_id = bamf.split("/")[-1].replace(".bam","")
        col_l.append(snv_id)
        with pysam.AlignmentFile(bamf, "rb") as inbam:
            for read in inbam.fetch(until_eof=True):
                qn_l.append(read.query_name)
                snv_l.append(snv_id)

    ### encode series of SNVs along gene
    qn_df = pd.DataFrame({"qn":qn_l,"snv":snv_id})
    qn_df["value"] = 1
    qn_df = pd.pivot(qn_df, index="qn", columns="snv", values="value").replace(np.nan,0)
    for col in col_l:
        if col not in qn_df.columns:
            qn_df[col] = 0

    snv_colnames = sorted(list(set([x[:-4] for x in qn_df.columns])))
    res_df = pd.DataFrame(index=qn_df.index, columns=snv_colnames)
    for c in snv_colnames:
        vals = [2 if x==1 and y==0 else \
                1 if x==0 and y==1 else \
                0 for x,y in zip(qn_df[f"{c}_ALT"],qn_df[f"{c}_REF"])]
        res_df[c] = vals
    res_df["code"] = res_df.astype(str).apply(''.join, axis=1)
    code_bank = pd.DataFrame(res_df["code"].value_counts())
    code_bank.reset_index(inplace=True)
    code_bank_ll = [[*x] for x in code_bank["index"]]
    code_bank = pd.DataFrame(code_bank_ll,index=code_bank["code"])
    code_bank = code_bank[(code_bank != 0).sum(1) > 1].astype(int)

    return code_bank,res_df


def make_haps(code_bank):
    """
    Parse encoded reads to determine haplotype.
    code_bank = dataframe, generated from `make_codebank` function that contains SNV-encoded reads
    """
    
    ### most abundant phased vars seeds graph
    hap1 = code_bank.iloc[0].tolist()
    hap_pos = list(range(len(hap1)))

    ### identify slots with no var observation to be filled
    z_pos_l = [i for i,x in enumerate(hap1) if x==0]
    nz_pos_l = [i for i,x in enumerate(hap1) if x!=0]

    ### create list of possible known to unknown slot combinations
    poss_combos = list(itertools.product(z_pos_l, nz_pos_l))

    ### loop over unfilled slots
    for z_pos in z_pos_l:

        ### grab all combinations for the unfilled slot
        c_l = [c for c in poss_combos if list(c)[0]==z_pos]

        ### create boolean statement for slicing reads that satify conditions
        bool_s=""
        for i,c in enumerate(c_l):
            if i == 0:
                ### TRUE if contains reads with known slot
                bool_s = code_bank[list(c_l[i])[1]] == hap1[list(c_l[i])[1]]
            else:
                ### update existing for remaining combinations using same logic
                bool_slc = list(code_bank[list(c_l[i])[1]] == hap1[list(c_l[i])[1]])
                bool_s = bool_s|bool_slc

        ### perform slice
        hap1_slc = code_bank[bool_s]       
        hap1_1_sum = sum([y for x,y in zip(hap1_slc[z_pos],hap1_slc.index) if x==1])
        hap1_2_sum = sum([y for x,y in zip(hap1_slc[z_pos],hap1_slc.index) if x==2])

        ### update haplotype based on most abundant allele at previously unknown slot: 2=ALT, 1=REF, 0=UNK
        if hap1_2_sum > hap1_1_sum:
            hap1[z_pos] = 2
        elif hap1_2_sum < hap1_1_sum:
            hap1[z_pos] = 1
        else:
            hap1[z_pos] = 0

    return hap1


def qn2bam(bam2qn_d, outfile):
    """
    Write read to BAM file if found in dictionary of input BAM paths and associated query names.
    bam2qn_d = dict, dictionary of input BAM paths and associated query names
    outfile = str, path to output BAM file
    """
    with pysam.AlignmentFile(outfile, "wb", template=pysam.AlignmentFile(list(bam2qn_d.keys())[0], "rb")) as outbam:
        for bamf in sorted(list(bam2qn_d.keys())): 
            with pysam.AlignmentFile(bamf, "rb") as inbam:
                for read in inbam.fetch(until_eof=True):
                    if read.query_name in bam2qn_d.get(bamf):
                        outbam.write(read)


def make_bam2qn(res_df, hap1, outdir, goi, hapid, bamdir):
    """
    Create intermeidate 
    res_df = dataframe, SNV encodings by query name
    hap1 = list, containing allele codes
    outdir = str, path
    goi = str, gene coordinates
    hapid = str, locus name (usually gene symbol)
    bamdir = str, path of ingested bam files
    """
    hap1_qn = []
    mixed_qn = []
    for i in res_df.drop("code",1).index:
        code = res_df.drop("code",1).loc[i,:].tolist()
        score1 = [q==r for q,r in zip(code,hap1) if q != 0]

        ### majority of observed alleles consistent with computed hap1
        if statistics.mode(score1) == True:
            hap1_qn.append(i)
        else:
            mixed_qn.append(i)

    for suffix,qn_idx in zip([hapid,f"not{hapid}"],[hap1_qn,mixed_qn]):

        ### map reads to bam file (no read dups)
        bam2qn_d = {}
        reads_df = res_df.drop("code",1).loc[qn_idx,:]
        qn_bank = set(reads_df.index.tolist())
        for col in reads_df.columns:
            if len(qn_bank) > 0:
                reads_slc = reads_df.loc[qn_bank,:][col]
                reads_l = reads_slc[reads_slc != 0].index.tolist()
                bam2qn_d[f"{bamdir}/{col}.bam"] = reads_l
                qn_bank = qn_bank-set(reads_l)
        try:
            outfile = f"{outdir}/{goi}_{suffix}.bam"
            qn2bam(bam2qn_d, outfile)
            pysam.sort("-o",outfile+".tmp", outfile)
            os.remove(outfile)
            os.rename(outfile+".tmp", outfile)
            pysam.index(outfile)
        except:
            ## somtimes no reads in key
            continue

    return mixed_qn

def hap_wrap(bam_l,outdir,goi):
    """
    Wrapper to execute haplotype determination for given locus.
    bam_l = list, list of str paths for a given locus
    outdir = str, path
    goi = str, genomic coordinates associated with bam_l
    """
    
    bamdir = bam_l[0].replace(bam_l[0].split("/")[-1],"")
    alt_bam_l = [x for x in bam_l if x.endswith("_ALT.bam")]
    
    ### determine H1
    code_bank,res_df = make_codebank(bam_l)
    hap1 = make_haps(code_bank)
    hap2_qn = make_bam2qn(res_df, hap1, outdir, goi,"H1",bamdir)
    
    ### determine H2
    if len(hap2_qn) > 0:
        ### create new repo of reads that were not mapped to H1
        res_df2 = res_df.loc[hap2_qn,:]
        code_bank = pd.DataFrame(res_df2["code"].value_counts())
        code_bank.reset_index(inplace=True)
        code_bank_ll = [[*x] for x in code_bank["index"]]
        code_bank = pd.DataFrame(code_bank_ll,index=code_bank["code"])
        code_bank = code_bank[(code_bank != 0).sum(1) > 1].astype(int)

        hap2 = make_haps(code_bank)
        hap2_qn = make_bam2qn(res_df2, hap2, outdir, goi,"H2",bamdir)
    else:
        hap2=[]
        
    ### shared SNV identity between H1 and H2 --> convert to 0 (UNK)
    if len(hap2) > 0:
        shared_l = [x==y for x,y in zip(hap1,hap2)]
        if any(shared_l):
            newH1 = [x if y==False else 0 for x,y in zip(hap1,shared_l)]
            newH2 = [x if y==False else 0 for x,y in zip(hap2,shared_l)]

            hap2_qn = make_bam2qn(res_df, hap1, outdir, goi,"H1",bamdir)
            res_df2 = res_df.loc[hap2_qn,:]
            hap2_qn = make_bam2qn(res_df2, hap2, outdir, goi,"H2",bamdir)
            
            
    return (goi,hap1,hap2,[x.split("/")[-1].replace("_ALT.bam","") for x in sorted(alt_bam_l)])


