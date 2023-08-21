import pysam
import pandas as pd
import tqdm
import sys

bamf = sys.argv[1]

qn_l = []
ql_l = []
ug_l = []
with pysam.AlignmentFile(bamf, "rb") as inbam:
    for read in tqdm.tqdm(inbam.fetch(until_eof=True)):
        tags_d = dict(read.get_tags())
        ug = tags_d.get("UG")
        qname = read.query_name
        qlen = read.query_length

        qn_l.append(qname)
        ql_l.append(qlen)
        ug_l.append(ug)

bam_df = pd.DataFrame({"qn":qn_l,
                       "ql":ql_l,
                       "ug":ug_l
                      })
bam_df = bam_df.sort_values("ql", ascending=False)
bam_df = bam_df.groupby(["ug"]).head(1)
outtsv = bamf.replace(bamf.split("/")[-1],"qname_umitools.txt")
bam_df.loc[:,["qn"]].to_csv(outtsv,
                            sep="\t", header=False, index=False
                           )
pysam.sort("-o", "bamf.bam", "bamf")

