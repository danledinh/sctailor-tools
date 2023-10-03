import pysam
import tqdm
import sys

bamf = sys.argv[1]
outbamf = bamf.replace(".bam",".clean.bam")

if bamf != outbamf:

    with pysam.AlignmentFile(bamf, "rb") as inbam:
        with pysam.AlignmentFile(outbamf, "wb", template=inbam) as outbam:
            for read in tqdm.tqdm(inbam.fetch(until_eof=True)):
                if not (read.is_secondary and read.is_supplementary):
                    tags_d = dict(read.get_tags())
                    if len(tags_d.get("UB")) == 12:
                        outbam.write(read)

    pysam.index(outbamf)
