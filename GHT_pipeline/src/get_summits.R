library(dplyr)
library(optparse)
library(data.table)
set.seed(1)
options(scipen = 999)

option_list = list(
    make_option(c("-i", "--in_bed"), type="character", metavar="character",
                default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/test1/test2/test3.bed",
                help=""),
    make_option(c("-o", "--out_bed"), type="character", metavar="character",
                default="/home/ahcorcha/repos/ahcorcha/Projects/P2_TF_Methyl/bin/codebook_ChIP_seq/data/03_process_GHT_SELEX/test1/test2/test3_out.bed",
                help=""));

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)

in_bed <- data.frame( read.csv(file= opt$in_bed, header = FALSE, sep = "\t") )

in_bed <- in_bed %>% group_by(V4) %>% top_n(1, V5)

in_bed <- setDT(in_bed)[, if(.N >1) .SD[ceiling(.N/2)] else .SD ,V4]

in_bed <- in_bed[, c(2,3,4,5) ]

write.table(x = in_bed, file = opt$out_bed, sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
