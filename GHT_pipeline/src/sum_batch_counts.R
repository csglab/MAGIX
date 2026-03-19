library(optparse)
set.seed(1)

########################################################   IN and load data ####
option_list = list(
    make_option(c("-a", "--file_with_filenames"), type="character", metavar="character",
                default="./ref/sample_list.txt",
                help=""),

    make_option(c("-b", "--header"), type="character", metavar="character",
                default="test",
                help=""),

    make_option(c("-c", "--outfile"), type="character", metavar="character",
                default="test.txt",
                help="")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser); rm(option_list, opt_parser)

files <- read.csv( file = opt$file_with_filenames, header = FALSE, sep = "\t" )

tmp <- read.csv(file = files[1, "V1"], header = TRUE, sep = "\t" )


for ( i in 2:length(files$V1) ) {
    cat(i); cat("\n")
    this_file <- read.csv(file = files[i, "V1"], header = TRUE, sep = "\t" )
    tmp <- cbind( tmp, this_file )
}

tmp$sum <- rowSums(tmp)

tmp <- data.frame( tmp[, c( "sum" ) ] )
colnames(tmp) <- opt$header

write.csv( x = tmp, file = opt$outfile, quote = FALSE, row.names = FALSE )
