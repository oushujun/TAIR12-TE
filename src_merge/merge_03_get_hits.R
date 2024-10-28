# This script gets gff file and the genome and return the set of candidate sequences for merging

suppressMessages({
  library(crayon)
  library(pannagram)
  library(optparse)
})

option_list = list(
  make_option("--path.out",      type="character", default="", help="Path to the output folder"),
  make_option("--file.fix",      type="character", default="", help="Path to the sequences file"),
  make_option("--file.fix.seqs", type="character", default="", help="Path to the sequences file")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# print(opt)

# ---- Parsing of parameters ----
opt = parse_args(opt_parser)

if (is.null(opt$path.out)) {
  stop("Error: --path.out is required.")
}
if (is.null(opt$file.fix)) {
  stop("Error: --file.fix is required.")
}

if (is.null(opt$file.fix.seqs)) {
  stop("Error: --file.fix.seqs is required.")
}


path.out = opt$path.out
file.fix = opt$file.fix
file.fix.seqs = opt$file.fix.seqs

# ---- Read the merged sequences

# pokaz("File name:", file.fix)
x = read.table(file.fix, stringsAsFactors = F)
# print(paste0('Number of rows in x', nrow(x)))

x.info = t(sapply(gsub("Chr", "", x$V2), function(s) as.numeric(strsplit(s, '\\|')[[1]][3:5])))
rownames(x.info) = NULL
colnames(x.info) = c('chr', 'beg', 'end')
x = cbind(x, x.info)

x = x[order(-x$end),]
x = x[order(x$beg),]
x = x[order(x$chr),]
rownames(x) = NULL

to.remove <- c()

for (i in 2:nrow(x)) {
  for (j in 1:(i-1)) {
    if ((x$chr[i] == x$chr[j]) &&
        (x$beg[j] <= x$beg[i]) && (x$end[i] <= x$end[j])) {
      to.remove <- cbind(to.remove, c(i, j))
      break
    }
  }
}


if(!is.null(to.remove)){
  to.remove = t(to.remove)
  x = x[-to.remove[,1],]
}

if(nrow(x) == 0){
  pokaz('No sequences as merged')
  return()
}

# Read all of the sequence files and take those sequences, which match in  names

file.seqs <- list.files(path = path.out, pattern = "^merged_seqs", full.names = TRUE)

seqs.all = c()

for(f in file.seqs){
  seqs.all = c(seqs.all,
               readFasta(f))
}

seqs.all = seqs.all[x$V2]

writeFasta(seqs.all, file.fix.seqs)






