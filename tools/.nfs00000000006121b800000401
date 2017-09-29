
library("protr", verbose = F, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
seq=args[1]

# seq="MILGAVFYIVFIALFFGIAVGIIFAIKSIKLI"

moreau <- extractMoreauBroto(seq)
moran <- extractMoran(seq)
paac <- extractPAAC(seq)
trans <- extractCTDT(seq)

out <- c(as.list(moreau), as.list(moran), as.list(paac), as.list(trans))
write.table(t(as.data.frame(out)), file='', sep="\t", quote=FALSE)