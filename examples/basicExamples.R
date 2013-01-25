setwd("~/Dropbox/QuasiAsterisk/QuasiViz/")

require("R.utils")
sourceDirectory("src/")

#readFasta("~/Dropbox/Projects/Martin/int/20130117/ALL/08-16621.fasta")
dir <- readDirectory("~/Dropbox/Projects/Martin/int/20130117/ALL/")

plotEntropy(dir,legendPosition="bottom")
plotHapDist(dir)