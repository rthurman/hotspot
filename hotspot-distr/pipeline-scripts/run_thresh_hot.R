#! /usr/bin/env Rscript
# -*- mode: R -*-
# Above tells emacs to go into R mode for editing (must be second line).

##
## Find FDR thresholds for hotspots (not peaks) and write out results.
## Run with R --no-save < run_thresh_hot.R
##

thisscr <- 'run_thresh_hot'
cat('\n', thisscr, '\n')

mrgwid <- _MERGE_DIST_
tags <- "_TAGS_"
fdrs <- _FDRS_

rand <- "_RANDIR_"
outdir <- "_OUTDIR_"

rootFindInterval <- c(3,35)
## Turn off warnings for a minute while we check if non-default root-finding interval has been specified.
warnopt <- getOption("warn")
options("warn" = -1)
if(!is.na(as.numeric("_FDR_ROOT_FIND_MIN_")))
  rootFindInterval[1] <- as.numeric("_FDR_ROOT_FIND_MIN_")
if(!is.na(as.numeric("_FDR_ROOT_FIND_MAX_")))
  rootFindInterval[2] <- as.numeric("_FDR_ROOT_FIND_MAX_")
options("warn" = warnopt)

proj <- gsub('.bed.starch$', '', gsub('.bam$', '', basename(tags)))
cat(proj, '\n')
ntags <- read.table(paste(outdir, '/', proj, '-pass1/', proj, '.stdout', sep = ''), as.is = T)[1,2]
ntagsr <- sprintf("%.0f", round(ntags/100000)*100000)
hotof <- paste(outdir, '/', proj, '-both-passes/', proj, '.hotspot.twopass.zscore.wig', sep = '')
hotrf <- Sys.glob(paste(rand, '/', ntagsr, '-ran*both-passes/', ntagsr, '-ran.*hotspot.twopass.zscore.wig', sep = ''))  
if(!file.exists(hotof)){
  warning(hotof, 'does not exist; skipping', immediate. = T)
  next
}
if(!file.exists(hotrf)){
  warning(hotrf, 'does not exist; skipping', immediate. = T)
  next
}
hoto <- read.table(hotof, skip = 1, as.is = T, col.names = c('Chr', 'Start', 'Stop', 'z'))
hotr <- read.table(hotrf, skip = 1, as.is = T, col.names = c('Chr', 'Start', 'Stop', 'z'))

for(fdr in unlist(strsplit(as.character(fdrs), split=" "))){
  if(fdr == "N")
    next
  else if(is.na(as.numeric(fdr))){
    warning('non-numeric value for fdr:', fdr, immediate. = T)
    next
  }
  fdr <- as.numeric(fdr)
  if(fdr == 0)
    thz <- max(hotr$z)
  else{
    res <- uniroot(function(th) sum(hotr$z > th)/sum(hoto$z > th) - fdr, interval = rootFindInterval, tol = .01)
    thz <- res$root
  }
  cat('fdr', fdr, 'z-score threshold =', thz, '-- number of thresholded hotspots =', sum(hoto$z > thz), '\n')
  nm <- paste(proj, '.hotspot.twopass.fdr', fdr, sep = '')
  trnm <- paste(proj, '.hotspot.fdr', fdr, sep ='')
  bed <- paste(outdir, '/', proj, '-both-passes/', nm, '.bed', sep = '')
  wig <- gsub('bed$', 'wig', bed)
  cat('track type=wiggle_0 visibility=full name=', trnm, '\n', sep = '', file = wig)
  mrg <- hoto[hoto$z > thz,]
  write.table(mrg, file = wig, quote = F, sep = '\t', row.names = F, col.names = F, append = T)
  ## For now, write a proper bed file, with score in 5th column.  We need this for the bedmap command to follow.
  write.table(cbind(mrg, mrg$z), file = bed, quote = F, sep = '\t', row.names = F, col.names = F)
  ## Merge nearby thresholded hotspots
  mrgwig <- gsub('bed$', 'merge.wig', bed)
  mrgnm <- paste(nm, '.merge', sep = '')
  cat('track type=wiggle_0 visibility=full name=', trnm, '\n', sep = '', file = mrgwig)
  system(paste('bedops --range ', mrgwid/2, ' -m ', bed, ' | bedops --range -', mrgwid/2, ' -m - | bedmap --delim "\t" --echo --max - ', bed, ' >> ', mrgwig, sep = ''))
  ## Replace previous bed file with 4-column one.
  write.table(mrg, file = bed, quote = F, sep = '\t', row.names = F, col.names = F)
}

