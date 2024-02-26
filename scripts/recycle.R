library(ggplot2)
library(data.table)
library(rjson)
library(groHMM)
library(dplyr)
library(cowplot)
library(pbmcapply)
library(default)

rm(list = ls())
source('meta.R')
source('events.R')
source('queue.R')

#########################################################
# inputdir <- '../data/C1D.matrix.gz'
# configdir <- '../data/config.txt'
# spikeinList <- c(19163158, 17625265, 17986917, 18285376, 20451950, 24631397)
# outputdir <- '../results'
#########################################################

config <- read.delim(configdir, header = F, stringsAsFactors = F, comment.char = "#")
for(i in 1:nrow(config)){
  assign(config[i, 1], config[i, 2])
}

inputMeta <- combineMeta(inputdir, spikeinList, outputdir)
metaGMerge <- inputMeta[[1]]
metaGPara <- inputMeta[[2]]
binsize <- inputMeta[[3]]

maxl <- max(metaGMerge$x) * binsize
maxPeak <- metaGPara[1, "peak"] * binsize
maxEnd <- metaGPara[1, "end"] * binsize

# reconstruct functions
default(fvx) <- list(phia = phia, phib = phib, phic = phic, keln = keln, maxPeak = maxPeak)
default(pvev) <- list(alpha = alpha, lambda = lambda, miu = (maxPeak+maxEnd)/2, sigma = (maxEnd-maxPeak)/3)
default(upho) <- list(theta = theta, gamma = gamma, taoa = taoa, taob = taob, u0a = u0a, u0b = u0b, u0c = u0c)
default(phoTrans) <- list(rho = rho)

# generate a cheat list of trajectory for saving time :)
t0 <- ifelse(tmin<=0, tmin - tstable, 0-tstable)
moveSlowest <- tmin-t0
moveTrack <- pbmclapply(seq(t0, tmax, tspan), function(ti){
  fpos(ti, (ti+moveSlowest), tspan)
}, mc.cores = ncores) %>%
  do.call(rbind, .) %>%
  as.data.frame()
rownames(moveTrack) <- as.character(round(seq(t0, tmax, tspan), 10))
colnames(moveTrack) <- as.character(round(seq(0, ceiling(moveSlowest), tspan), 10))


# recycle
ReRes <- ReSingleCell(pfree = pfree, thres = thres, moveTrack = moveTrack)

# normalize and combine track
ReTrack <- ReTrackNorm(ReRes[[1]], metaGMerge)
p1 <- plotGMergeRe(ReTrack, metaGMerge)
p1
ggsave(filename = paste0(outputdir,"/model.tune.theoTrack.pdf"), plot = p1, width = 5, height = 10)

# recycling factor
ReDf <- data.frame(t = seq(tmin, tmax, tspan), pool = ReRes[[2]], nload = ReRes[[3]], ndrop = ReRes[[4]])
ReDf$RF <- ReDf$pool / ReDf$pool[1]
p2 <- ggplot(ReDf, aes(x = t, y = RF))+
  geom_line() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  ylab("Recycling factor") +
  scale_x_continuous(breaks = seq(tmin, 180, by = 30))+
  theme_classic()+
  theme(axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"))
p2
ggsave(filename = paste0(outputdir,"/model.recycle.pdf"), plot = p2, width = 5, height = 3)


