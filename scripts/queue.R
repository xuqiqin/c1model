### given Pol II loads at ti, its position at t
fpos <- function(ti, t, tspan, parms=c()){ # given x(ti)=0, solve x(t)
  txMod <- function (t, x, parms) {
    with(as.list(c(x, parms)), {
      dx <- fvx(t, x)
      return(list(dx))
    })
  }
  
  times <- seq(ti, t, by = tspan)
  sol <- deSolve::ode(y = c(x=0), times = times, func = txMod, parms = parms)
  
  sol[,"x"]
}

fposCheat <- function(ti, t, moveTrack){
  mapply(function(ti, t) {
    moveTrack[as.character(round(ti, 5)), as.character(round(t-ti, 5))]
  }, ti, t)
}

### number of Pol II that drops at [x1, x2]
ReCdfvev <- function(x1, x2, t, ...){
  mapply(function(x1, x2, t, ...) {
    integrate(pvev, x1, x2, t, ...)$value
  }, x1, x2, t)
}


### recycling
ReSingleCell <- function(pfree, thres, moveTrack){
  
  t <- tmin-tspan
  
  # load existing Pol2
  t0 <- ifelse(tmin<=0, tmin - tstable, 0-tstable)
  preLoadTime <- seq(t0+tspan, tmin-tspan, tspan) %>% round(., 10)
  preLoadPol2 <- rep(vini*tspan, length(preLoadTime))
  
  preLoadPos <- fposCheat(preLoadTime, t, moveTrack) #mclapply(pExist, function(x){fpos(x, t)}, mc.cores = 10) %>% unlist() 
  
  pol2t <- data.frame(n = 1:length(preLoadTime), nPol2 = preLoadPol2, pos = preLoadPos, ti = preLoadTime, nDrop = -1)
  nPol2 <- max(pol2t$n)
  
  pol2t <- pol2t[pol2t$pos<maxl, ]  
  rownames(pol2t) <- NULL
  
  combinedf <- data.frame(pol2t, t=t)
  
  # initialize pool, loadn
  pfreeIni <- pfree
  loadn <- NULL
  poolList <- c()
  loadList <- c()
  dropList <- c()
  
  while (t < tmax) {
    cat("\r","t=", t)
    
    t <- round(t + tspan, 10)
    
    # update existed Pol II
    newPos <- fposCheat(pol2t$ti, t, moveTrack) #mclapply(pol2t$ti, function(x){fpos(x, t)}, mc.cores = 10) %>% unlist()
    newDrop <- pol2t$nPol2 * ReCdfvev(pol2t$pos, newPos, t = pol2t$t)
    pol2t$pos <- newPos
    pol2t$nDrop <- newDrop
    
    # generate new Pol II according to pool size and update pool
    loadn <- case_when(
      pfree>pfreeIni ~ vini * tspan * pfree/pfreeIni,
      pfree<=pfreeIni & pfree>=thres ~ vini * tspan,
      pfree<thres ~ vini * tspan * pfree/thres
    )
    newdf <- data.frame(n = nPol2+1,
                        nPol2 = loadn,
                        pos = 0,
                        ti = t,
                        nDrop = 0)
    nPol2 <- nPol2 + 1
    pol2t <- rbind(pol2t, newdf)
    combinedf <- rbind(combinedf, data.frame(pol2t, t = t))
    
    dropn <- sum(pol2t$nDrop)
    pfree <- pfree - round(loadn,5) + phoTrans(t) * round(dropn,5)
    poolList <- c(poolList, pfree)
    loadList <- c(loadList, loadn)
    dropList <- c(dropList, dropn)
    
    # remove dissotiated Pol II
    pol2t <- pol2t[pol2t$pos<maxl, ]
    
  }
  
  list(combinedf, poolList, loadList, dropList)
}

ReCreateTrack <- function(df){
  
  df$trackV <-  fvx(t = df$t, x = df$pos)
  df$trackKeep <-  1- cdfvev(x = df$pos, t = df$t)
  df$trackPho <- upho(x = df$pos, t=df$t)
  df$trackPosDelta <- c(abs(diff(df$pos)),NA)
  df$trackDensity <- df$nPol2 / df$trackPosDelta
  df$trackTheo <- df$trackDensity  * df$trackKeep * df$trackPho
  
  df <- df[!is.na(df$trackTheo), ]
  
  df
}

ReTrackNorm <- function(ReList, metaGMerge, timePlot = c(0,20,40,60,90,180)){
  # create track
  ReTrack <- ReList[ReList$nDrop!= -1, ] %>%
    split(., .$t) %>%
    .[as.character(timePlot)]  %>%
    lapply(., ReCreateTrack) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  
  # norm track
  normBin <- min(ReTrack$pos[ReTrack$t==0])
  ReTrack$signal <- ReTrack$trackTheo/(ReTrack$trackTheo[ReTrack$t==0 & ReTrack$pos == min(ReTrack$pos[ReTrack$t==0])] / metaGMerge$spline[1] )
  
  ReTrack
  
}

