### given Pol II loads at ti, its position at t
fposSto <- function(ti, t, tspan, parms=c()){ # given x(ti)=0, solve x(t)
  txMod <- function (t, x, parms) {
    with(as.list(c(x, parms)), {
      dx <- fvx(t, x)
      return(list(dx))
    })
  }
  
  if ( abs(t-ti) < 1e-5){
    return(0)
  }else{
    times <- seq(ti, t, by = tspan)
    sol <- deSolve::ode(y = c(x=0), times = times, func = txMod, parms = parms)
    
    solx <- sol[,"x"]
    # return(solx)
    return(solx[length(solx)])
  }
}

### stochastic loading of Pol II
randPol2Load <- function(tstart, tend, meanIni, aseed){
  ot <- c()
  t <- tstart
  i <- 0
  while (t < tend) {
    i <- i + 1
    set.seed(aseed + i)
    randt <- rexp(n = 1, rate = meanIni)
    t <- t + randt
    ot <- c(ot, t)
  }
  ot[ot<=tend]
}

### drop distrubution generator
cdfvevSolve <- Vectorize(function(r, t, ...){
  #nleqslv(0, function(x) cdfvev(x, t, ...)-r)$x 
  uniroot(f = function(x){cdfvev(x, t, ...)-r}, interval = c(-1, maxl+1))$root
})

pol2DropPos <- function(loadtime, aseed, ...){
  set.seed(aseed)
  randn <- runif(length(loadtime))
  dropPos <- cdfvevSolve(randn, loadtime, ...)
  
  dropPos
}


### stochastic process for single cell
singleCell <- function(loady = 0, dropy = -1, aseed){
  
  # initialize
  t <- tmin - tspanSto
  t0 <- ifelse(tmin<=0, tmin - tstable, 0-tstable)
  preLoadTime <- randPol2Load(t0, t, viniSto, aseed)

  # load existing Pol2
  preLoadPos <- lapply(preLoadTime, function(x) fposSto(x, t, tspan = (t-x)/20)) %>% unlist()
  preDropPos <- pol2DropPos(preLoadTime, aseed)
  
  pol2t <- data.frame(n = 1:length(preLoadTime), pos = preLoadPos, y = loady, ti = preLoadTime, posd = preDropPos, isNew = F)
  pol2t$upho <- upho(pol2t$pos, t)
  nPol2 <- max(pol2t$n)
  
  pol2t$pos <- ifelse(pol2t$pos>=pol2t$posd, NA, pol2t$pos)
  pol2t <- na.omit(pol2t)
  rownames(pol2t) <- NULL
  
  combinedf <- data.frame()
  
  # initialize pool, loadn
  pfreeIni <- pfreeSto
  loadn <- NULL
  paraList <- data.frame()

  while (t < tmax) {

    t <- t + tspanSto
    bseed <- aseed + t*111
    
    # update existed Pol II
    pol2t <- lapply(1:nrow(pol2t), function(i){
      newPos <- fposSto(pol2t[i, "ti"], t, tspan = (t-pol2t[i, "ti"])/20)
      if (pol2t[i, "y"] == dropy ) { # delete previous dropped pol2
        pol2t[i, "pos"] <- NA
      }else if ( newPos >= pol2t[i, "posd"]) { # change y for dropped pol2
        pol2t[i, "pos"] <- newPos
        pol2t[i, "y"] <- dropy
      }else{
        pol2t[i, "pos"] <- newPos
        pol2t[i, "upho"] <- upho(newPos, t)
      }
      pol2t[i,]
    }) %>% rbindlist() %>% as.data.frame()
    pol2t <- na.omit(pol2t)
    pol2t$isNew <- F
    
    # generate new Pol II according to pool size and update pool
    viniUpdate <- case_when(
      pfreeSto>=thresSto ~ viniSto,
      pfreeSto<thresSto ~ viniSto * pfreeSto/thresSto
    )
    newPol2 <- randPol2Load(t-tspanSto, t, viniUpdate, aseed = bseed)
    
    if (length(newPol2) != 0) {
      newPos <- lapply(newPol2, function(x) fposSto(x, t, tspan = (t-x)/20)) %>% unlist()
      newDrop <- pol2DropPos(newPol2, aseed = bseed)
      loadn <- length(newPol2)
      
      newdf <- data.frame(n = (nPol2+1):(nPol2+loadn), pos = newPos, y = loady, ti = newPol2, posd = newDrop, isNew = T)
      newdf$upho <- upho(newdf$pos, t)
      newdf$y <- ifelse(newdf$pos>=newdf$posd, dropy, newdf$y) # there's also a probability for the newly load pol2 to drop
      
      nPol2 <- nPol2 + loadn
      pol2t <- rbind(pol2t, newdf)
    }
    combinedf <- rbind(combinedf, data.frame(pol2t, t = t))
    
    dropn <- sum(pol2t$y == dropy)
    set.seed(bseed)
    isdrop <- runif(dropn) <= phoTrans(t) # % return to pool
    pfreeSto <- pfreeSto - loadn + sum(isdrop)
    paraList <- rbind(paraList, data.frame(t=t, pool = pfreeSto, loadn=loadn, dropn=dropn, vini=viniUpdate))
    
  }
  
  list(combinedf, paraList)
}

### count by bins
countPol2 <- function(df) {
  df <- df[df$y>=0, ] # remove droppped Pol2
  posHist <- wtd.hist(df$pos, weight = df$upho, breaks = seq(0, maxl, by = binsize), plot = F,right = F)
  o <- posHist$counts
  names(o) <- posHist$breaks[-1]
  
  o/ncell
}

