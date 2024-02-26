
### velocity field
fvx <- function(t, x, phia, phib, phic, w1=3, w2=30, kel1 = 2, keln, maxPeak){
  t <- ifelse(t<=0, 0, t)
  
  phi <-  phic*(phia*exp(phib*t)+1-phia)
  keln <-  keln*(phia*exp(phib*t)+1-phia)
  b1 <- (kel1-phi) / ((-w1)^2)
  b2 <- (keln-phi) / ((maxPeak-w2)^2)
  
  fv <- case_when(
    x<=w1 ~ b1*(x-w1)^2+phi,
    x>w1 & x<=w2 ~ phi,
    x>w2 & x<=maxPeak ~ b2*(x-w2)^2+phi,
    x > maxPeak ~ b2*(maxPeak-w2)^2+phi
  )
  
  fv
}

### processivity
pvev <-  function(x, t, alpha, lambda, miu, sigma){
  t <- ifelse(t<=0, 0, t)
  
  pvev <- alpha * dexp(x, lambda) + (1-alpha)* dnorm(x, miu, sigma)
  
  pvev
}

cdfvev <- Vectorize(function(x, t, ...){
  mapply(function(x, t, ...) {
    integrate(pvev, 0, x, t, ...)$value
  }, x, t)
})

### pSer2 state per Pol II
upho <- function(x, t, theta, gamma, taoa, taob, u0a, u0b, u0c){
  t <- ifelse(t<=0, 0, t)
  tao <- taoa * exp(taob*t) + 1- taoa

  u0 <- ifelse(t<=u0a, u0b*t+1, u0c*t+u0b*u0a+1-u0c*u0a)
  c <- 1-exp(-theta*gamma)
  upho <- ifelse(x <= maxPeak, 
                 u0 * (tao*(exp(theta*(x-gamma)) + c) + 1-tao), 
                 u0 * (tao*(exp(theta*(maxPeak-gamma)) + c) + 1-tao))

  upho
}

### r_depho
phoTrans <- function(t, rho){
  t <- ifelse(t<=0, 0, t)
  r <- exp(-rho*t)
  
  r
}

### plot functions
plotParaT <- function(timePoint=c(0,20,40,60,90,180), fn, prefix, ...){
  dfMerge <- list()
  for (t in timePoint) {
    df <- data.frame(x = 1:maxl, para = fn(1:maxl, t = t, ...), prefix = paste0(t,"m"))
    dfMerge <- append(dfMerge, list(df))
  }
  
  dfMerge <- rbindlist(dfMerge)
  
  ggplot(dfMerge, aes(x = x*binsize, y = para, color = prefix)) +
    geom_line() +
    scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D")) +
    xlab("position") + ylab(prefix) +
    theme_classic()
}


plotGMergeRe <- function(ReTrack, metaGMerge, timePlot = c(0,20,40,60,90,180)){
  ReTrackT <- ReTrack[ReTrack$t %in% timePlot, ]
  
  p1 <- metaGMerge %>%
    ggplot(., aes(x = x * binsize, y = spline, color = sample)) +
    geom_line() +
    scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D"),
                       breaks = unique(metaGMerge$sample),
                       labels = paste0(timePlot," min")) +
    scale_x_continuous(breaks = c(0,35,maxl), labels = c("TSS","TES","TES+15kb")) +
    xlab(NULL) + ylab("pSer2 ChIP-Seq") + ylim(c(0,32)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    theme(legend.position = "none")
  
  p1.5 <- ggplot(ReTrackT, aes(x = pos, y = signal, color = factor(t))) +
    geom_line() +
    scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D"),
                       breaks = timePlot) +
    scale_x_continuous(breaks = c(0,35,maxl), labels = c("TSS","TES","TES+15kb")) +
    xlab(NULL) + ylab("Overall pSer2 state signals") + ylim(c(0,32)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    theme(legend.position = "none")

  # p2 <- ggplot(ReTrackT, aes(x = pos, y = trackV, color = factor(t))) +
  #   geom_line() +
  #   scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D")) +
  #   scale_x_continuous(breaks = c(0,35,maxl), labels = c("TSS","TES","TES+15kb")) +
  #   xlab(NULL) +
  #   #ylim(c(-0.5,0.5))+
  #   theme_classic() +
  #   theme(axis.text = element_text(color = "black"),
  #         axis.ticks = element_line(color = "black")) +
  #   theme(legend.position = "none")

  p4 <- ggplot(ReTrackT, aes(x = pos, y = trackPho, color = factor(t))) +
    geom_line() +
    scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D"),
                       breaks = timePlot) +
    scale_x_continuous(breaks = c(0,35,maxl), labels = c("TSS","TES","TES+15kb")) +
    xlab(NULL) + ylab("pSer2 state per Pol II") +
    ylim(c(1,2)) +
    theme_classic()+
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    theme(legend.position = "none")
  
  p5 <- ggplot(ReTrackT, aes(x = pos, y = signal / trackPho, color = factor(t))) +
    geom_line() +
    scale_color_manual(values = c("#564592", "#354E87", "#6A9778","#F8DB51", "#db7c26", "#AE433D"),
                       breaks = timePlot) +
    scale_x_continuous(breaks = c(0,35,maxl), labels = c("TSS","TES","TES+15kb")) +
    xlab(NULL) + ylab('Active Pol II signals') +
    #ylim(c(0.5,NA)) +
    theme_classic()+
    theme(axis.text = element_text(color = "black"),
          axis.ticks = element_line(color = "black")) +
    theme(legend.position = "none")
  
  pLegend <- get_legend(
    p1 + guides(color = guide_legend(nrow = 1),
                ) +
      theme(legend.position = "bottom")
  )
  
  pg <- plot_grid(p1, p1.5, p4, p5, ncol = 1, align="hv")
  pgLgd <- plot_grid(pg, pLegend, ncol = 1, rel_heights = c(1, .1))
  
  pgLgd
}
