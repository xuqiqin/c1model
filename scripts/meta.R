
gatherMeta <- function(df, info, mode = "mean"){
  outList <- list()
  nskip <- info$upstream[1]/info$`bin size`[1]
  
  for (i in 1:(length(info$sample_boundaries)-1)  ) {
    sample_name <- info$sample_labels[i]
    
    processData <- df[, (7+info$sample_boundaries[i]):(6+info$sample_boundaries[i+1])]
    if (mode == "mean") {
      processDataMedian <- apply(processData, 2, mean)
    }else if(mode == "median"){
      processDataMedian <- apply(processData, 2, median)
    }
    processDataMedian <- processDataMedian[(nskip+1):length(processDataMedian)]
    
    processDataDF <- data.frame(x = 1:length(processDataMedian),
                                gcounts = processDataMedian,
                                sample = sample_name, stringsAsFactors = F)
    rownames(processDataDF) <- NULL
    outList[[sample_name]] <- processDataDF
  }
  
  outList
}

calcPeakEnd <- function(signal, size, genePABin, downstream=10000){
  # calculate peak and gene wave end
  # using groHMM for wave end
  signal <- unlist(signal)
  g <- list()
  g[[1]] <- signal
  
  ## Fit transition and initial probabilities.
  iTrans <- length(signal) - downstream/size
  tProb  <- as.list(data.frame(
    log(c(1-1/iTrans, 1/iTrans) ),
    log(c(0, 1)) ) )  # transition probabilities
  iProb  <- as.double(log(c(1, 0))) # initial probabilities
  
  ## Fit initial distribution paremeters for emission probabilities.
  ePrDist <- c("norm", "norm") 
  parPsi  <- Rnorm(signal[1:iTrans])
  parBas  <- Rnorm(signal[(iTrans+1):length(signal)])
  ePrVars <- data.frame(c(parPsi$mean, sqrt(parPsi$var), -1, -1), 
                        c(parBas$mean, sqrt(parBas$var), -1, -1))
  
  sink("/dev/null")
  BWem <- .Call("RBaumWelchEM", as.integer(2), g, as.integer(1),
                ePrDist, ePrVars, tProb, iProb, 0.01, c(TRUE, TRUE), c(TRUE, TRUE), 
                as.integer(1), TRUE, PACKAGE="groHMM")
  sink()
  
  ansVitervi <- BWem[[3]][[1]]
  DTe <- max(which(ansVitervi==0))
  
  # peak
  closestPeak <- which(signal == max(signal[genePABin:DTe]))
  
  rl <- list(closestPeak = closestPeak, closestEnd = DTe)
  rl
  
}

combineMeta <- function(inputdir, spikeinList, outputdir){
  
  dfC <- fread(inputdir, header = F, stringsAsFactors = F, nrows = 1, sep = "", data.table=F )[1,1]
  info <- fromJSON(gsub("@", "", dfC))
  df <- fread(inputdir, stringsAsFactors = F, sep = "\t", skip = 1, data.table=F)
  dfList <- gatherMeta(df = df, info = info)
  
  size <- info$`bin size`[1]
  genePABin <- info$body[1] / size
  
  if (length(spikeinList) != length(dfList)) {
    stop("Length of the spikein list should be the same as the sample number.")
  }else{
    spikeinList <- spikeinList/spikeinList[1]
  }
  
  metaGRes <- list()
  metaGPara <- list()
  
  for (sample_i in 1:length(dfList)) {
    prefix <- names(dfList)[sample_i]
    
    ### 1. read track
    metaG <- dfList[[sample_i]]
    metaG$gcounts <- metaG$gcounts / spikeinList[sample_i]
    
    ### 2. smooth track: this is optional
    metaGExclude <- metaG[c(1:(info$body[1]/size-5),(info$body[1]/size+1):nrow(metaG)), ]
    expNLS <- smooth.spline(x = metaGExclude$x, y = metaGExclude$gcounts)
    metaG$spline <- predict(expNLS, x=metaG$x)$y
    metaGRes[[sample_i]] <- metaG
    
    ### 3. peak and wave end
    paraOthers <- calcPeakEnd(signal = metaG["spline"], size = size, genePABin = genePABin)
    closestPeak <- paraOthers$closestPeak
    closestEnd <- paraOthers$closestEnd
    metaGPara[[sample_i]] <- data.frame(sample=prefix,
                                        peak=closestPeak,
                                        end=closestEnd,
                                        stringsAsFactors = F)
    
    ### 4. for visualization
    # metaG[,c("x","gcounts","spline")] %>%
    #   reshape2::melt(id.vars = "x", variable.name = "type", value.name = "signal") %>%
    #   ggplot(., aes(x = x, y = signal, color = type)) +
    #   geom_line() +
    #   scale_color_manual(values = c("grey60","blue","red")) +
    #   scale_x_continuous(breaks = c(0, genePABin, closestPeak, closestEnd),
    #                      labels = c("TSS","polyA", "peak", "end")) +
    #   theme_classic()
  }
  
  metaGMerge <- rbindlist(metaGRes) %>% as.data.frame()
  metaGPara <- rbindlist(metaGPara) %>% as.data.frame()

  return(list(metaGMerge, metaGPara, size/1000))
}


