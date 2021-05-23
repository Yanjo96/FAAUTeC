library("stringr")
library("reshape2")
library("tidyverse")

library("devtools")
install_github("vqv/ggbiplot")
require(ggbiplot)

########################################################################################################
# Remove special characters (* and s) from the data frame
########################################################################################################
removeSpecials <- function(x){
  as.numeric(str_replace_all(x, pattern = "[*s]", replacement = ""))
}

########################################################################################################
# Read in the CSV tables created by FAAUTeC
########################################################################################################
readRun <- function(path){
  run <- read.csv(path, header=T)
  genes <- run$gene
  #genes[43] <- "c125spp"
  rownames(run) = genes
  run = run[,-1]
  run = lapply(run, removeSpecials)
  run <- data.frame(run)
  rownames(run) = genes
  return(run)
}

########################################################################################################
# Remove additional information from the gene names such like "_ge and sorts the data frame
########################################################################################################
kuemmerDich <- function(df, marker = "_test"){
  rownames(df) <- str_replace_all(rownames(df),marker,"")
  rownames(df) <- str_replace_all(rownames(df),"_auto","")
  df$genes <- rownames(df)
  df <- subset(df[order(df$genes),], select = -genes)

  return(df)
}

########################################################################################################
# Sort the dataframes for hypotheses
########################################################################################################
hypoN <- function(theData, nRuns, nAU, hypo){
  HypoN <- list()
  for (k in 0:(nAU-1)){
    for (l in 1:nRuns){
      HypoN[[l+k*5]] <- data.frame(theData[[l]][,hypo+k])
    }
  }
  HypoN <- do.call(cbind, HypoN)
  rownames(HypoN) <- rownames(theData[[1]])
  return( as.data.frame(t(HypoN)))

}

########################################################################################################
# Plottet any number of genes in a direct AU Test comparison plot
########################################################################################################
plotter <- function(genes, dataIQTree, dataRAxML, nHypo, path, nPlots, nCols, hypoName, w, h, showTitle = FALSE, showLegend = FALSE, textsize=12){
  nRuns <- length(dataIQTree)

  if(all(genes == "all")) genes <- rownames(dataIQTree[[1]])
  if(nPlots == "all")  nPlots <- nrow(dataIQTree[[1]])

  IQTree_Hypos <- list()
  RAxML_Hypos <- list()

  ## create dataframes with Hypotheses
  for (i in 1:nHypo){
    IQTree_Hypos[[i]] <- hypoN(dataIQTree, nRuns, 3, (i-1)*3+1)
    RAxML_Hypos[[i]] <- hypoN(dataRAxML, nRuns, 3, (i-1)*3+1)
  }

  Hypos <- list()
  for(j in 1:length(genes)){
    test <- list()
    for (i in 1:nHypo){
      test[[i]] <- data.frame(IQTree_Hypos[[i]][genes[j]],
                              RAxML_Hypos[[i]][genes[j]]
                   )
    }

    test <- do.call(cbind, test)

    colnames(test) <- paste(genes[j], paste0(c("IQTree ","RAxML "), rep(hypoName, each=2)))

    test$AU_Test <- c(rep("CONS", nRuns), rep("IQ1", nRuns), rep("IQ2", nRuns))

    for (i in 1:(nHypo*2)){
      if(sum(test[which(test$AU_Test == "CONS"),i]) == 0){
        test[which(test$AU_Test == "CONS"),i] <- NaN
      }
    }

    test <- melt(test, variable.name = "Hypo")

    Hypos[[j]] <- test
  }

  if(nCols >= 2){
    h1 <- paste(rep(genes[1:nCols], each=2*nHypo), rep(c("IQTree","RAxML"), each=1), rep(hypoName, each=2))
    h2 <- paste(rep(genes[(nCols+1):length(genes)], each=2*nHypo), rep(c("IQTree","RAxML"), each=1), rep(hypoName, each=2))

    lvl <- c()
    for(i in seq(1,length(h1), 2)){lvl <- c(lvl, h1[i:(i+1)], h2[i:(i+1)])}


    Hypos <- do.call(rbind,Hypos)
    Hypos$Hypo <- factor(Hypos$Hypo, levels = lvl)

    t <- ggplot(Hypos) + aes(x=AU_Test,y=value, color=AU_Test) + geom_boxplot(width = 0.6, show.legend = showLegend) + facet_wrap(~Hypo, dir = "h", ncol=nCols*2) + geom_hline(yintercept=0.05, color="darkred") +
      labs(x = "", y = "p-value") + theme_bw() + theme(text = element_text(size = textsize), plot.margin = unit(c(0.05, 0.05, -0.5, 0), "cm")) + scale_color_manual(values=colortheme)

    if(!showLegend) t <- t + theme(legend.position = "none")

    if(showTitle) t <- t + labs(title = paste(genes, collapse=" "))

    ggsave(filename = paste0(path,paste(genes, collapse="_"),".pdf"), plot = t, width=w, height=h, units="cm")
  }else{
    for(i in 1:nPlots){
      first <- (i-1)*floor(length(genes)/nPlots)+1
      second <- floor(length(genes)/nPlots * i)

      #w = 5
      #h = nHypo*length(genes[first:second])+1

      plotHypos <- do.call(rbind, Hypos[first:second])

      t <- ggplot(plotHypos) + aes(x=AU_Test,y=value, color=AU_Test) + geom_boxplot(width = 0.6, show.legend = showLegend) + facet_wrap(~Hypo, dir = "h", ncol=2) + geom_hline(yintercept=0.05, color="darkred") +
                                labs(x = "", y = "p-value") + theme_bw() + theme(text = element_text(size = textsize)) + scale_color_manual(values=colortheme)

      if (!showLegend) t <- t + theme(legend.position = "none")

      if(showTitle) t <- t + labs(title = paste(genes[first:second], collapse=" "))


      ggsave(paste0(path,paste(genes[first:second], collapse="_"),".pdf"),t, width=w, height=h, units="cm")
    }
  }
}

########################################################################################################
# Plotted for all genes each AU Test outcome separated for RAxMl and IQTree
########################################################################################################
derAllesInEinemPlotter <- function(dataIQTree, dataRAxML, nHypo, path, h, w, textsize = 12){
  nRuns <- length(dataIQTree)
  rnames <- rownames(dataRAxML[[1]])
  rnames[which(rnames=="concatenated_125spp")] <- "c125spp"
  for (j in seq(1,3*nHypo,3)){
    CONS_RAxML <- list()
    IQT1_RAxML <- list()
    IQT2_RAxML <- list()

    CONS_IQTree <- list()
    IQT1_IQTree <- list()
    IQT2_IQTree <- list()

    for (i in 1:nRuns){
      CONS_RAxML[[i]] <- dataRAxML[[i]][j]
      IQT1_RAxML[[i]] <- dataRAxML[[i]][j+1]
      IQT2_RAxML[[i]] <- dataRAxML[[i]][j+2]

      CONS_IQTree[[i]] <- dataIQTree[[i]][j]
      IQT1_IQTree[[i]] <- dataIQTree[[i]][j+1]
      IQT2_IQTree[[i]] <- dataIQTree[[i]][j+2]
    }

    CONS_RAxML <- do.call(cbind, CONS_RAxML)
    IQT1_RAxML <- do.call(cbind, IQT1_RAxML)
    IQT2_RAxML <- do.call(cbind, IQT2_RAxML)

    CONS_IQTree <- do.call(cbind, CONS_IQTree)
    IQT1_IQTree <- do.call(cbind, IQT1_IQTree)
    IQT2_IQTree <- do.call(cbind, IQT2_IQTree)

    rownames(CONS_RAxML) = rnames
    CONS_RAxML <- as.data.frame(t(CONS_RAxML))

    for (k in 1:length(CONS_RAxML)){
      if(sum(CONS_RAxML[,k]) == 0){
        CONS_RAxML[,k] <- NaN
      }
    }
    CONS_RAxML$au <- rep("CONS", nRuns)
    CONS_RAxML <- reshape2::melt(CONS_RAxML)

    rownames(CONS_IQTree) = rnames
    CONS_IQTree <- as.data.frame(t(CONS_IQTree))

    for (k in 1:length(CONS_IQTree)){
      if(sum(CONS_IQTree[,k]) == 0){
        CONS_IQTree[,k] <- NaN
      }
    }
    CONS_IQTree$au <- rep("CONS", nRuns)
    CONS_IQTree <- melt(CONS_IQTree)

    rownames(IQT1_RAxML) = rnames
    IQT1_RAxML <- as.data.frame(t(IQT1_RAxML))
    IQT1_RAxML$au <- rep("IQT1", nRuns)
    IQT1_RAxML <- melt(IQT1_RAxML)

    rownames(IQT1_IQTree) = rnames
    IQT1_IQTree <- as.data.frame(t(IQT1_IQTree))
    IQT1_IQTree$au <- rep("IQT1", nRuns)
    IQT1_IQTree <- melt(IQT1_IQTree)

    rownames(IQT2_RAxML) = rnames
    IQT2_RAxML <- as.data.frame(t(IQT2_RAxML))
    IQT2_RAxML$au <- rep("IQT2", nRuns)
    IQT2_RAxML <- melt(IQT2_RAxML)

    rownames(IQT2_IQTree) = rnames
    IQT2_IQTree <- as.data.frame(t(IQT2_IQTree))
    IQT2_IQTree$au <- rep("IQT2", nRuns)
    IQT2_IQTree <- melt(IQT2_IQTree)

    DF_IQTree <- rbind(CONS_IQTree,IQT1_IQTree,IQT2_IQTree)
    DF_IQTree$mlcalc <- "IQTree"

    DF_RAxML <- rbind(CONS_RAxML,IQT1_RAxML,IQT2_RAxML)
    DF_RAxML$mlcalc <- "RAxML"

    DF <- rbind(DF_IQTree, DF_RAxML)
    t <- ggplot(DF) + aes(x=variable,y=value, color=au) + geom_boxplot() + theme_bw() +
                    theme(legend.position = "top", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),
                          legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.05, 0.5, 0), "cm"), axis.title.y = element_blank(), panel.spacing.x = unit(0.2, "lines")) +
                    geom_hline(yintercept=0.05, color="darkred") +
                    labs(x = "", y = "p-value", fill="AU Test:") + coord_flip() +scale_color_manual(values=colortheme)  + facet_wrap(~mlcalc, dir = "h", ncol=2)


    ggsave(filename = paste0(path,"Hypo_",ceiling(j/3),".pdf"), plot = t, width=w, height=h, units="cm")

  }
}

########################################################################################################
# returns the standard deviation, for more than 1 run it returns the standard deviation of the mean
########################################################################################################
getSD <- function(data, gene, dings, nRuns){
  means <- c()
  for (i in dings){
    values <- c()
    for (j in 1:nRuns){
      if (data[[j]][gene, i] > 0){
        values <- c(values, data[[j]][gene, i])
      }
    }
    if(!is.null(values))
      means <- c(means, mean(values))
  }
  return(sd(means))
}

########################################################################################################
# plotted the standard deviation
########################################################################################################
sdPlotter <- function(dataRAxML, dataIQTree, hypo, nRuns, path, dataset){
  r <- list()
  data <- list(dataRAxML, dataIQTree)
  mlcalc <- c("RAxML","IQTree")
  nHypo = length(hypo)
  for (i in seq(1,length(dataRAxML[[1]][,1])*(2*nHypo),(2*nHypo))){
    j <- ceiling(i/(2*nHypo))
    genename <- rownames(dataRAxML[[1]][j,])
    for (k in 1:2){
      for (l in 1:nHypo) {
        r[[i+((l-1)+((k-1)*nHypo))]] <- data.frame(value = getSD(data = data[[k]], gene = j, dings = seq((l*3)-2, l*3), nRuns),
                                                   gene = genename,
                                                   mlcalc = mlcalc[k],
                                                   hypo = hypo[l])
      }
    }
  }
  r <- do.call(rbind,r)
  t <- ggplot(r) + aes(x=hypo,y=value, color=hypo) + geom_boxplot() +
    facet_wrap(~mlcalc, dir = "h", ncol=2) + theme_bw() + scale_color_manual(values=colortheme) +
    labs(x = "", y = "sd of the p-values", fill="AU Test:") +
    theme(legend.position = "none", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.05, 0.5, 0), "cm"), panel.spacing.x = unit(0.2, "lines"))

  if(nRuns == 1){
    t <- t + ggtitle(paste("p-value standard deviation of", dataset, "over", nRuns,"run"))
  }else{
    t <- t + ggtitle(paste("p-value standard deviation of", dataset, "over", nRuns,"runs"))
  }

  ggsave(filename = paste0(path,"variance_",dataset,"_",nRuns,"Runs.pdf"), plot = t, width=15, height=8, units="cm")

  for (i in mlcalc){
    for (j in hypo){
      x <- r[which(r$mlcalc == i & r$hypo == j),]
      print(paste(i,j))
      print(paste("Median:", median(x$value), "Mean:", mean(x$value)))
      print(paste(unlist(subset(x[order(x$value, decreasing = TRUE)[1:10],], select=gene))))
      print(paste(unlist(subset(x[order(x$value, decreasing = TRUE)[1:10],], select=value))))
      print("")
    }
  }

}

########################################################################################################
# returns the means for each AU Tester
########################################################################################################
meaner <- function(data, gene, values, nRuns){
  t <- list()
  for (i in 1:nRuns){
    t[[i]] <- data[[i]][gene,values]
  }
  t <- do.call(rbind,t)
  return(unlist(lapply(t, mean)))
}

########################################################################################################
# returns the means for each AU Tester
########################################################################################################
medianer <- function(data, gene, values, nRuns){
  t <- list()
  for (i in 1:nRuns){
    t[[i]] <- data[[i]][gene,values]
  }
  t <- do.call(rbind,t)
  return(unlist(lapply(t, median)))
}

########################################################################################################
# Plots an heatmap for all genes with most likely hypothesis
########################################################################################################
whichHypo <- function(dataRAxML, dataIQTree, nRuns, hypos, path, w, h, color, threshold){
  nHypo <- length(hypos)
  genes <- rownames(dataRAxML[[1]])
  
  hypos <- c(hypos, "reject")
  color <- c(color, "black")
  
  whichHypo <- list()
  for (i in seq(1,length(genes)*6,6)){
    gene <- genes[ceiling(i/6)]
    IQTree_CONS <- meaner(dataIQTree, gene, seq(1, nHypo*3, 3), nRuns)
    IQTree_IQT <- meaner(dataIQTree, gene, seq(2, nHypo*3, 3), nRuns)
    IQTree_IQT2 <- meaner(dataIQTree, gene, seq(3, nHypo*3, 3), nRuns)

    RAxML_CONS <- meaner(dataRAxML, gene, seq(1, nHypo*3, 3), nRuns)
    RAxML_IQT <- meaner(dataRAxML, gene, seq(2, nHypo*3, 3), nRuns)
    RAxML_IQT2 <- meaner(dataRAxML, gene, seq(3, nHypo*3, 3), nRuns)

    if(any(c(IQTree_CONS, IQTree_IQT, IQTree_IQT2, RAxML_CONS, RAxML_IQT, RAxML_IQT2) >= threshold)) {
    
      whichHypo[[i+0]] <- data.frame(gene = gene, au = "IQTree CONS", hypo = hypos[which.max(c(IQTree_CONS, threshold-0.001))])
      whichHypo[[i+1]] <- data.frame(gene = gene, au = "IQTree IQT1", hypo = hypos[which.max(c(IQTree_IQT, threshold-0.001))])
      whichHypo[[i+2]] <- data.frame(gene = gene, au = "IQTree IQT2", hypo = hypos[which.max(c(IQTree_IQT2, threshold-0.001))])

      whichHypo[[i+3]] <- data.frame(gene = gene, au = "RAxML CONS", hypo = hypos[which.max(c(RAxML_CONS, threshold-0.001))])
      whichHypo[[i+4]] <- data.frame(gene = gene, au = "RAxML IQT1", hypo = hypos[which.max(c(RAxML_IQT, threshold-0.001))])
      whichHypo[[i+5]] <- data.frame(gene = gene, au = "RAxML IQT2", hypo = hypos[which.max(c(RAxML_IQT2, threshold-0.001))])
    }
  }

  
  whichHypo <- do.call(rbind, whichHypo)
  whichHypo[which(whichHypo$gene=="concatenated_125spp"),1] <-"c125spp"
  #print(whichHypo)

  t<- ggplot(whichHypo) + aes(au, gene, fill=hypo) + geom_tile() + theme_bw() + scale_fill_manual(values=color) +
    labs(x = "", y = "", fill="Hypo:") +
    theme(legend.position = "top", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.2, 0.5, 0), "cm"), panel.spacing.x = unit(0.2, "lines"))

  ggsave(filename = path, plot = t, width=w, height=h, units="cm")

}


########################################################################################################
# Plots the biplots
########################################################################################################
biplotter <- function(dataRAxML, dataIQTree, nRuns, dataset, hypos, path, w = 17, h = 10, threshold = 0.05){
  nHypo = length(hypos)
  genes <- c()
  rax <- list()
  iqt <- list()
  for(gene in 1:nrow(dataRAxML[[1]])){
    rax_gene <- meaner(data = dataRAxML, gene = gene, values = 1:(nHypo*3), nRuns = nRuns)
    iqt_gene <- meaner(data = dataIQTree, gene = gene, values = 1:(nHypo*3), nRuns = nRuns)
    
    if (any(c(rax_gene,iqt_gene) > threshold)){
      rax[[gene]] <- rax_gene
      iqt[[gene]] <- iqt_gene

      genes <- c(genes, rownames(dataRAxML[[1]])[gene])
    }
  }
  rax <- do.call(rbind, rax)
  iqt <- do.call(rbind, iqt)
  
  for(i in seq(1,nHypo*3,3)){
    hypo <- do.call(cbind, list(iqt[,i:(i+2)], rax[,i:(i+2)]))
    colnames(hypo) <- c("IQT_CONS", "IQT_IQT1", "IQT_IQT2", "RAx_CONS", "RAx_IQT1", "RAx_IQT2")
    
    hypo_iqt <- iqt[,i:(i+2)]
    hypo_rax <- rax[,i:(i+2)]
    
    colnames(hypo_iqt) <- c("CONSEL", "IQTree1", "IQTree2")
    colnames(hypo_rax) <- c("CONSEL", "IQTree1", "IQTree2")
    
    rownames(hypo) <- genes
    rownames(hypo)[which(rownames(hypo) == "concatenated_125spp")] <- "c125spp"
    
    rownames(hypo_iqt) <- genes
    rownames(hypo_iqt)[which(rownames(hypo_iqt) == "concatenated_125spp")] <- "c125spp"
    
    rownames(hypo_rax) <- genes
    rownames(hypo_rax)[which(rownames(hypo_rax) == "concatenated_125spp")] <- "c125spp"
    
    
    pca <- prcomp(hypo)
    pca_iqt <- prcomp(hypo_iqt)
    pca_rax <- prcomp(hypo_rax)
    
    bi <-  ggbiplot(pcobj = pca,
                    choices = c(1,2),
                    obs.scale = 1, var.scale = 1,
                    labels = rownames(hypo)
                    )
    
    bi <- bi + labs(title = paste(dataset, "Hypothesis: ", hypos[ceiling(i/3)])) +
      theme(legend.position = "top", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.2, 0.5, 0), "cm"), panel.spacing.x = unit(0.2, "lines")) +
      theme_bw()
    
    bi_iqt <-  ggbiplot(pcobj = pca_iqt,
                    choices = c(1,2),
                    obs.scale = 1, var.scale = 1,
                    labels = rownames(hypo_iqt))
    
    bi_iqt <- bi_iqt + labs(title = paste(dataset, "Hypothesis: IQTree, ", hypos[ceiling(i/3)])) +
      theme(legend.position = "top", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.2, 0.5, 0), "cm"), panel.spacing.x = unit(0.2, "lines")) +
      theme_bw()
    
    bi_rax <-  ggbiplot(pcobj = pca_rax,
                    choices = c(1,2),
                    obs.scale = 1, var.scale = 1,
                    labels = rownames(hypo_rax)
                    )
    
    bi_rax <- bi_rax + labs(title = paste(dataset, "Hypothesis: RAxML, ", hypos[ceiling(i/3)])) +
      theme(legend.position = "top", text = element_text(size = textsize), legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,-10,-10,-10), plot.margin = unit(c(0.5, 0.2, 0.5, 0), "cm"), panel.spacing.x = unit(0.2, "lines")) +
      theme_bw()
    
    
    png(file = paste0(path, "/hypo_", hypos[ceiling(i/3)], ".png"), width=w, height=h, unit="cm", res=200)
    print(bi)
    dev.off()
    
    png(file = paste0(path, "/hypo_", hypos[ceiling(i/3)], "_ca", ".png"), width=w, height=h, unit="cm", res=200)
    plot(ca(hypo), main = paste(dataset, hypos[ceiling(i/3)]))
    dev.off()
    
    png(file = paste0(path, "/IQT_hypo_", hypos[ceiling(i/3)], ".png"), width=w, height=h, unit="cm", res=200)
    print(bi_iqt)
    dev.off()
    
    png(file = paste0(path, "/IQT_hypo_", hypos[ceiling(i/3)], "_ca", ".png"), width=w, height=h, unit="cm", res=200)
    plot(ca(hypo_iqt), main = paste(dataset, hypos[ceiling(i/3)], "IQTree"))
    dev.off()
    
    png(file = paste0(path, "/RAx_hypo_", hypos[ceiling(i/3)], ".png"), width=w, height=h, unit="cm", res=200)
    print(bi_rax)
    dev.off()
    
    png(file = paste0(path, "/RAx_hypo_", hypos[ceiling(i/3)], "_ca", ".png"), width=w, height=h, unit="cm", res=200)
    plot(ca(hypo_rax), main = paste(dataset, hypos[ceiling(i/3)], "RAxML"))
    dev.off()
    
  }
}

colortheme <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#a6cee3","#1f78b4")
colortheme1 <- c("#a6cee3","#1f78b4","#b2df8a","#33a02c")
textsize <- 11
