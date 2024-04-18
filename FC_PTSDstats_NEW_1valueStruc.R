################################################################################
#################   PTSD - Functional Connectivity Stats   #####################
#############   SS -  started March 2024, last edit:    ###############
################################################################################

# 1. upload packages and functions
rm(list=ls())
library(lme4)
library(ggplot2)
library(grid)
library(lmerTest) # show p values in anova of lmer
library(reshape2)
library(wesanderson) # palette
library(paletteer) # 
library(RColorBrewer)
library(emmeans) # posthoc test avec ajustement des pvalues
library(tidyverse)
library(gapminder)
library(MuMIn)
library(ggpubr)
library(readxl)
library("stringr")
library(rstatix)
#library(gghalves)

source(paste("C:/Users/simula/OneDrive - Aix-Marseille Université/PhD/MyProjects/R", '/summary_data.R', sep = ""), 
       encoding = 'UTF-8', echo=TRUE) # charge cette fonction pour les stats descriptives


############################# 2. read the datasets  #############################    
setwd("//dynaserv/meg/nicolas/PTSD/analysis/FC/FC_tables/New_mean")
path_data <- getwd()
path_RES <- "//dynaserv/meg/nicolas/PTSD/analysis/FC/Rstats/New_mean_by_structure" 

# bands <- c('broad', 'delta', 'theta', 'alpha', 'beta', 'lowgamma')
band <- "broad"
raw <- read_excel(paste(path_data, "FCtable-DMN_CEN-broad-mean-29-Mar-2024.xlsx", sep = "/"), sheet = 1, col_names = TRUE)


for (ii in 1:length(bands)) {
  band <- bands[ii]
  if (band == 'broad') {
    raw <- read_excel(paste(path_data, "FCtable-inv_ni-broad-20-Mar-2024.xlsx", sep = "/"), sheet = 1, col_names = TRUE)
  } else {
    raw <- read_excel(paste(path_data, paste("FCtable-inv_ni", band, "20-Mar-2024.xlsx", sep = "-"), sep = "/"), sheet = 1, col_names = TRUE)
  }
  
  #############    define groups and pre-process data    #########################
  thr_GAD <- 7    # anxiety
  thr_NDDI <- 15  # depression in epilepsy
  thr_PCL <- 31   # PTSD score
  
  raw <- cbind(raw, PTSDgroup=NA)
  raw$PTSDgroup[raw$PCL_5 >= thr_PCL & is.na(raw$PCL_5) == FALSE] <- "PTSD+"
  raw$PTSDgroup[raw$PCL_5 < thr_PCL & is.na(raw$PCL_5) == FALSE] <- "PTSD-"
  
  raw <- cbind(raw, GADgroup=NA)
  raw$GADgroup[raw$GAD_7 >= thr_GAD & is.na(raw$GAD_7) == FALSE] <- "anxiety"
  raw$GADgroup[raw$GAD_7 < thr_GAD & is.na(raw$GAD_7) == FALSE] <- "no_anxiety"
  
  raw <- cbind(raw, DEPRgroup=NA)
  raw$DEPRgroup[raw$NDDI_E >= thr_NDDI & is.na(raw$NDDI_E) == FALSE] <- "depression"
  raw$DEPRgroup[raw$NDDI_E < thr_NDDI & is.na(raw$NDDI_E) == FALSE] <- "no_depression"
  
  ##### melt table
  psychosubj <- "9af5a2bd4d70"
  fc <- raw[raw$subject != psychosubj,]
  fc$involved[str_detect(fc$brain_area, "Anterior-cingulate-cortex")] <- "inv"
  #fc <- fc[is.na(fc$involved) == FALSE,]
  my_colnames <- colnames(fc)    
  my_colnames
  
  fc <- melt(as.data.frame(fc), id=my_colnames[c(1:17,28:31)], measure =my_colnames[seq(19, 27, by = 2)], value.name = "node_strength", variable.name = "chan_selection")
  # also out: fc <- melt(as.data.frame(fc), id=my_colnames[c(1:17,28:31)], measure =my_colnames[seq(19:27)], value.name = "node_strength", variable.name = "chan_selection")
  colnames(fc)  
  fc[, c(9:17,23)] <- lapply(fc[, c(9:17,23)], as.numeric)
  fc[, c(1:8,18:22)] <- lapply(fc[, c(1:8,18:22)], as.factor)
  
  unique(fc$chan_selection)
  fc <- fc[fc$chan_selection != "TOTrest_inv_ni", ]
  fc <- fc[fc$chan_selection != "TOTrest_ni_inv", ]

  # Filter for "inv" channels
  inv_channels <- fc %>% filter(involved == "inv")
  inv_channels <- inv_channels[inv_channels$chan_selection != "TOTrest_ni_ni", ]
  
  # filter for ipsilateral channels
  ipsiL_channels <- fc %>% filter(EpilepsySide == "L" & str_detect(brain_area, "Left"))
  ipsiR_channels <- fc %>% filter(EpilepsySide == "R" & str_detect(brain_area, "Right"))
  bilat_channels <- fc %>% filter(EpilepsySide == "Bilateral")
  ipsi_all <- rbind(ipsiL_channels, ipsiR_channels, bilat_channels)
  
  # Filter for contralateral channels
  contraL_channels <- fc %>% filter(EpilepsySide == "L" & str_detect(brain_area, "Right"))
  contraR_channels <- fc %>% filter(EpilepsySide == "R" & str_detect(brain_area, "Left"))
  contra_all <- rbind(contraL_channels, contraR_channels)
  
  ######### FC normalization(optional)
  ## Calculate the maximum TOTrest_all value for each subject
  # max_totrest_per_subject <- inv_channels %>% group_by(subject) %>%
  #   summarize(max_totrest_all = max(TOTrest_all, na.rm = TRUE))
  # 
  # # Join the max value back to the original dataframe and normalize
  # inv_channels <- inv_channels %>% left_join(max_totrest_per_subject, by = "subject")
  # inv_channels$TOTrest_all_normalized <- inv_channels$TOTrest_all / inv_channels$max_totrest_all
  
  
  #############  correlation with indexes   #############
  ipsi_all <- cbind(ipsi_all, ipsi="ipsi")
  contra_all <- cbind(contra_all, ipsi="contra")
  all <- rbind(ipsi_all, contra_all)
  allregions <- c("Orbito-frontal-cortex", "Insula", "Amygdala", "Hippocampus", "Rhinal-cortex", "Anterior-cingulate-cortex")
  
  for (i in 1:length(allregions)) {
    cur_region <- allregions[i]
    
    rm(sel)
    
      sel <- all[str_detect(all$brain_area, cur_region),]
    
    tmp <- sel[sel$chan_selection == "TOTrest_all",]
    my_colnames <- colnames(sel)
    scores <- melt(as.data.frame(tmp), id=my_colnames[c(1:8,18:23)], measure =my_colnames[c(9:17)], value.name = "score_value", variable.name = "score_name")
    
    ymin <- -0.1
    ymax <- 0.3
    
    titre = paste("Correlation of node strength with different scores - ", cur_region)
    AEcorr <- ggplot(data = scores, aes(x = score_value, y = node_strength)) + 
      xlab("score value") + ylab("Node strength h2") + ggtitle(titre)+
      geom_smooth(method = "lm", se = T, colour="#71AA8E", linewidth = 0.5) +
      facet_wrap(score_name~ipsi, ncol =6) +
      theme_light() +
      theme(plot.title = element_text(size=15, hjust = 0.5),
            axis.title = element_text(size = 12, face = "bold"), 
            axis.text = element_text(colour = "black", size = 12),
            strip.text = element_text(size = 12), strip.background = element_rect(fill = "skyblue2")) + 
      geom_point(size = 0.4) +
      stat_cor(method = "spearman", label.x = 0, label.y = -0.08, p.accuracy = 0.001, r.accuracy = 0.01, size = 4)+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.1), limits = c(ymin, ymax))
    
    print(AEcorr)
    # png(paste(path_RES, "Correlation_with_scores", paste(cur_region, ".png", sep = ""), sep = "/"), width = 30, height = 15, units = "cm",
    #     res = 300)
    # print(AEcorr)
    # dev.off()
    
  }
  
  ####################  All inv pooled together  ########################
  colmag <- 8
  sel_regions <- c("Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")
  DMN <- c("Parahippocampal", "Temporal", "Hippocampus", "Posterior-cingulate-cortex", "Precuneus", "Insula", "prefrontal", "Angular")
  DMN_notPTSD <- c("Parahippocampal", "Temporal", "Posterior-cingulate-cortex", "Precuneus", "prefrontal", "Angular")
  CEN <- c("F2", "SFS")  # "Heschl", "Putamen", "SMA","Precentral",
  neutral <- c("STS", "ITS", "T1", "T2", "T3")
  
  pattern <- paste(neutral, collapse = "|")
  
  
  rm(pwc)
  rm(cur_data)
  rm(pattern)
  pattern <- paste(DMN, collapse = "|")
  cur_data <- fc %>%
    filter(str_detect(brain_area, pattern) ) #& chan_selection == "TOTrest_all"
  
  cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  ymin = 0
  
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  
  titre = paste("DMN regions")
  DMN_graph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                        fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  print(DMN_graph)
  
  ############ DMN - pTSD regions
  rm(pwc)
  rm(cur_data)
  rm(pattern)
  pattern <- paste(DMN_notPTSD, collapse = "|")
  cur_data <- fc %>%
    filter(str_detect(brain_area, pattern))
  cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  ymin = 0
  
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  titre = paste("DMN without PTSD regions")
  DMN_notptsd_graph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                   fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  print(DMN_notptsd_graph)
  
  ############   CEN
  rm(pwc)
  rm(cur_data)
  rm(pattern)
  pattern <- paste(CEN, collapse = "|")
  cur_data <- fc %>%
    filter(str_detect(brain_area, pattern)) # & chan_selection == "TOTrest_all"
  cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  ymin = 0
  
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  titre = paste("CEN")
  CEN_graph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                   fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  print(CEN_graph)
  
  ############ neutral regions
  rm(pwc)
  rm(cur_data)
  rm(pattern)
  pattern <- paste(neutral, collapse = "|")
  cur_data <- fc %>%
    filter(str_detect(brain_area, pattern) & chan_selection == "TOTrest_all")
  cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
  ymax = round(max(cur_data$node_strength), 2) +0.03
  ymin = 0
  
  pwc <- cur_data %>% 
    group_by(chan_selection) %>%
    pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
  pwc <- pwc %>% add_xy_position(scales = "free_y")
  
  pwc$y.position <- rep(ymax-0.02, nrow(pwc))
  coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
    group_by(chan_selection) %>%
    cohens_d(node_strength ~ PTSDgroup, paired = F)  
  pwc <- cbind(pwc, coco[,4], coco[,colmag])
  pwc$p <- round(pwc$p, 3)
  pwc$effsize <- round(pwc$effsize, 1)
  
  titre = paste("Neutral regions")
  neutral_graph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                   fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
    #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
    geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
    stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
    xlab("group") + ylab("Node strength") + ggtitle(titre) +
    theme_classic() + 
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    theme(axis.title = element_text(size = 10, face = "bold"), 
          axis.text = element_text(colour = "black", size = 10),
          axis.text.x = element_blank()) +
    theme(legend.text = element_text(size = 10),
          legend.title=element_blank(),
          legend.margin = margin(1,1,1,1, unit = "mm"),
          legend.direction = "horizontal",
          legend.position = "bottom")+
    scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
  # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
  print(neutral_graph)
  

  
  figure <- ggarrange(DMN_graph, DMN_notptsd_graph, CEN_graph, neutral_graph, nrow = 2, ncol=2, labels = c("A", "B", "C", "D"),
                      common.legend = FALSE, legend = "bottom")
  title_fig <- paste("FC in different networks -", band)
  figure <- annotate_figure(figure, top = text_grob(title_fig, size = 12))
  
  print(figure)
  svg(paste(path_RES, band, paste("FC in DMN, CEN, neutral regions -", band, ".svg", sep = ""), sep = "/"), width = 10, height = 10)
  print(figure)
  dev.off()
  print(figure)
  png(paste(path_RES, band, paste("FC in DMN, CEN, neutral regions -", band, ".png", sep = ""), sep = "/"), width = 25, height = 25, units = "cm", res = 700)
  print(figure)
  dev.off()
  
  ########################## IPSILATERAL vs CONTRALATERAL ########################
  colmag <- 8
  sel_regions <- c("Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula")
  DMN <- c("Parahippocampal", "Temporal", "Hippocampus", "Posterior-cingulate-cortex", "Precuneus", "Insula", "prefrontal", "Angular")
  DMN_notPTSD <- c("Parahippocampal", "Temporal", "Posterior-cingulate-cortex", "Precuneus", "prefrontal", "Angular")
  CEN <- c("F2", "SFS")  # "Heschl", "Putamen", "SMA","Precentral",
  neutral <- c("STS", "ITS", "T1", "T2", "T3")
  
  pattern <- paste(neutral, collapse = "|")
  

    rm(pwc)
    rm(cur_data)
    
    ## can choose only between selected inv regions if wanted
    cur_data <- ipsi_all %>%
      filter(str_detect(brain_area, pattern) & chan_selection == "TOTrest_all")
    
    
    # cur_data <- ipsi_all %>%
    #   filter(str_detect(brain_area, pattern) & chan_selection == "TOTrest_all")
    
    #### INVOLVED AND NON-INVOLVED######
    cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
    ymax = round(max(cur_data$node_strength), 2) +0.03
    ymin = 0
    
    pwc <- cur_data %>% 
      group_by(chan_selection) %>%
      pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
    pwc <- pwc %>% add_xy_position(scales = "free_y")
    
    pwc$y.position <- rep(ymax-0.02, nrow(pwc))
    coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
      group_by(chan_selection) %>%
      cohens_d(node_strength ~ PTSDgroup, paired = F)  
    pwc <- cbind(pwc, coco[,4], coco[,colmag])
    pwc$p <- round(pwc$p, 3)
    pwc$effsize <- round(pwc$effsize, 1)
    
    
    titre = paste("Ipsilateral regions")
    invgraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                           fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
      #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
      geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
      stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
      xlab("group") + ylab("Node strength") + ggtitle(titre) +
      theme_classic() + 
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(size = 10, face = "bold"), 
            axis.text = element_text(colour = "black", size = 10),
            axis.text.x = element_blank()) +
      theme(legend.text = element_text(size = 10),
            legend.title=element_blank(),
            legend.margin = margin(1,1,1,1, unit = "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom")+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
    # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
    print(invgraph)
    
    rm(pwc)
    rm(cur_data)
    
    ## can choose only between selected inv regions if wanted
    cur_data <- contra_all %>%
      filter(str_detect(brain_area, pattern) & chan_selection == "TOTrest_all")
    
    # cur_data <- ipsi_all %>%
    #   filter(str_detect(brain_area, pattern) & chan_selection == "TOTrest_all")
    
    #### INVOLVED AND NON-INVOLVED######
    cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
    #ymax = round(max(cur_data$node_strength), 2) +0.03
    ymin = 0
    
    pwc <- cur_data %>% 
      group_by(chan_selection) %>%
      pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
    pwc <- pwc %>% add_xy_position(scales = "free_y")
    
    pwc$y.position <- rep(ymax-0.02, nrow(pwc))
    coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
      group_by(chan_selection) %>%
      cohens_d(node_strength ~ PTSDgroup, paired = F)  
    pwc <- cbind(pwc, coco[,4], coco[,colmag])
    pwc$p <- round(pwc$p, 3)
    pwc$effsize <- round(pwc$effsize, 1)
    
    
    titre = paste("Contralateral regions")
    invgraph_contra <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                          fill = "PTSDgroup", facet.by = c("chan_selection"), outlier.shape = NA) +
      #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
      geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
      stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
      xlab("group") + ylab("Node strength") + ggtitle(titre) +
      theme_classic() + 
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(size = 10, face = "bold"), 
            axis.text = element_text(colour = "black", size = 10),
            axis.text.x = element_blank()) +
      theme(legend.text = element_text(size = 10),
            legend.title=element_blank(),
            legend.margin = margin(1,1,1,1, unit = "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom")+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
    # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
    print(invgraph_contra)
    
    
    
    
    figure <- ggarrange(invgraph, invgraph_contra,  nrow = 1, labels = c("A", "B"),
                        common.legend = TRUE, legend = "bottom")
    title_fig <- paste("FC by PTSD in", pattern, "-", band)
    figure <- annotate_figure(figure, top = text_grob(title_fig, size = 12))
    
    print(figure)
    svg(paste(path_RES, band, paste("neutral regions connectivity ipsicontra -", band, ".svg", sep = ""), sep = "/"), width = 8, height = 5)
    print(figure)
    dev.off()
    print(figure)
    png(paste(path_RES, band, paste("neutral regions connectivity ipsicontra -", band, ".png", sep = ""), sep = "/"), width = 28, height = 20, units = "cm", res = 700)
    print(figure)
    dev.off()
    
    
    ################ EACH REGION #################
  
  colmag <- 8
  allregions <- c("all channels", "Orbito-frontal-cortex", "Amygdala", "Hippocampus", "Anterior-cingulate-cortex", "Rhinal-cortex", "Insula", "Thalamus")
  
  for (i in 1:length(allregions)) {
    cur_region <- allregions[i]
    rm(pwc)
    rm(cur_data)
    
    if (cur_region == "all channels") {
      cur_data <- fc
    } else{
      cur_data <- fc[str_detect(fc$brain_area, cur_region),]
    }
    
    #### INVOLVED AND NON-INVOLVED  ######
    cur_data <- cur_data[is.na(cur_data$node_strength) == FALSE,]
    ymax = round(max(cur_data$node_strength), 2) +0.03
    ymin = 0
    
    pwc <- cur_data %>% 
      group_by(chan_selection) %>%
      pairwise_wilcox_test(node_strength ~ PTSDgroup, p.adjust.method = "fdr", paired = FALSE)
    pwc <- pwc %>% add_xy_position(scales = "free_y")
    
    pwc$y.position <- rep(ymax-0.02, nrow(pwc))
    coco <- cur_data %>%  #|d| < 0.2 "negligible", |d| < 0.5 "small", |d| < 0.8 "medium", otherwise "large"
      group_by(chan_selection) %>%
      cohens_d(node_strength ~ PTSDgroup, paired = F)  
    pwc <- cbind(pwc, coco[,4], coco[,colmag])
    pwc$p <- round(pwc$p, 3)
    pwc$effsize <- round(pwc$effsize, 1)
    
    
    titre = paste(cur_region)
    regiongraph <- ggboxplot(cur_data, x = "PTSDgroup", y = "node_strength",
                           fill = "PTSDgroup", facet.by = "chan_selection", outlier.shape = NA) +
      #geom_line(aes(group = channel), alpha = 0.5, colour = "grey") +
      geom_jitter(position = position_jitter(width = 0.2, height = 0), size = 0.7, alpha = 0.6, colour = 'gray20')+
      stat_pvalue_manual(pwc, hide.ns = FALSE, label =  "{p.adj.signif}: p={p.adj}, d={effsize}", bracket.nudge.y = 0, tip.length = 0.005) +
      xlab("group") + ylab("Node strength") + ggtitle(titre) +
      theme_classic() + 
      theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
      theme(axis.title = element_text(size = 10, face = "bold"), 
            axis.text = element_text(colour = "black", size = 10),
            axis.text.x = element_blank()) +
      theme(legend.text = element_text(size = 10),
            legend.title=element_blank(),
            legend.margin = margin(1,1,1,1, unit = "mm"),
            legend.direction = "horizontal",
            legend.position = "bottom")+
      scale_y_continuous(breaks = seq(ymin, ymax, 0.05), limits = c(ymin, ymax)) 
    # scale_fill_manual(values = c( "#cce3de", "#f6d5f7"))
    print(regiongraph)
    
  }
}
