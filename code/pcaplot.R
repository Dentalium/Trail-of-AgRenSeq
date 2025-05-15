library(ggplot2)
library(dplyr)

setwd("proj/junks/50412_agrenseq/result/05_pca_manual/")

value <- read.table("scores.tsv", header = T)
subpop <- read.table("subpop_list.txt", header = T)

colnames(value)[1] <- "accession"
value <- value %>% 
  inner_join(subpop, by = "accession")

ggplot(data = value) +
  geom_point(aes(x=PC1, y=PC2, color=line)) +
  scale_color_manual(name=NULL,
                     values = c("L2"="cyan3", "Inter"="darkgreen", "L1"="magenta2"),
                     labels = c( "Intermediate", "ssp. tauschii", "ssp. strangulata")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.3,0.8)
  )

ggplot(data = value) +
  geom_point(aes(x=PC3, y=PC2, color=line)) +
  scale_color_manual(name=NULL,
                     values = c("L2"="cyan3", "Inter"="darkgreen", "L1"="magenta2"),
                     labels = c( "Intermediate", "ssp. tauschii", "ssp. strangulata"))
