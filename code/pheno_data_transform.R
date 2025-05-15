library(tidyr)
library(dplyr)

setwd("F://projects/34_agrenseq/pheno/")

#s_ttksk <- read.table("score_ttksk.txt", header = T, fill = T)

for (i in c("ttksk","rkqqc","qthjc")) {
  message(i)
  input <- read.table(paste0("score_",i,".txt"), header = T, fill = T)
  
  transformed <- input |> 
    pivot_longer(cols=r1:r14, values_to = "it", names_to = NULL) |> 
    drop_na(it) |> 
    filter(it!="") |> 
    mutate(score=case_match(it,
                            "0;"~2,
                            ";"~1.67,
                            "1-"~1.33,
                            "1"~1,
                            "1+"~0.67,
                            "2-"~0.33,
                            "2"~0,
                            "2+"~-0.33,
                            "3-"~-0.67,
                            "3"~-1,
                            "3+"~-1.33,
                            "4"~-2)) |> 
    group_by(Accession) |> 
    summarise(mean_score=mean(score))
  
  write.table(transformed, paste0("score_",i,".mean.txt"), quote = F, row.names = F, sep = "\t")
}




