#library(ggplot2)
#library(reshape2)
#library( RColorBrewer )
#library( devtools )
#library(circlize)
#library(dplyr)



plotjci_heatmap <- function(a, b){
  
  
  jaccard <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
  }
  
  
  
  ##set working directory
  setwd("~/R/TCR_reparetoire_analysis/stress_TCR/plotfancyvj_sonaf_TRA")
  
  ##import data sets
  #ctrl
  sonaf_TRA_ctrl1 <- read.table("filterd_CTRL1_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_ctrl2 <- read.table("filterd_CTRL2_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_ctrl3 <- read.table("filterd_CTRL3_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_ctrl4 <- read.table("filterd_CTRL4_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_ctrl5 <- read.table("filterd_CTRL5_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_ctrl6 <- read.table("filterd_CTRL6_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  #imm
  sonaf_TRA_imm1 <- read.table("filterd_Imm1_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_imm2 <- read.table("filterd_Imm2_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_imm3 <- read.table("filterd_Imm3_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_imm4 <- read.table("filterd_Imm4_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_imm5 <- read.table("filterd_Imm5_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  sonaf_TRA_imm6 <- read.table("filterd_Imm6_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  
  ##rename_colnames(データの列名を変更する)
  #ctrl
  sonaf_TRA_ctrl1<- rename(sonaf_TRA_ctrl1, "j.region" = .)
  sonaf_TRA_ctrl2<- rename(sonaf_TRA_ctrl2, "j.region" = .)
  sonaf_TRA_ctrl3<- rename(sonaf_TRA_ctrl3, "j.region" = .)
  sonaf_TRA_ctrl4<- rename(sonaf_TRA_ctrl4, "j.region" = .)
  sonaf_TRA_ctrl5<- rename(sonaf_TRA_ctrl5, "j.region" = .)
  sonaf_TRA_ctrl6<- rename(sonaf_TRA_ctrl6, "j.region" = .)
  #imm
  sonaf_TRA_imm1<- rename(sonaf_TRA_imm1, "j.region" = .)
  sonaf_TRA_imm2<- rename(sonaf_TRA_imm2, "j.region" = .)
  sonaf_TRA_imm3<- rename(sonaf_TRA_imm3, "j.region" = .)
  sonaf_TRA_imm4<- rename(sonaf_TRA_imm4, "j.region" = .)
  sonaf_TRA_imm5<- rename(sonaf_TRA_imm5, "j.region" = .)
  sonaf_TRA_imm6<- rename(sonaf_TRA_imm6, "j.region" = .)
  
  ##melting data for merge(データをマージさせるために形式を変更する)
  #ctrl
  TRA_CTRL1m <- melt(sonaf_TRA_ctrl1,id.vars = "j.region")
  TRA_CTRL2m <- melt(sonaf_TRA_ctrl2,id.vars = "j.region")
  TRA_CTRL3m <- melt(sonaf_TRA_ctrl3,id.vars = "j.region")
  TRA_CTRL4m <- melt(sonaf_TRA_ctrl4,id.vars = "j.region")
  TRA_CTRL5m <- melt(sonaf_TRA_ctrl5,id.vars = "j.region")
  TRA_CTRL6m <- melt(sonaf_TRA_ctrl6,id.vars = "j.region")
  #imm
  TRA_Imm1m <- melt(sonaf_TRA_imm1,id.vars = "j.region")
  TRA_Imm2m <- melt(sonaf_TRA_imm2,id.vars = "j.region")
  TRA_Imm3m <- melt(sonaf_TRA_imm3,id.vars = "j.region")
  TRA_Imm4m <- melt(sonaf_TRA_imm4,id.vars = "j.region")
  TRA_Imm5m <- melt(sonaf_TRA_imm5,id.vars = "j.region")
  TRA_Imm6m <- melt(sonaf_TRA_imm6,id.vars = "j.region")
  
  
  ##set working directory
  setwd("~/R/TCR_reparetoire_analysis/stress_TCR/plotfancyvj_4wks_TRA")
  
  ##import data sets
  #ctrl
  af4wks_TRA_ctrl1 <- read.table("filterd_cont1_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_ctrl2 <- read.table("filterd_cont2_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_ctrl3 <- read.table("filterd_cont3_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_ctrl4 <- read.table("filterd_cont4_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_ctrl5 <- read.table("filterd_cont5_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  #imm
  af4wks_TRA_imm1 <- read.table("filterd_imm1_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_imm2 <- read.table("filterd_imm2_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_imm3 <- read.table("filterd_imm3_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_imm4 <- read.table("filterd_imm4_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  af4wks_TRA_imm5 <- read.table("filterd_imm5_TRA.fancyvj.wt.txt", header = TRUE, sep="\t")
  
  ##rename_colnames(データの列名を変更する)
  #ctrl
  af4wks_TRA_ctrl1<- rename(af4wks_TRA_ctrl1, "j.region" = .)
  af4wks_TRA_ctrl2<- rename(af4wks_TRA_ctrl2, "j.region" = .)
  af4wks_TRA_ctrl3<- rename(af4wks_TRA_ctrl3, "j.region" = .)
  af4wks_TRA_ctrl4<- rename(af4wks_TRA_ctrl4, "j.region" = .)
  af4wks_TRA_ctrl5<- rename(af4wks_TRA_ctrl5, "j.region" = .)
  #imm
  af4wks_TRA_imm1<- rename(af4wks_TRA_imm1, "j.region" = .)
  af4wks_TRA_imm2<- rename(af4wks_TRA_imm2, "j.region" = .)
  af4wks_TRA_imm3<- rename(af4wks_TRA_imm3, "j.region" = .)
  af4wks_TRA_imm4<- rename(af4wks_TRA_imm4, "j.region" = .)
  af4wks_TRA_imm5<- rename(af4wks_TRA_imm5, "j.region" = .)
  
  ##melting data for merge(データをマージさせるために形式を変更する)
  #ctrl
  TRA_CTRL1m4 <- melt(af4wks_TRA_ctrl1,id.vars = "j.region")
  TRA_CTRL2m4 <- melt(af4wks_TRA_ctrl2,id.vars = "j.region")
  TRA_CTRL3m4 <- melt(af4wks_TRA_ctrl3,id.vars = "j.region")
  TRA_CTRL4m4 <- melt(af4wks_TRA_ctrl4,id.vars = "j.region")
  TRA_CTRL5m4 <- melt(af4wks_TRA_ctrl5,id.vars = "j.region")
  #imm
  TRA_Imm1m4 <- melt(af4wks_TRA_imm1,id.vars = "j.region")
  TRA_Imm2m4 <- melt(af4wks_TRA_imm2,id.vars = "j.region")
  TRA_Imm3m4 <- melt(af4wks_TRA_imm3,id.vars = "j.region")
  TRA_Imm4m4 <- melt(af4wks_TRA_imm4,id.vars = "j.region")
  TRA_Imm5m4 <- melt(af4wks_TRA_imm5,id.vars = "j.region")
  
  for (i in 1:6){
    x <- paste0("TRA_CTRL",i,"m")
    p <- paste0("TRA_CTRL",i,"m <- TRA_CTRL",i,"m%>%mutate(vjseg = paste(TRA_CTRL",i,"m[,1],TRA_CTRL",i,"m[,2])) ; 
              sonaf_C",i," <- TRA_CTRL",i,"m")
    eval(parse(text = p))
    
    y <- paste0("TRA_Imm",i,"m")
    p2 <- paste0("TRA_Imm",i,"m <- TRA_Imm",i,"m%>%mutate(vjseg = paste(TRA_Imm",i,"m[,1],TRA_Imm",i,"m[,2])) ; 
              sonaf_I",i," <- TRA_Imm",i,"m")
    eval(parse(text = p2))
  }
  
  for (i in 1:5){
    x <- paste0("TRA_CTRL",i,"m4")
    p <- paste0("TRA_CTRL",i,"m4 <- TRA_CTRL",i,"m4%>%mutate(vjseg = paste(TRA_CTRL",i,"m4[,1],TRA_CTRL",i,"m4[,2])) ; 
              wk4_C",i," <- TRA_CTRL",i,"m4")
    eval(parse(text = p))
    
    y <- paste0("TRA_Imm",i,"m4")
    p2 <- paste0("TRA_Imm",i,"m4 <- TRA_Imm",i,"m4%>%mutate(vjseg = paste(TRA_Imm",i,"m4[,1],TRA_Imm",i,"m4[,2])) ; 
              wk4_I",i,"<- TRA_Imm",i,"m4")
    eval(parse(text = p2))
  }
  
  
  ###jaccard matrix作成
  t0<- NULL
  t0<- matrix(nrow = 22,ncol = 22)
  
  ##1_[i,x]x固定でiが先に変わる
  ##1列目の1,2,3行から値が入る
  for (x in 1:6){
    for (i in 1:6){
      p <- paste0("t0[",i,",",x,"] <- jaccard(sonaf_C",x,"$vjseg,sonaf_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##2_[y(i),z(x)]x固定でiが先に変わる
  ##7列目の7,8,9行から値が入る
  for (x in 1:6){
    for (i in 1:6){
      y <- i+6
      z <- x+6
      p <- paste0("t0[",y,",",z,"] <- jaccard(sonaf_I",x,"$vjseg,sonaf_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##3_[i,z(x)]x固定でiが先に変わる
  ##7列目の1,2,3行から値が入る -> sonaf_I1 x sonaf_C1,2,3
  for (x in 1:6){
    for (i in 1:6){
      z <- x+6
      p <- paste0("t0[",i,",",z,"] <- jaccard(sonaf_I",x,"$vjseg,sonaf_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##4_[i,z(x)]x固定でiが先に変わる
  ##1列目の7,8,9行から値が入る -> sonaf_C1 x sonaf_I1,2,3
  for (x in 1:6){
    for (i in 1:6){
      y <- i+6
      p <- paste0("t0[",y,",",x,"] <- jaccard(sonaf_C",x,"$vjseg,sonaf_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  
  ##5_[i,x]x固定でiが先に変わる
  ##13列目の13,14,15行から値が入る
  for (x in 1:5){
    for (i in 1:5){
      y <- i+12
      z <- x+12
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_C",x,"$vjseg,wk4_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##[y(i),z(x)]x固定でiが先に変わる
  ##18列目の18,19,20行から値が入る
  for (x in 1:5){
    for (i in 1:5){
      y <- i+17
      z <- x+17
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_I",x,"$vjseg,wk4_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##[y(i),z(x)]x列固定でi行が先に変わる
  ###13列目の18,19,20行から値が入る -> wk4_C1 x wk4_I1,2,3
  for (x in 1:5){
    for (i in 1:5){
      y <- i+17
      z <- x+12
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_C",x,"$vjseg,wk4_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##[y(i),z(x)]x固定でiが先に変わる
  ##18列目の13,14,15行から値が入る -> wk4_I1 x wk4_C1,2,3
  for (x in 1:5){
    for (i in 1:5){
      y <- i+12
      z <- x+17
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_I",x,"$vjseg,wk4_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  
  
  ##[y(i),z(x)]x固定でiが先に変わる
  ##1列目の13,14,15行から値が入る -> wk4_I1 x wk4_C1,2,3
  for (x in 1:6){
    for (i in 1:5){
      y <- i+12
      p <- paste0("t0[",y,",",x,"] <- jaccard(sonaf_C",x,"$vjseg,wk4_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##[y(i),z(x)]x固定でiが先に変わる
  ##1列目の18,19,20行から値が入る -> wk4_I1 x wk4_C1,2,3
  for (x in 1:6){
    for (i in 1:5){
      y <- i+17
      p <- paste0("t0[",y,",",x,"] <- jaccard(sonaf_C",x,"$vjseg,wk4_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##[y(i),z(x)]x固定でiが先に変わる
  ##7列目の13,14,15行から値が入る -> wk4_I1 x wk4_C1,2,3
  for (x in 1:6){
    for (i in 1:5){
      y <- i+12
      z <- x+6
      p <- paste0("t0[",y,",",z,"] <- jaccard(sonaf_I",x,"$vjseg,wk4_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##[y(i),z(x)]x固定でiが先に変わる
  ##7列目の18,19,20行から値が入る -> wk4_I1 x wk4_C1,2,3
  for (x in 1:6){
    for (i in 1:5){
      y <- i+17
      z <- x+6
      p <- paste0("t0[",y,",",z,"] <- jaccard(sonaf_I",x,"$vjseg,wk4_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  
  
  
  
  
  ##[y(i),z(x)]x固定でiが先に変わる
  ##13列目の1,2,3行から値が入る -> wk4_C1 x sonaf_C1,2,3
  for (x in 1:5){
    for (i in 1:6){
      z <- x+12
      p <- paste0("t0[",i,",",z,"] <- jaccard(wk4_C",x,"$vjseg,sonaf_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##[y(i),z(x)]x固定でiが先に変わる
  ##13列目の7,8,9行から値が入る -> wk4_C1 x sonaf_I1,2,3
  for (x in 1:5){
    for (i in 1:6){
      y <- i+6
      z <- x+12
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_C",x,"$vjseg,sonaf_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  ##[y(i),z(x)]x固定でiが先に変わる
  ##18列目の1,2,3行から値が入る -> wk4_I1 x sonaf_C1,2,3
  for (x in 1:5){
    for (i in 1:6){
      z <- x+17
      p <- paste0("t0[",i,",",z,"] <- jaccard(wk4_I",x,"$vjseg,sonaf_C",i,"$vjseg)")
      eval(parse(text = p))
    }}
  ##[y(i),z(x)]x固定でiが先に変わる
  ##18列目の7,8,9行から値が入る -> wk4_C1 x sonaf_I1,2,3
  for (x in 1:5){
    for (i in 1:6){
      y <- i+6
      z <- x+17
      p <- paste0("t0[",y,",",z,"] <- jaccard(wk4_I",x,"$vjseg,sonaf_I",i,"$vjseg)")
      eval(parse(text = p))
    }}
  
  
  t <- as.matrix(t0)
  
  rownames(t) <- c("sonaf_C1","sonaf_C2","sonaf_C3","sonaf_C4","sonaf_C5","sonaf_C6",
                   "sonaf_I1","sonaf_I2","sonaf_I3","sonaf_I4","sonaf_I5","sonaf_I6",
                   "wk4_C1","wk4_C2","wk4_C3","wk4_C4","wk4_C5","wk4_I1","wk4_I2","wk4_I3","wk4_I4","wk4_I5")
  
  colnames(t) <- c("sonaf_C1","sonaf_C2","sonaf_C3","sonaf_C4","sonaf_C5","sonaf_C6",
                   "sonaf_I1","sonaf_I2","sonaf_I3","sonaf_I4","sonaf_I5","sonaf_I6",
                   "wk4_C1","wk4_C2","wk4_C3","wk4_C4","wk4_C5","wk4_I1","wk4_I2","wk4_I3","wk4_I4","wk4_I5")
  
  
  paletteLength <- 100
  myColor <- colorRampPalette(c("green", "white", "red"))(paletteLength)
  
  
  ComplexHeatmap::Heatmap(t,name="JCI",
                          cluster_columns = F, 
                          cluster_rows = F,
                          width = ncol(t)*unit(5, "mm"),
                          height = nrow(t)*unit(5, "mm"),
                          col = colorRamp2(breaks = c(0.9,0.95,1), colors = c("green", "white", "red")))
  
}