##function(x,y)
##x = CTRL file number 
##y = Imm file number 
#library(ggplot2)
#library(reshape2)
#library( RColorBrewer )
#library( devtools )
#library(circlize)
#library(dplyr)

plotfilteredTCR <- function(x, y) {

setwd("~/R/TCR_reparetoire_analysis/stress_TCR/sonaf_filtered_TRA")
for (a in 1:x){
  p <- paste0("filC",a," <- read.table('filterd_CTRL",a,"_TRA.txt', header = TRUE, sep='\t')")
  eval(parse(text = p))
}
for (b in 1:y){
  p <- paste0("filI",b," <- read.table('filterd_Imm",b,"_TRA.txt', header = TRUE, sep='\t')")
  eval(parse(text = p))
}

setwd("~/R/TCR_reparetoire_analysis/stress_TCR/sonaf_TRA")
for (a in 1:x){
  p <- paste0("C",a," <- read.table('CTRL",a,"_TRA.txt', header = TRUE, sep='\t')")
  eval(parse(text = p))
}
for (b in 1:y){
  p <- paste0("I",b," <- read.table('Imm",b,"_TRA.txt', header = TRUE, sep='\t')")
  eval(parse(text = p))
}


t0 <-NULL
t0 <-as.data.frame(t0)
for (a in 1:x){
  p <- paste0("t0[",a,",1] <- nrow(C",a,")%>%as.character()")
  eval(parse(text = p))
}
for (a in 1:x){
  p <- paste0("t0[",a,",2] <- nrow(filC",a,")%>%as.character()")
  eval(parse(text = p))
}

w <- (x+1)
z <- (x+y)
for (a in w:z){
  b <- a-x
  p <- paste0("t0[",a,",1] <- nrow(I",b,")%>%as.character()")
  eval(parse(text = p))
}
for (a in w:z){
  b <- a-6
  p <- paste0("t0[",a,",2] <- nrow(filI",b,")%>%as.character()")
  eval(parse(text = p))
}
t0$fil <- as.numeric(t0$V1)-as.numeric(t0$V2)
name1 <- paste0("C",1:x)
name2 <- paste0("I",1:y)
row.names(t0) <- c(name1,name2)
colnames(t0) <- c("nonfil","fil","deb")
t0$sample <- rep(c("CTRL", "Imm"), each = x)

ggplot(t0)+geom_point(aes(x = sample, y =deb,fill = sample,color =sample))
}