# Zhulab database Script
# Modified by Dai Yulong

## step1：根据.RData文件提取化合物名称和离子信息。
options(digits = 6)
parse_data <- function(data_path, name, outputpath) {
  # rm(list=ls())
  print(data_path)
  load(data_path)
  
  data <- get(name)

  # View(data)
  pos <- data$compound$pos
  neg <- data$compound$neg
  names_pos <- names(pos)
  names_neg <- names(neg)
  
  for (pos_name in names_pos) {
    # print(pos_name)
    if (as.numeric(length(pos[[pos_name]]$`30`))==0){
      next;
    }
    
    df_30 <- pos[[pos_name]]$`30`
    
    
    output_dir <-  paste0(outputpath, "/compound/pos/")
    if (!dir.exists(output_dir)){
      dir.create(output_dir,recursive = TRUE)
    }
    write.csv(df_30, file = paste0(outputpath, "/compound/pos/", "pos_", pos_name, ".csv"), row.names = FALSE)
  }
  
  for (neg_name in names_neg) {
    # print(pos_name)
    if (as.numeric(length(neg[[neg_name]]$`30`))==0){
      next;
      
    }
    
    df_30 <- neg[[neg_name]]$`30`
    output_dir <-  paste0(outputpath, "/compound/neg/")
    if (!dir.exists(output_dir)){
      dir.create(output_dir,recursive = TRUE)

    }
    write.csv(df_30, file = paste0(outputpath, "/compound/neg/", "neg_", neg_name, ".csv"), row.names = FALSE)
  }

  # meta process
  
  pos <- data$meta$pos$`30`
  neg <- data$meta$neg$`30`

  output_neg <-  paste0(outputpath, "/meta")
  if (!dir.exists(output_neg)){
    dir.create(output_neg,recursive = TRUE)
  }
  
  output_pos <-  paste0(outputpath, "/meta")
  if (!dir.exists(output_pos)){
    dir.create(output_pos,recursive = TRUE)
  }
  
  write.csv(neg, file = paste0(output_neg, "/neg_30.csv"), row.names = FALSE)
  write.csv(pos, file = paste0(output_pos, "/pos_30.csv"), row.names = FALSE)
  write.csv(data$meta$compound, file = paste0(outputpath, "/meta/compound", ".csv"), row.names = FALSE)
}

parse_data("E:/projectResearch/03Zhulab/zhulib/NISTlib.Rdata","NISTlib","E:/projectResearch/03Zhulab/output/NISTlib")
parse_data("E:/projectResearch/03Zhulab/zhulib/zhulib.Rdata", "zhulib", "E:/projectResearch/03Zhulab/output/zhulib")
parse_data("E:/projectResearch/03Zhulab/zhulib/zhuMetlib.Rdata","zhuMetlib","E:/projectResearch/03Zhulab/output/zhuMetlib")
parse_data("E:/projectResearch/03Zhulab/zhulib/zhuNISTlib.Rdata","zhuNISTlib","E:/projectResearch/03Zhulab/output/zhuNISTlib")



##step2：根据compond和meta结果，提取msp文件需要的"NAME","PRECURSORMZ","PRECURSORTYPE","FORMULA","IONMODE"等信息。

zhulab_data <- function(name,Ionmode,outpath){
  
  indir <- paste0(outpath,"/",name)  
  
  comp <- read.csv(paste0(indir,"/meta/compound.csv"),header = T,sep = ",")
  neg_30 <- read.csv(paste0(indir,"/meta/neg_30.csv"),header = T,sep = ",")
  pos_30 <- read.csv(paste0(indir,"/meta/pos_30.csv"),header = T,sep = ",")
  
  comp_neg <- merge(comp,neg_30[,c(1,4)],by = "labid")
  comp_pos <- merge(comp,pos_30[,c(1,4)],by = "labid")
  
  # pos msp file
  
  dfindex <- data.frame(x=NA,y=NA)
  
  if (Ionmode == 'Positive') {
    
    
    dat_all <- data.frame()
    for (i in 1:nrow(comp_pos)) {
      dat_block <- data.frame(lab=c("NAME:","PRECURSORMZ:","PRECURSORTYPE:","FORMULA:","IONMODE:","Num Peaks:"),val=0)
      
      dat_block$val[3] <- "[M+H]+"
      dat_block$val[5] <- Ionmode
      
      indat1 <- read.csv(paste0(indir,"/compound/pos/pos_",comp_pos[1][[1]][i],".csv"),header = T,sep = ",",check.names = F)
      indat1 <- indat1[,c(1,ncol(indat1))]
      dat_block$val[6] <- nrow(indat1)
      dat_block$val[1] <- comp_pos[2][[1]][i]
      dat_block$val[2] <- comp_pos[3][[1]][i]+1
      dat_block$val[4] <- comp_pos[4][[1]][i]
      dat_block_new <- as.data.frame(rbind(as.matrix(dat_block),as.matrix(indat1),as.matrix(dfindex)))
      
      dat_all <- rbind(dat_all,dat_block_new)
      dat_all[is.na(dat_all)] <- ""
    }
    write.table(dat_all,paste0(outpath,"/",name,"/pos_meta.txt"),row.names = F,col.names = F,quote = F,sep = " ")
    
  }else{
    
    # neg msp file 
    dat_all_neg <- data.frame()
    for (i in 1:nrow(comp_neg)) {
      dat_block_neg <- data.frame(lab=c("NAME:","PRECURSORMZ:","PRECURSORTYPE:","FORMULA:","IONMODE:","Num Peaks:"),val=0)
      
      dat_block_neg$val[3] <- "[M-H]-"
      dat_block_neg$val[5] <- Ionmode
      
      indat1_neg <- read.csv(paste0(indir,"/compound/neg/neg_",comp_neg[1][[1]][i],".csv"),header = T,sep = ",",check.names = F)
      indat1_neg <- indat1_neg[,c(1,ncol(indat1_neg))]  # 有些物质的强度表中有2列intensity值，只取其中1个intensity值。
      dat_block_neg$val[6] <- nrow(indat1_neg)
      dat_block_neg$val[1] <- comp_neg[2][[1]][i]
      dat_block_neg$val[2] <- comp_neg[3][[1]][i]-1
      dat_block_neg$val[4] <- comp_neg[4][[1]][i]
      dat_block_neg_new <- as.data.frame(rbind(as.matrix(dat_block_neg),as.matrix(indat1_neg),as.matrix(dfindex)))
      
      dat_all_neg <- rbind(dat_all_neg,dat_block_neg_new)
      dat_all_neg[is.na(dat_all_neg)] <- ""
    }
    write.table(dat_all_neg,paste0(outpath,"/",name,"/neg_meta.txt"),row.names = F,col.names = F,quote = F,sep = " ")
    
  }
  
}

zhulab_data("NISTlib","Positive","E:/projectResearch/03Zhulab/output")
zhulab_data("NISTlib","Negative","E:/projectResearch/03Zhulab/output")
zhulab_data("zhulib","Positive","E:/projectResearch/03Zhulab/output")
zhulab_data("zhulib","Negative","E:/projectResearch/03Zhulab/output")
zhulab_data("zhuMetlib","Positive","E:/projectResearch/03Zhulab/output")
zhulab_data("zhuMetlib","Negative","E:/projectResearch/03Zhulab/output")
zhulab_data("zhuNISTlib","Positive","E:/projectResearch/03Zhulab/output")
zhulab_data("zhuNISTlib","Negative","E:/projectResearch/03Zhulab/output")



##step3： 将4个数据库分别整理好的neg和pos模式的数据库文件按照离子模式进行数据合并
outpath <- "E:/projectResearch/03Zhulab/output"

for (mod in c("pos","neg")) {
  
  dat <- data.frame()
  for (labname in c("NISTlib","zhulib","zhuMetlib","zhuNISTlib")) {
    infile <- paste0(outpath,"/",labname,"/",mod,"_meta.txt")
    df1 <- read.delim(infile,check.names = F,header = F)
    dat <- rbind(dat,df1)
  }
  write.table(dat,"E:/projectResearch/03Zhulab/output/",mod,"_all_meta.txt",row.names = F,col.names = F,quote = F) 
}

