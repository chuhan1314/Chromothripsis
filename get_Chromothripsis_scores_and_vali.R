.Library="~/library"
.libPaths("~/library")
library(plyr)
library(dplyr)
library(hdp)
library(stringr)
library('GenomicFeatures')
library(BSgenome.Hsapiens.UCSC.hg38)
library(IRanges)
library("GenomicRanges")
library(RColorBrewer)
library(mclust)
library(pheatmap)
library('pROC')
set.seed(1)


snps = read.table("hg38.chrom.sizes.txt", header=F, sep=",")
levels(snps$V2) <- c(levels(snps$V2), "chr23","chr24")
snps$V2[snps$V2 == "chrX"] <- "chr23"
snps$V2[snps$V2 == "chrY"] <- "chr24"
chr_list <- as.character(1:24) 

cyto <- list()
for  (i in (1:24)){
  #i=1
  chr_num = i
  chr_filter <- paste("chr",chr_num,sep="")
  limit <- snps[snps$V2 == chr_filter,]
  post <- limit$V3
  x <- seq(from = 0, to = (post - 1000000), by = 1000000)
  p <- seq(from = 1000000, to = post, by = 1000000)
  alfa<- length(p)
  tail <-  c(p[alfa],post)
  y <- data.frame(x,p)
  y <- rbind(y,tail)
  y$chr <- rep(chr_filter,nrow(y))
  colnames(y)<-c("start","end","chr")  
  cyto[[chr_filter]] <- y
}
caryo <- do.call("rbind", cyto) # merge results from all chromosomes  
caryo$chr<-gsub("chr", "",caryo$chr)

gr1 = with(caryo, GRanges(chr, IRanges(start=start, end=end)))


cyto_10mb <- list()
for  (i in c(1:22,"X")){
  chr_filter = i
  limit <- snps[snps$V1 == chr_filter,]
  post <- limit$V3
  x <- seq(from = 0, to = (post - 10000000), by = 10000000)
  p <- seq(from = 10000000, to = post, by = 10000000)
  alfa<- length(p)
  tail <-  c(p[alfa],post)
  y <- data.frame(x,p)
  y <- rbind(y,tail)
  y$chr <- rep(chr_filter,nrow(y))
  colnames(y)<-c("start","end","chr")  
  cyto_10mb[[chr_filter]] <- y
}
cyto_10mb_all <- do.call("rbind", cyto_10mb)
gr1_10mb = with(cyto_10mb_all, GRanges(chr, IRanges(start=start, end=end)))



type=list.files("./TCGA/")
cnv_mm=list()
for (i in type) {
  tmp=read.csv(paste0("./TCGA/",i,"/",i,"-ALLele-specificCopyNumberSegment.csv"),row.names = 1)
  tmp$type=i
  cnv_mm[[i]]=tmp
}
#全部都是somatic cnv，用ascat2，肿瘤和血液比较获得的结果。
cnv_mm=do.call(rbind,cnv_mm)
#length(unique(cnv_mm$Sample))
cnv_mm$sample=substr(cnv_mm$Sample,1,12)
cnv_mmrf_sel=cnv_mm
length(unique(cnv_mmrf_sel$sample))

cnv_mmrf_sel$Chrom=sapply(strsplit(cnv_mmrf_sel$Chromosome,split="hr"),function(x){return(x[2])})
###该算法中，major就是totol cnv
colnames(cnv_mmrf_sel)= c("GDC_Aliquot","Chromosome","start","end","major",
                          "Major_Copy_Number","minor","Sample","type","sample","Chrom")  
#cnv_mmrf_sel=cnv_mmrf_sel[!cnv_mmrf_sel$sample %in% c("TCGA-FE-A237","TCGA-EL-A4JV","TCGA-DJ-A3UT"),]
#cnv_mmrf_sel$Chrom=as.numeric(cnv_mmrf_sel$Chrom)
cnv_mmrf_sel$start=as.numeric(cnv_mmrf_sel$start)
cnv_mmrf_sel$end=as.numeric(cnv_mmrf_sel$end)
cnv_mmrf_sel$minor=as.numeric(cnv_mmrf_sel$minor)
cnv_mmrf_sel$major=as.numeric(cnv_mmrf_sel$major)
#cnv_mmrf_sel_bak=cnv_mmrf_sel
#cnv_mmrf_sel$major=cnv_mmrf_sel$major+cnv_mmrf_sel$minor
phe = read.table("./tcga_chromth_phe.txt")
phe$sample = gsub("[.]","-",phe$sample)
cnv_mmrf_sel = cnv_mmrf_sel[cnv_mmrf_sel$sample %in% phe$sample, ]
length(unique(cnv_mmrf_sel$sample))


#########################################################################
#########################################################################
#########################################################################
#cnv_mmrf_sel=cnv_mmrf_sel[grep("^p",cnv_mmrf_sel$sample),]#zhao
#cnv_mmrf_sel=cnv_mmrf_sel_bak[grep("^Pt",cnv_mmrf_sel_bak$sample),]#hugo
#########################################################################
#########################################################################
#########################################################################


cnv_mmrf_sel$code_row<- 1:nrow(cnv_mmrf_sel)
# Remove IgH, IgK and IgL loci to reduce artefacts due to VDJ-rearrangement / class-switch recombination
igh_cnv <- cnv_mmrf_sel[cnv_mmrf_sel$Chrom == 14 & cnv_mmrf_sel$start >106032614 &  cnv_mmrf_sel$end< 108288051 | 
                          cnv_mmrf_sel$Chrom == 22 & cnv_mmrf_sel$start >21080474. &  cnv_mmrf_sel$end< 26065085 | 
                          cnv_mmrf_sel$Chrom == 2 & cnv_mmrf_sel$start >87090568 &  cnv_mmrf_sel$end< 93274235,]
cnv_mmrf_no_igh<- cnv_mmrf_sel[!cnv_mmrf_sel$code_row %in% igh_cnv$code_row ,]

# Remove chromosome X to avoid overestimation of deletion
cnv_mmrf_no_igh_no_x<- cnv_mmrf_no_igh[cnv_mmrf_no_igh$Chrom!="X",]

# Remove segments < 50Mb in size
cnv_mmrf_final<- cnv_mmrf_no_igh_no_x[(cnv_mmrf_no_igh_no_x$end-cnv_mmrf_no_igh_no_x$start)>50000,]

# Collapse adjacent segments with the same copy number
cnv_mmrf2<- list()
sample_list<- unique(cnv_mmrf_final$sample)
for(j in (1:length(sample_list))){
  #j=1
  cna_mmrf_single<- cnv_mmrf_final[cnv_mmrf_final$sample == sample_list[j],]
  sam_cnv_list<- list()
  chr_list<- unique(cna_mmrf_single$Chrom)
  for(i in (1:length(chr_list))){
    #i=1
    cna_mmrf_single_chr<- cna_mmrf_single[cna_mmrf_single$Chrom == chr_list[i],]
    
    list_chr<- list()
    vec<- rle((paste(cna_mmrf_single_chr$major)))$length 
    for(w in (1:length(vec))){
      if(w==1){
        int<- cna_mmrf_single_chr[1:vec[w],]
        cna_mmrf_single_row<- c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
      }else{
        int<- cna_mmrf_single_chr[(sum(vec[1:(w-1)])+1):sum(vec[1:(w)]),]
        cna_mmrf_single_row<- c(int$sample[1], int$Chrom[1], int$start[1], int$end[nrow(int)], int$major[1], max(int$minor))
      }
      list_chr[[w]]<- cna_mmrf_single_row
    }
    list_chr2<- do.call("rbind",list_chr)
    sam_cnv_list[[i]]<- list_chr2
  }
  sam_cnv_list2<- do.call("rbind", sam_cnv_list)
  cnv_mmrf2[[j]]<-  sam_cnv_list2
}
cnv_mmrf<- do.call("rbind", cnv_mmrf2)
cnv_mmrf<- as.data.frame.matrix(cnv_mmrf)
colnames(cnv_mmrf)<-c("sample","Chrom","start","end", "major","minor")
cnv_mmrf$sample<- as.character(as.character(cnv_mmrf$sample))
cnv_mmrf$Chrom<- as.character(as.character(cnv_mmrf$Chrom))
cnv_mmrf[,3:ncol(cnv_mmrf)]<- apply(cnv_mmrf[,3:ncol(cnv_mmrf)], 2, function(x){as.numeric(as.character(x))})

mat_sig_cnv_final<- list()
max_10mb<- 100
max_copy_number<- max(cnv_mmrf$major)

###############样本修剪###############
cnv_mmrf = cnv_mmrf[cnv_mmrf$sample != sample_list[w],]
length(unique(cnv_mmrf$sample))
#### Create list for each copy number feature
count_10mb_all<-list()
size_all<- list()
count_jump_all<- list()
count_cnv_all<- list()
band_rate_all<-list()
osci_all<- list()
sample_list<- unique(cnv_mmrf$sample)

# Start loop for each sample
for(w in (1:length(sample_list))){
  #w=1
  cnv_mmrf2<- cnv_mmrf[cnv_mmrf$sample ==sample_list[w],]
  
  # Segment size
  cnv_mmrf2$seg_size<- (cnv_mmrf2$end - cnv_mmrf2$start) 
  size_all[[w]]<-cnv_mmrf2[,c("sample","seg_size")]
  
  # Number of breaks per 10Mb segment
  cnv_temp_brk<- cnv_mmrf2[,c(1,2,3,5,6,7)]
  cnv_mmrf2_second<- cnv_mmrf2[,c(1,2,4,5,6,7)]
  
  # Remove diploid whole chromosome regions
  int_dipl<- as.data.frame.matrix(table(cnv_mmrf2_second$Chrom, cnv_mmrf2_second$major))
  diploid_chr<- rownames(int_dipl[int_dipl$`2`==1 & rowSums(int_dipl)==1,])
  
  # Remove second, so each break is only counted once
  cnv_temp_brk<- cnv_mmrf2_second[! cnv_mmrf2_second$Chrom %in% diploid_chr,]
  
  gr_cna_comm = with(cnv_temp_brk, GRanges(Chrom, IRanges(start=end, end=end)))
  values(gr_cna_comm) <- DataFrame(sample = cnv_temp_brk$sample, major = cnv_temp_brk$major, minor= cnv_temp_brk$minor, seg_size= cnv_temp_brk$seg_size)
  range_10mb <- merge(as.data.frame(gr_cna_comm),as.data.frame(gr1_10mb),by="seqnames",suffixes=c("A","B"))
  range_dri_10mb <- range_10mb[with(range_10mb, startB <= startA & endB >= endA),]
  count_brk_10mb<- as.data.frame(table(paste(range_dri_10mb$seqnames, range_dri_10mb$startB)))
  
  # Use the reference for 10Mb without CNV breaks
  options("scipen"=100, digits=10) 
  cyto_10mb_all$Var1<- paste(cyto_10mb_all$chr, cyto_10mb_all$start)
  count_10mb_file <-join(cyto_10mb_all, count_brk_10mb, by="Var1")
  count_10mb_file[is.na(count_10mb_file)]<-0
  count_10mb_file$sample<- sample_list[w]
  count_10mb_df<- count_10mb_file[,c("sample","Freq")]
  colnames(count_10mb_df)[2]<-"count"
  count_10mb_all[[w]]<-count_10mb_df 
  
  # Assess adjacent change in copy number (jump)
  cnv_mmrf2_second$jump<- NA
  summary_jump<- as.data.frame(table(cnv_mmrf2_second$Chrom))
  cnv_mmrf2_second_jump<- cnv_mmrf2_second
  if(length(unique(cnv_mmrf2_second_jump$Chrom))!=0){
    chr_list_jump<- unique(cnv_mmrf2_second_jump$Chrom)
    cnv_mmrf2_second_jump$jump<- NA
    all_chr_jump<- list()
    for(jj in chr_list_jump){
      #jj="1"
      cnv_mmrf2_second_jump_int<- cnv_mmrf2_second_jump[cnv_mmrf2_second_jump$Chrom == jj,]
      for(z in (1:nrow(cnv_mmrf2_second_jump_int))){
        if(z==1){
          #z=1
          cnv_mmrf2_second_jump_int$jump[1]=NA
        }else{
          cnv_mmrf2_second_jump_int$jump[z]<- abs((cnv_mmrf2_second_jump_int$major[z]) - 
                                                    (cnv_mmrf2_second_jump_int$major[z-1]))
        }
      }
      all_chr_jump[[jj]]<-cnv_mmrf2_second_jump_int
    }
    all_chr_jump2<- do.call("rbind", all_chr_jump)
  }else{
    all_chr_jump2<- cnv_mmrf2_second[1,]
    all_chr_jump2$jump<-0
  }
  all_chr_jump2<- all_chr_jump2[! is.na(all_chr_jump2$jump),]
  temp_jump<- all_chr_jump2[,c("sample","jump")]
  colnames(temp_jump)<- c("sample","count")
  count_jump_all[[w]]<-temp_jump 
  
  # Count absolute copy number of each segment
  count_cnv_final_df<- cnv_mmrf2[,c("sample","major")]
  colnames(count_cnv_final_df)[2]<-"count"
  count_cnv_all[[w]]<- count_cnv_final_df[,c("sample","count")]
  
  # Count breakpoints per chromosome arm
  chrom_arms<- read.delim("CentromerePosition_hg38.txt")##################################################
  chrom_arms$chrom<- gsub("chr", "",chrom_arms$chrom)
  cnv_mmrf2_second<- cnv_mmrf2[,c(1,2,4,5,6,7)]
  cnv_temp_brk_arm <- cnv_mmrf2_second
  
  if(nrow(cnv_temp_brk_arm)!=0){
    gr_cna_comm = with(cnv_temp_brk_arm, GRanges(Chrom, IRanges(start=end, end=end)))
    values(gr_cna_comm) <- DataFrame(sample = cnv_temp_brk_arm$sample, 
                                     major = cnv_temp_brk_arm$major, minor= cnv_temp_brk_arm$minor, seg_size= cnv_temp_brk_arm$seg_size)
    
    gr_band = with(chrom_arms, GRanges(chrom, IRanges(start=chromStart, end=chromEnd)))
    
    range_arm <- merge(as.data.frame(gr_cna_comm),as.data.frame(gr_band),by="seqnames",suffixes=c("A","B"))
    range_arm$arm<- NA
    range_arm$arm[range_arm$startA> range_arm$endB]<-"q_arm"
    range_arm$arm[range_arm$startA< range_arm$startB]<-"p_arm"
    range_arm$arm[range_arm$startB <= range_arm$startA & range_arm$endB >= range_arm$startA]<-"centro"
    
    table(paste(range_arm$seqnames, range_arm$arm))
    db_arm_counts<- as.data.frame(table(paste(range_arm$seqnames, range_arm$arm)))
  }else{
    db_arm_counts<- matrix(c("13 q_arm", 0), nrow=1)
    db_arm_counts<- as.data.frame(db_arm_counts)
    colnames(db_arm_counts)<- c("Var1","Freq")
  }
  
  file_int_band<- as.data.frame(c(paste0(c(1:22), (" p_arm")), paste0(c(1:22), (" q_arm"))))
  colnames(file_int_band)<-"Var1"
  band_rate<- join(file_int_band, db_arm_counts, by="Var1")
  band_rate[is.na(band_rate)]<-0
  band_rate$sample<- sample_list[w]
  colnames(band_rate)[2]<-"count"
  band_rate_all[[w]]<-band_rate[,c("sample","count")]
  
  # Assess oscillating copy number length 
  out<-c()
  chrs<-unique(cnv_mmrf2$Chrom)
  cnv_mmrf2$tot<- cnv_mmrf2$major
  oscCounts<-c()
  for(c in chrs){
    #c="3"
    currseg<-cnv_mmrf2[cnv_mmrf2$Chrom==c,"tot"]
    currseg<-round(as.numeric(currseg))
    #currseg
    
    if(length(currseg)>3){
      prevval<-currseg[1]
      count=0
      for(j in 3:length(currseg)){
        if(j==length(currseg)){
          oscCounts<-rbind(oscCounts,c(c,count))
          count=0
        }else{
          if(abs(currseg[j]-prevval)<=1 & currseg[j]!=currseg[j-1]){
            count<-count+1
          }else{
            oscCounts<-rbind(oscCounts,c(c,count))
            count=0
          }
        }
        prevval<-currseg[j-1]
      }
    }else{
      oscCounts<- rbind(oscCounts,c(c,0))
    }
  }
  
  oscCounts_df<- as.data.frame(oscCounts)
  oscCounts_df$sample<- sample_list[w]
  oscCounts_df$V2<-as.numeric(as.character(oscCounts_df$V2))
  osci_all[[w]]<-oscCounts_df
}


set.seed(999)
# mclust for segment size classification
size_all22<-do.call("rbind", size_all)
size_all22<- as.data.frame(size_all22)
size_all22$seg_size<- as.numeric(as.character(size_all22$seg_size))
size_all22_alt<- size_all22
myMclust_size <- Mclust(size_all22_alt$seg_size,G=2:10,verbose=FALSE)
size_all22_alt$size_code <- myMclust_size$classification
size_coordinate<- aggregate(seg_size ~ size_code, data = size_all22_alt, max)

# define absolute CNV state
count_cnv_all2<-do.call("rbind", count_cnv_all)
count_cnv_all2$cnv_count_code<- NA
count_cnv_all2$cnv_count_code[count_cnv_all2$count==0]<- 1
count_cnv_all2$cnv_count_code[count_cnv_all2$count==1]<- 2
count_cnv_all2$cnv_count_code[count_cnv_all2$count==2]<- 3
count_cnv_all2$cnv_count_code[count_cnv_all2$count==3]<- 4
count_cnv_all2$cnv_count_code[count_cnv_all2$count>=4]<- 5
count_cnv_coordinate<- aggregate(count ~ cnv_count_code, data = count_cnv_all2, max)

# mclust for CNV breakpoints per 10Mb
count_10mb_all2<-do.call("rbind", count_10mb_all)
count_10mb_all2<- as.data.frame(count_10mb_all2)
count_10mb_all2$count<- as.numeric(as.character(count_10mb_all2$count))
count_10mb_all2_alt<- count_10mb_all2[count_10mb_all2$count!=0,]
myMclust_10mb <- Mclust(count_10mb_all2_alt$count,G=4:4,verbose=FALSE)#######################raw4:4################？
count_10mb_all2_alt$mb_code <- myMclust_10mb$classification
count_cnv_coordinate_10mb <- aggregate(count ~ mb_code, data = count_10mb_all2_alt, max)

# mclust for jumps between adjacent segments 
count_jump_all2<-do.call("rbind", count_jump_all)
count_jump_all2<- as.data.frame(count_jump_all2)
count_jump_all2$count<- as.numeric(as.character(count_jump_all2$count))
count_jump_all2_alt<- count_jump_all2
myMclust_jump<- Mclust(count_jump_all2_alt$count,G=2:3,verbose=FALSE)
count_jump_all2_alt$cnv_count_code <- myMclust_jump$classification
jump_coordinate_10mb <- aggregate(count ~ cnv_count_code, data = count_jump_all2_alt, max)

# mclust for breakpoints per chromosome arm
band_rate_all2<-do.call("rbind", band_rate_all)
band_rate_all2<- as.data.frame(band_rate_all2)
band_rate_all2$Freq<- as.numeric(as.character(band_rate_all2$count))
band_rate_all2_alt<- band_rate_all2[band_rate_all2$count!=0,]
myMclust_band<- Mclust(band_rate_all2_alt$count,G=5,verbose=FALSE)
band_rate_all2_alt$band_code <- myMclust_band$classification
band_coordinate <- aggregate(count ~ band_code, data = band_rate_all2_alt, max)

# mclust for oscillation
osci_all2<-do.call("rbind", osci_all)
osci_all2<- as.data.frame(osci_all2[,c(3,2)])
colnames(osci_all2)[2]<-"count"
osci_all2$count<- as.numeric(as.character(osci_all2$count))
osci_all2_alt<- osci_all2
myMclust_osci<- Mclust(osci_all2_alt$count,G=4:5,verbose=FALSE) 
osci_all2_alt$osci_code <- myMclust_osci$classification
osci_coordinate <- aggregate(count ~ osci_code, data = osci_all2_alt, max)

count_cnv_coordinate_10mb$type<- "10_MB"
colnames(count_cnv_coordinate_10mb)<-c("code","count_limit","type")
count_cnv_coordinate$type<-"cnv_count"
colnames(count_cnv_coordinate)<-c("code","count_limit","type")
jump_coordinate_10mb$type<- "jump"
colnames(jump_coordinate_10mb)<-c("code","count_limit","type")
band_coordinate$type<-"band"
colnames(band_coordinate)<-c("code","count_limit","type")
osci_coordinate$type<-"osci"
colnames(osci_coordinate)<-c("code","count_limit","type")
size_coordinate$type<-"size"
colnames(size_coordinate)<-c("code","count_limit","type")

final_classification_value<- rbind.data.frame(count_cnv_coordinate_10mb,
                                              count_cnv_coordinate,
                                              jump_coordinate_10mb,
                                              band_coordinate,
                                              osci_coordinate,
                                              size_coordinate
)


osci_tab<- as.data.frame.matrix(table(osci_all2_alt$sample, osci_all2_alt$osci_code))
colnames(osci_tab)<- paste0("osci_", colnames(osci_tab))
osci_tab$sample<- rownames(osci_tab)

band_tab<- as.data.frame.matrix(table(band_rate_all2_alt$sample, band_rate_all2_alt$band_code))
colnames(band_tab)<- paste0("band_", colnames(band_tab))
band_tab$sample<- rownames(band_tab)

jump_tab<- as.data.frame.matrix(table(count_jump_all2_alt$sample, count_jump_all2_alt$cnv_count_code))
colnames(jump_tab)<- paste0("jump_", colnames(jump_tab))
jump_tab$sample<- rownames(jump_tab)

# Jump category has a small number of missing values- need to input 0 
miss_jump_sam<- osci_tab$sample[!osci_tab$sample %in% jump_tab$sample]
##################################################
############需要修改列数，与jump匹配##############
##################################################
miss_jmp_file<- cbind(rep(0, length(miss_jump_sam)), rep(0, length(miss_jump_sam)), 
                      miss_jump_sam)
##################################################
miss_jmp_file<- as.data.frame(miss_jmp_file)
colnames(miss_jmp_file)<-colnames(jump_tab)
jump_tab<- rbind.data.frame(jump_tab, miss_jmp_file)

mb_10_tab<- as.data.frame.matrix(table(count_10mb_all2_alt$sample, count_10mb_all2_alt$mb_code))
colnames(mb_10_tab)<- paste0("mb_10_", colnames(mb_10_tab))
mb_10_tab$sample<- rownames(mb_10_tab)

count_tab<- as.data.frame.matrix(table(count_cnv_all2$sample, count_cnv_all2$cnv_count_code))
colnames(count_tab)<- paste0("count_cnv_", colnames(count_tab))
count_tab$sample<- rownames(count_tab)

size_tab<- as.data.frame.matrix(table(size_all22_alt$sample, size_all22_alt$size_code))
colnames(size_tab)<- paste0("size_cnv_", colnames(size_tab))
size_tab$sample<- rownames(size_tab)

hdp_final<-Reduce(merge, list(mb_10_tab, count_tab, jump_tab, 
                              band_tab, osci_tab, size_tab))

#####

head(hdp_final)
rownames(hdp_final)<- hdp_final$sample
hdp_final2<- hdp_final[,-1] # remnove column "sample"
hdp_final2[,1:ncol(hdp_final2)]<- apply(hdp_final2[,1:ncol(hdp_final2)], 2, function(x){as.numeric(as.character(x))})
pheatmap(t(hdp_final2), show_colnames = FALSE)

mat_final<- hdp_final2
n=30
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(n)

cnv_colors<- cols[c(1,5,10, 13,17,24)]
cnv_colors_final<- c(rep(cnv_colors[1], 3), rep(cnv_colors[2], 5), rep(cnv_colors[3], 2), 
                     rep(cnv_colors[4], 2), rep(cnv_colors[5], 3), rep(cnv_colors[6], 9))

channel_names2 <- c("1", "2","3","1", "2", "3", "4", "5", "1", "2", "1", "2", "1", "2", "3", "1", "2", "3", "4", "5", "6", "7", "8","9")

x <-barplot(colSums(mat_final), las=2,
            col=cnv_colors_final, border = NA, xaxt = "n", cex.axis = 1.15)
axis(1, at=x,  label=rep("",26), mgp= c(3,1,0.2))
mtext(1, at=x, text=c(channel_names2), col=cnv_colors_final, padj = 1.5, cex = 1.15)

names_pts<- rownames(mat_final)
channel_names<- colnames(mat_final)
genomicData<- (mat_final) 
n<- ncol(genomicData)
shape<- 1
invscale<- 1
hdp<- hdp_init(ppindex=0, #index of the parent DP for initial DP
               cpindex=1, #index of alphaa and alphab for initial DP
               hh=rep(1/n,n), #params for base distn (uniform Dirichlet)
               alphaa=shape,
               alphab=invscale)

hdp<- hdp_adddp(hdp,
                numdp=nrow(genomicData),
                pp=1,
                cp=1)

hdp<- hdp_setdata(hdp= hdp,dpindex=1:nrow(genomicData)+1,data=genomicData)
hdp<- dp_activate(hdp,1:(nrow(genomicData)+1),10)

chlist <- vector("list", 4)
for (i in 1:4){
  chlist[[i]] <- hdp_posterior(hdp,
                               burnin=20000,
                               n=100,
                               space=50,
                               cpiter=3,
                               seed=i*1e4)
}

mut_example_multi <- hdp_multi_chain(chlist)
saveRDS(mut_example_multi, "tcga_ctlpS9200_to_mut_example_multi_4_chain.RDS")


mut_example_multi=readRDS("tcga_ctlpS9200_to_mut_example_multi_4_chain.RDS")###################read  data#####################
mut_example_multi_0.85_10 <- hdp_extract_components(mut_example_multi, cos.merge = 0.9, min.sample = 10) 
mut_example_multi <- mut_example_multi_0.85_10
par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
plot_comp_size(mut_example_multi, bty="L")

mut_example_multi_plot <- mut_example_multi
posteriorMeans_plot<- t(comp_categ_distn(mut_example_multi_plot)[[1]])
rownames(posteriorMeans_plot)<-  channel_names2
plotnames<- c("offset", "CN-SIG1", "CN-SIG2", "CN-SIG3", "CN-SIG4", "CN-SIG5")

for  (i in (2:6)){
  x <- barplot(posteriorMeans_plot[,i], las=2,
               col=cnv_colors_final, border = NA, xaxt = "n",  cex.axis = 1.5, main= plotnames[i], ylim = c(0, 0.25), cex.main = 2)
  axis(1, at=x,  label=rep("",24), mgp= c(3,1,0.2))
  mtext(1, at=x, text=c(channel_names2), col=cnv_colors_final, padj = 1.5, cex = 1.5)
}

x<-((mut_example_multi@comp_dp_distn))
kk<- x[["mean"]]
rownames(kk)<-c("offset", rownames(genomicData))
colnames(kk)<-c("offset", "CN-SIG1", "CN-SIG2", "CN-SIG3", "CN-SIG4", "CN-SIG5")
head(kk)
kk = kk[-1,-1] 
kk =as.data.frame(kk)
kk$sample=rownames(kk)
kk1=merge(kk,phe,by.y="sample",by.x="sample")
kk1$chromothripsis_code=ifelse(kk1$ctlp=="chromth",1,0)
kk2=kk1[,c(colnames(kk),"chromothripsis_code")]
colnames(kk2)=c("CN_SIG1","CN_SIG2","CN_SIG3","CN_SIG4","CN_SIG5","sample","chromothripsis_code")


rownames(kk2)=kk2$sample
kk2=kk2[,colnames(kk2)!="sample"]
datafr_genomes=kk2

AUCV = NULL
ROCVSENS = NULL
ROCVSPENS = NULL
RUNV =NULL

Len = length(datafr_genomes[,1]) 
SS = 1:Len
start=1
end= 930 # if the dataset contains 752 samples, then 10x cross-fold validation requires 75 for this value
num_add <- 923 # and again here

for (i in 1:10){ 
  i=1
  kk = SS[-(start:end)]
  start=start + num_add 
  end=end + num_add
  
  aacov1 = data.frame(chromth=datafr_genomes$chromothripsis_code[kk],
                      CN_SIG1 = (datafr_genomes$CN_SIG1)[kk],
                      CN_SIG2 = (datafr_genomes$CN_SIG2)[kk],
                      CN_SIG3 = (datafr_genomes$CN_SIG3)[kk],
                      CN_SIG4 = (datafr_genomes$CN_SIG4)[kk],
                      CN_SIG5 = (datafr_genomes$CN_SIG5)[kk]
  )
  
  aacov = glm(chromth~.,data = aacov1,family='binomial')
  
  aacov2 = data.frame(chromth=datafr_genomes$chromothripsis_code[-kk],
                      CN_SIG1 = (datafr_genomes$CN_SIG1)[-kk],
                      CN_SIG2 = (datafr_genomes$CN_SIG2)[-kk],
                      CN_SIG3 = (datafr_genomes$CN_SIG3)[-kk],
                      CN_SIG4 = (datafr_genomes$CN_SIG4)[-kk],
                      CN_SIG5 = (datafr_genomes$CN_SIG5)[-kk]
                      
  )
  
  predpr <- predict(aacov,newdata=aacov2,type=c("response"))
  roccurve <- roc(aacov2$chromth,predpr)
  aa = auc(roccurve)  
  aa = as.numeric(aa) 
  
  AUCV = c(AUCV,aa)
  ROCVSENS = c(ROCVSENS,roccurve$sensitivities)
  ROCVSPENS = c(ROCVSPENS,roccurve$specificities)
  RUNV =c(RUNV,rep(i,length(roccurve$specificities)))
}

aa<- mean(AUCV)
aa

Y = ROCVSENS
X1  = 1 - ROCVSPENS 

X11 = 1-roccurve$specificities 
Y11 = roccurve$sensitivities 
X11= sort(unique(X1))

par(mar = c(5, 5,3,2) +0.01)
fit <- smooth.spline(X1, Y, nknots = 10)
pred <- stats:::predict.smooth.spline(fit, X11)$y  

plot(X1[RUNV==1], Y[RUNV==1],lwd=1.7,type='l',col='blue',main=paste0('(TCGA)Chromothripsis from genomic sigs ', "AUC=",aa),
     ylab='True Positive',xlab='False Positive', cex.lab = 1.5, cex.axis = 1.3, mgp = c(3, 1,0))
for (kk in 1:10){
  points(X1[RUNV==kk], Y[RUNV==kk],type='l',col='blue',lwd=1.7)
}

lines(X11, pred, lwd = 2, col = 2) 
abline(a=0,b=1)
save(kk2, file = "TCGA_cnvsigmodel_res.rda")

