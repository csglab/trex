suppressPackageStartupMessages({
    library(getopt)
    library(optparse)
    library(dplyr)
    library(data.table)
    library(DESeq2)
})

# Input args

option_list <- list(
  make_option(c("-s", "--sample"), type='character',dest="sampleid",
              help="Sample identifier used to name the output files"),
  make_option(c("-e", "--event"), type='character',dest="evtype",
              help="Name of AS event type to quantify. Must be one of A3,A5,AF,AL,MX,RI,SE,AL"),
  make_option(c("-r","--refmap"), type="character",dest = "isomapfile",
              help="Path to the ioe file from SUPPA2"),
  make_option(c("-c","--counts"), type="character",dest = "countsfile",
              help="Path to the transcript counts tximport object file"), 
  make_option(c("-d", "--directory"), type='character',dest="dir",
              help="Directory with evenCounts")
)

options<-parse_args(OptionParser(option_list = option_list))
cancer <- options$sampleid
outdir <- options$dir
evtype<-options$evtype 
trex.map <- options$isomapfile
ev_cts<-options$countsfile

min_ecount<-10 # Min event counts
flags<-list("cancer"=cancer,"event_type"=evtype)

message("------------ Analyzing ",cancer,"------------")
message("---- Fitting models for ",evtype," events -----")

# Load data
map <- fread(trex.map,data.table = F,stringsAsFactors = F) %>% as.data.frame() %>% filter(event_type==evtype)
load(ev_cts,verbose = F)
bads<-colnames(trex.meta)[colSums(is.na(assays(trex.meta)[["counts"]]))!=0 | colSums(assays(trex.meta)[["counts"]])==0] %>%
      c(.,sub("A_","C_",.)) 
trex.meta<-trex.meta[,!colnames(trex.meta)%in%bads] # remove samples in wich either count type equals exactly zero in all samples or has NAs
ns<-ncol(trex.meta)/2

flags[["initial_num_samples"]]<-ns
flags[["initial_num_tumor"]]<-sum(trex.meta@colData$condition=="tumor")/2
flags[["initial_num_normal"]]<-sum(trex.meta@colData$condition=="normal")/2

# Remove outliers 
message("Found ",ns," samples")
tout<-which(trex.meta@colData$is_outlier.lof==TRUE & as.character(trex.meta@colData$condition)=="tumor")
flags[["num_tumor_outliers"]]<-length(tout)

if(length(tout)!=0){
    trex.meta<-trex.meta[,-1*tout]
    message("Removed ",ns-ncol(trex.meta)/2," tumor outliers")
}else{
    message("No tumor outliers present")
}

# Process covariates
trex.meta@colData$condition <- factor(as.character(trex.meta@colData$condition),levels=c("normal","tumor"))
trex.meta@colData$sex <- factor(ifelse(trex.meta@colData$gender=="male",0,1),levels = c(0,1))
trex.meta@colData$age <- trex.meta@colData$age_at_index
trex.meta@colData$age[is.na(trex.meta@colData$age)] <- mean(trex.meta@colData$age,na.rm=T)

trex.meta<-trex.meta[,!(is.na(trex.meta@colData$condition))]
message("Removed ",ns-ncol(trex.meta)/2," samples without metadata")
flags[["num_samples_no_metadata"]]<-ns-ncol(trex.meta)/2

# Remove samples without impurity estimations
ns<-sum(trex.meta@colData$condition=="tumor")/2
trex.meta<-trex.meta[,(!is.na(trex.meta@colData$impurity) & trex.meta@colData$condition=="tumor") | trex.meta@colData$condition=="normal"]
nt<-sum(trex.meta@colData$condition=="tumor")/2
nn<-sum(trex.meta@colData$condition=="normal")/2

message("Removed ",ns-nt," tumor samples without purity estimations")
flags[["num_samples_no_purity"]]<-ns-nt

# Filter events with low counts 
cts<-trex.meta@assays@data$counts
totcts<-cts[,grepl("A_",colnames(cts))]+cts[,grepl("C_",colnames(cts))]
trex.meta<-trex.meta[rowSums(totcts >= min_ecount) >= nn,]

message("Removed ",nrow(cts)-nrow(trex.meta)," events with low counts")
flags[["num_events_low_counts"]]<-nrow(cts)-nrow(trex.meta)
flags[["num_events_valid"]]<-nrow(trex.meta)
flags[["final_num_tumor"]]<-sum(trex.meta@colData$condition=="tumor")/2
flags[["final_num_normal"]]<-sum(trex.meta@colData$condition=="normal")/2

if(nt<3 | nn<3 ){
    
    flags[["models_fitted"]]<-NA
    print(flags)
    outfile<-paste0(outdir,"/",evtype,".flags.tsv")
    write.table(as.data.frame(flags),file=outfile,sep='\t',row.names = F)
    stop("Not enough samples, cannot fit models.")
}

flags[["models_fitted"]]<-c()

# Analysis with covariates
#######################################

message("Fitting model with covariates...")
analysis<-"condition"

nf_n<-sum(trex.meta@colData$gender == "female" & trex.meta@colData$condition=="normal",na.rm=T)
nm_n<-sum(trex.meta@colData$gender == "male" & trex.meta@colData$condition=="normal",na.rm=T)
nf_t<-sum(trex.meta@colData$gender == "female" & trex.meta@colData$condition=="tumor",na.rm=T)
nm_t<-sum(trex.meta@colData$gender == "male" & trex.meta@colData$condition=="tumor",na.rm=T)

if(nf_n>3 & nm_n>3 & nf_t>3 & nm_t>3){
    message("Creating model matrix including sex")
    m <- model.matrix(~ trex_count + impurity + condition + sex + age + (impurity + condition + age + sex):trex_count, trex.meta@colData)
    flags[["models_fitted"]]<-paste(flags[["models_fitted"]],analysis,sep=";")
}else{
    message("Not enough samples to include sex in the model, creating model matrix without sex information")
    m <- model.matrix(~ trex_count + impurity + condition + age + (impurity + condition + age):trex_count, trex.meta@colData)
    flags[["models_fitted"]]<-paste(flags[["models_fitted"]],paste0(analysis,".nosex"),sep=";")
}

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = trex.meta@assays@data$counts[,rownames(m)], 
                              colData = colData(trex.meta)[rownames(m),],
                              design = m)
dds <- estimateSizeFactors(dds,type='poscounts')
dds <- DESeq(dds, quiet = FALSE,fitType = 'local')

# Save results
message("Saving dds object")
outfile<-paste0(outdir,"/",evtype,".dds.",analysis,".RData")
save(dds,file=outfile)

message("Saving res files")
contrasts<-resultsNames(dds)[c(-1,-2)]
for(contrast in contrasts){
    
    res <- results(dds,name=contrast) %>%
           as.data.frame() %>%
           tibble::rownames_to_column("event_id") %>%
           left_join(.,map %>% select(-total_transcripts),by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t')
    
    res <- lfcShrink(dds,coef=contrast,format = "DataFrame",quiet=TRUE,type = "normal") %>%
                as.data.frame() %>%
                tibble::rownames_to_column("event_id") %>%
                left_join(.,map %>% select(-total_transcripts),by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.lfcShrink.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t',row.names = F)
}

message("Writing flags...")
outfile<-paste0(outdir,"/condition.",evtype,".flags.tsv")
write.table(as.data.frame(flags),file=outfile,sep='\t',row.names = F)
message("Completed analysis!")