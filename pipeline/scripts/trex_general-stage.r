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

options <- parse_args(OptionParser(option_list = option_list))
cancer <- options$sampleid
outdir <- options$dir
evtype <- options$evtype 
trex.map <- options$isomapfile
ev_cts<-options$countsfile


min_ecount<-10 # Min event counts
flags<-list("cancer"=cancer,"event_type"=evtype)

message("------------ Analyzing ",cancer,"------------")
message("---- Fitting models for ",evtype," events -----")

# Load data
message("Loading event counts...")
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
message("Checking for outliers...")
tout<-which(trex.meta@colData$is_outlier.lof==TRUE & as.character(trex.meta@colData$condition)=="tumor")
flags[["num_tumor_outliers"]]<-length(tout)

if(length(tout)!=0){
    trex.meta<-trex.meta[,-1*tout]
    message("Removed ",ns-ncol(trex.meta)," tumor outliers")
}else{
    message("No tumor outliers present")
}

# Process covariates
trex.meta@colData$sex <- factor(ifelse(trex.meta@colData$gender=="male",0,1),levels = c(0,1))
trex.meta@colData$age <- trex.meta@colData$age_at_index

# Remove samples without impurity estimations
ns<-sum(trex.meta@colData$condition=="tumor")
message("Found ",ns," tumor samples")
trex.meta<-trex.meta[,(!is.na(trex.meta@colData$impurity) & trex.meta@colData$condition=="tumor") | trex.meta@colData$condition=="normal"]
nt<-sum(trex.meta@colData$condition=="tumor")

message("Removed ",ns-nt," tumor samples without purity estimations")
flags[["num_samples_no_purity"]]<-ns-nt

# Filter events with low counts 
cts<-trex.meta@assays@data$counts
totcts<-cts[,grepl("A_",colnames(cts))]+cts[,grepl("C_",colnames(cts))]
nn<-sort(table(nn<-sum(trex.meta@colData$stage)))[1]
message("Number of samples in the smallest group:", nn)
trex.meta<-trex.meta[rowSums(totcts >= min_ecount) >= nn,]

message("Removed ",nrow(cts)-nrow(trex.meta)," events with low counts")
flags[["num_events_low_counts"]]<-nrow(cts)-nrow(trex.meta)
flags[["num_events_valid"]]<-nrow(trex.meta)
flags[["final_num_tumor"]]<-sum(trex.meta@colData$condition=="tumor")/2
flags[["final_num_normal"]]<-sum(trex.meta@colData$condition=="normal")/2

if(nt<3){
    stop("Not enough samples with estimated purity, cannot fit models.")
}

# Remove normal samples
message("Discarding normal samples...")
trex.meta.org<-trex.meta # Save original object for later analysis
ns<-ncol(trex.meta.org)
trex.meta<-trex.meta[,trex.meta@colData$condition=="tumor"]
message("Removed ",ns-ncol(trex.meta)," normal samples")

# Flag number of samples per group 

message("Defining flags for sample sizes...")
include_sex<-TRUE
include_normal<-TRUE

n_stages<-table(trex.meta@colData$stage)
n_sex<-table(trex.meta@colData$gender)
n_normal<-sum(trex.meta.org@colData$condition=="normal")

message("Summary of group sample sizes")
message(" - Condition: ",paste(" ",names(table(trex.meta.org@colData$condition)),"=",table(trex.meta.org@colData$condition),sep=""))
message(" - Continous stage: ",paste(" S",names(n_stages),"=",n_stages,sep=""))
message(" - Sex: ",paste(" ",names(n_sex),"=",n_sex,sep=""))

# Check sample sizes

valid_stages<-n_stages>=1 # At least one sample per stage
if(sum(valid_stages)>=2){
    message(" ....Passed required number of samples per stage")
}else{
    stop("All tumor samples are the same stage, cannot fit this series of models.")
    message("Exiting")
}

valid_sex<-n_sex>=1 # At least three samples per sex group
if(sum(valid_sex)==2){
    message("...Passed required number of samples per sex group")
}else{
    include_sex<-FALSE
    message("...Failed required number of samples per sex group")
}

if(n_normal>=3){
    message("...Passed required number of samples per cond group")
}else{
    include_normal<-FALSE
    message("...Failed required number of samples per cond group")
}

message("Summary of flags:")
message("Include sex = ",include_sex)
message("Include normal = ",include_normal)

# Analysis with covariates - Continuous
#######################################

message("Fitting continuous stage model with covariates...")
analysis<-"tumorStageContinuous"

if(include_sex){
    message("Creating model matrix with sex")
    m <- model.matrix(~ trex_count + impurity + stage + sex + age + (impurity + stage + age + sex):trex_count, trex.meta@colData)
    flags[["models_fitted"]]<-paste(flags[["models_fitted"]],analysis,sep=";")

}else{
    message("Not enough samples to include sex in the model, creating model matrix without sex information")
    m <- model.matrix(~ trex_count + impurity + stage + age + (impurity + stage + age):trex_count, trex.meta@colData)
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
           left_join(.,map,by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t')
    
    message(":: Shrinking ",contrast)
    res <- lfcShrink(dds,coef=contrast,format = "DataFrame",quiet=TRUE,type = "normal") %>%
                as.data.frame() %>%
                tibble::rownames_to_column("event_id") %>%
                left_join(.,map %>% select(-total_transcripts),by = "event_id")
    outfile<-paste(outdir,"/",evtype,".res.lfcShrink.",analysis,".",contrast,".tsv",sep="")
    write.table(as.data.frame(res),file=outfile,quote=F,sep='\t',row.names = F)
}

message("Writing flags")
outfile<-paste0(outdir,"/stage.",evtype,".flags.tsv")
write.table(as.data.frame(flags),file=outfile,sep='\t',row.names = F)
message("Finished successfully!")
