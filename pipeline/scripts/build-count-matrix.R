suppressPackageStartupMessages({
	library("readr")
	library("data.table")
	library("dplyr")
	library("tximport")
	library("GenomicFeatures")
})

file.not.empty<-function(filenames){file.info(filenames)$size != 0}

# IPUTS

args = commandArgs(trailingOnly=TRUE)

inpdir<-args[1]
outdir<-args[2]
gtf<-args[3]
metadat<-args[4]
cancer<-args[5]

# MAIN 

print("Reading sample metadata...")
samp.info<-read.csv(metadat,header=T,row.names=1)


#Writing txi objects
################################################


print("=============Input directory============")
print(inpdir)

#Reading sample metadata

print("Reading sample files")
dirs<-list.dirs(inpdir,recursive=FALSE,full.names=F)
valid_ids<-lapply(dirs,function(d,...){
				file<-paste(inpdir,as.character(d),"quant.sf",sep="/")
				if(file.exists(file)&file.not.empty(file)){
					return(d)
				}
			})%>%unlist()
print(length(valid_ids))

stats<-data.frame(cancer=cancer,n_samples=nrow(samp.info),n_samples_salmon=length(valid_ids))
samp.info<-samp.info[valid_ids,]

# write metadata of valid files
metadat.out<-sub(".csv",".salmon.csv",metadat)
write.csv(samp.info,file=metadat.out)

#files<-file.path(inpdir,rownames(samp.info),"quant.sf")

#Build tximport object
#print("Reading files from input directory...")
#cmode<-"raw"
#txi<-tximport(files,type="salmon",txOut=TRUE,countsFromAbundance="no",dropInfReps=TRUE)
#colnames(txi$counts)<-rownames(samp.info)

#Save tximport object

#objpath<-paste(outdir,"/",cancer,".txiobject.RData",sep="")
#save(txi,file=objpath)
#print(paste("Writing object to:",objpath))

#sumpath<-paste(outdir,"/",cancer,".counts.summary.csv",sep="")
#write.csv(colSums(txi$counts),file=sumpath)

# Write summary stats
#statfile<-paste("../data/tcga_metadata/",cancer,".stats.csv",sep="")
#write.csv(stats,file=statfile,quote=FALSE,row.names=FALSE)



