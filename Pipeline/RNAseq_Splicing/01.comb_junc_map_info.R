##combine all junction read mapping info for multi samples 
source(paste("../SharedCodes/Rfunc.inc.R",sep=""))

args<- commandArgs(TRUE)
print(args)

args_v<-NA
if(length(args)>1){
	args_v <- args[c(FALSE, TRUE)]
	names(args_v) <- sub("\\-+","", args[c(TRUE, FALSE)])
	print(args_v)
}


study_name<- ifelse(is.na(args_v["study_name"]), "example_project", args_v["study_name"])
geno<- ifelse(is.na(args_v["geno"]), "hg19", args_v["geno"])
samples<- unlist(strsplit(as.character(args_v["samples"])," "))
id_header=ifelse(is.na(args_v["id_header"]), "refseqid", args_v["id_header"])

indir<- ifelse(is.na(args_v["indir"]), paste("04.GeneReadNum/",id_header,".",sep=""), args_v["indir"])
infext<- ifelse(is.na(args_v["infext"]), ".ReadNum.tbl.temp.junc2gene.tbl", args_v["infext"])
out_cbf<- ifelse(is.na(args_v["out_cbf"]), paste(indir,"combine.junc2gene.tbl",sep=""), args_v["out_cbf"])

if(id_header=="refseqid"){
	gene_ano_f=paste("../ReferenceDB/gene/02transcript_gene_ano/",geno,".refflat.desc.nonRedundant.txt",sep=""); 
	gene_info_headers=c("gene_symbol","gene_id","contig","strand","gene_desc")
}else if(id_header=="ensid"){
	gene_ano_f=paste("../ReferenceDB/gene/02transcript_gene_ano/",geno,".ensGene.desc.nonRedundant.txt",sep=""); 
	gene_info_headers=c("gene_symbol","gene_Biotype","contig","strand","gene_desc")
}
if(!is.na(args_v["gene_ano_f"])){ gene_ano_f=args_v["gene_ano_f"]}

juncIdHeaders=c("readType",id_header,"row","juncpos5","juncpos3")
juncIdHeaders2=c("contig","strand","juncpos5","juncpos3")


mkdir_if_not_exist(out_cbf)
#1 load all data, generate all junction IDs
comb_d=NULL
for(sample in samples){
	inf=paste(indir,sample,infext,sep=""); print(inf)
	d=read.table(inf, header=T, sep="\t", quote="", stringsAsFactors=F)
	if(sample==samples[1]){
		comb_d=d
	}else{
		comb_d=merge(comb_d,d, all.x=T, all.y=T, by.y=juncIdHeaders)
	}
}


readnum_names=grep("^num_",names(comb_d))
comb_d[readnum_names][is.na(comb_d[readnum_names])]=0
colSums(comb_d[readnum_names])

##2, add gene annotation
gene_ano_d=read.table(gene_ano_f, header=T, sep="\t", stringsAsFactors=F, comment.char="", quote="")
comb_d=comb_d[(comb_d$row-1) %in% gene_ano_d$ori_rowID, ] #only study junctions mapped to non-redundant genes
comb_d[,gene_info_headers]=gene_ano_d[match(comb_d$row-1,gene_ano_d$ori_rowID), gene_info_headers]
ifAllIdMatched=all(comb_d[,id_header]==gene_ano_d[match(comb_d$row-1,gene_ano_d$ori_rowID), id_header])
if(!ifAllIdMatched){print("Error! ID not match with gene_ano_f. !!!")}


ifsel=!is.na(comb_d$gene_symbol) & comb_d$gene_symbol !=""; table(ifsel)
comb_d=comb_d[ifsel, ]

##3, output
write.table(comb_d, file=out_cbf, sep="\t", col.names=T, row.names=F, quote=F)
##re-open
#comb_d=read.table(out_cbf, header=T, sep="\t", quote="", comment.char="", stringsAsFactors=F)

