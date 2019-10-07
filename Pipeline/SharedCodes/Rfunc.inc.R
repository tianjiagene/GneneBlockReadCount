
geno_len <- list(hg18=3107677273, hg19=3137161264, mm8=2664455088, mm9=2725765481, rn4=2834127293)


tax_ids        <- c(9606,9606,9606,9606,9606,                  10090,10090,10090,10090,10090
	,10116,10116,10116,10116,              9031,9031,9031,9031
	,9544,9544,9544
	,9615,9615,9615)
names(tax_ids) <- c("hs","human","Homo_sapiens","hg18","hg19", "mm","mouse","Mus_musculus","mm8","mm9"
	,"rn","rat","Rattus_norvegicus","rn4", "galGal","galGal3","chicken","Gallus_gallus"
	,"monkey","rheMac2","Macaca_mulatta"
	,"dog","canFam3","Canis_lupus_familiaris")
spe_common_names <- c("human","mouse","rat","chicken","monkey","dog")
names(spe_common_names) <- tax_ids[spe_common_names]

spe_latin_names <- c("Homo_sapiens","Mus_musculus","Rattus_norvegicus","Macaca_mulatta","Canis_lupus_familiaris")
names(spe_latin_names) <- tax_ids[spe_latin_names]

mkdir_if_not_exist<-function(out_root){
	out_root_dir=sub("[^/]*$","",out_root)
	system(paste("mkdir -p \'",out_root_dir,"\'", sep=""))
}

string2list<-function(s){
	l<-list()
	for(v1 in unlist(strsplit(s,";"))){ #eg: "REGU:UP DN;XXXX:XX XXX
		v2<-unlist(strsplit(v1, ":")) #PA:P0 P1 P2 P3
		l[[v2[1]]] <- unlist(strsplit(v2[2], " "))
	}
	return(l)
}


change_values<-function(inv, from_patt, to_txt){ #change values in a vector
	out<-as.character(unlist(inv))
	for(i in 1:length(from_patt)){
		out<-gsub(from_patt[i], to_txt[i], out)
	}
	return(out)
}



gid2_gene_desc <- function(gids, spe="human", gene_info_file=NA){
	geneinfo_data=open_geneInfo_data(spe, gene_info_file=gene_info_file)
	outtb <- geneinfo_data[match(gids, geneinfo_data$gene_id), c("gene_symbol","gene_desc")]
}
gsb2_gene_desc <- function(gsbs, spe="human", gene_info_file=NA){
	geneinfo_data=open_geneInfo_data(spe, gene_info_file=gene_info_file)
	outtb <- geneinfo_data[match(gsbs, geneinfo_data$alias), c("gene_id","gene_symbol","gene_desc")]
}
open_geneInfo_data<-function(spe="human", gene_info_file=NA){
	tax_id <- tax_ids[spe]
	spe_latin_name=spe_latin_names[as.character(tax_id)]
	if(is.na(gene_info_file)){ gene_info_file <- paste("01gene_alias2id/",spe_common_names[as.character(tax_id)],".",tax_id,".tbl",sep="") }
	if(! file.exists(gene_info_file)){
		gene_info_file <- paste("../ncbi/gene/geneinfo/All_Mammalia.gene_info",sep="")
		print(gene_info_file)
		geneinfo_data <- read.table(file=gene_info_file,header=F,quote="",comment.char="", sep="\t", skip=1,stringsAsFactors=F)
		names(geneinfo_data)[c(1,2,3,5,9)] <- c("tax_id","gene_id","gene_symbol","alias","gene_desc")
		geneinfo_data <- geneinfo_data[,c(1,2,3,5,9)]
		if(!is.na(tax_id)){
			geneinfo_data <- geneinfo_data[geneinfo_data$tax_id == tax_id,]
		}
		geneinfo_data[nrow(geneinfo_data)+1,] <- rep("-",ncol(geneinfo_data)) 
	}else{
		print(gene_info_file)
		geneinfo_data <- read.table(file=gene_info_file,header=T,quote="",comment.char="", sep="\t", stringsAsFactors=F)
	}
}

#id mapping
geneidmapping <- function(invec, intype="alias", inspe="Mus_musculus", outtype="gene_id", outspe="Homo_sapiens", outspetype="gid"){
	#type should be alias, gene_symbol or gene_id (header in the gsb2gid_file)
	#spe should be Homo_sapiens Mus_musculus ... (header in the homgene.64.tbl hid     gid_Mus_musculus        gsb_Mus_musculus        gid_Rattus_norvegicus   gsb_Rattus_norvegicus   gid_Magnaporthe_grisea...)
	gsb2gid_file <- paste(root_dir,"analyze/club/9gene_alias2id/",spe_common_names[as.character(tax_ids[inspe])],".",tax_ids[inspe],".tbl",sep="")
	if(intype==outtype){
		outvec1 <- invec
	}else{
		gsb2gid_data <- read.table(gsb2gid_file, header=T,sep="\t",quote="",comment.char="",stringsAsFactors=F)
		outvec1 <- gsb2gid_data[[outtype]][match(toupper(invec), toupper(gsb2gid_data[[intype]]) )]
		outvec1[is.na(invec)] <- NA
		print(table(is.na(invec)))
		print(table(is.na(outvec1)))
	}
	if(inspe==outspe){
		return(outvec1)
	}else{ #tranfer species
		homologene_file <- paste(root_dir,"analyze/club/8HomoloGene_tb/HomoloGene.build68.tbl",sep="")
		print(homologene_file)
		homologene_data <- read.table(homologene_file, header=T,sep="\t",quote="",comment.char="")
		s_inout_type <- sub("geneid|gene_id","gid",c(outtype,outspetype))
		s_inout_type <- sub("alias|genesymbol|gene_symbol","gsb",s_inout_type)
		print(s_inout_type)
		inheader <- paste(s_inout_type[1],inspe,sep="_")
		outheader <- paste(s_inout_type[2],outspe,sep="_")
		outvec2 <- homologene_data[[outheader]][match( toupper(outvec1), toupper(homologene_data[[inheader]]) )]
		print("table(is.na(outvec2))")
		print(table(is.na(outvec2)))
		return(outvec2)
	}
}

re_format_tb<-function(data, delete_headers=NULL, del_pattern=NULL, Inf_headers=NULL, front_headers=NULL, back_headers=NULL, ch_header_from=NULL,ch_header_to=NULL ){ #reformat table, delete columns, change Inf to 999, change headings, etc
	if(!is.null(delete_headers)){
		print (paste(c("delete:",delete_headers)) )
		data=data[,setdiff(names(data), delete_headers)]
	}
	if(!is.null(del_pattern)){
		print (paste(c("delete column pattern:",del_pattern)) )
		del_cols=grep(del_pattern,names(data))
		if(length(del_cols)>0){
			data=data[,-del_cols]
		}
	}
	if(!is.null(Inf_headers)){
		print (paste(c("Inf headers=:",Inf_headers)) )
		data[Inf_headers][data[Inf_headers]==Inf]= 999
		data[Inf_headers][data[Inf_headers]==-Inf]= -999
	}
	if(!is.null(front_headers)){
		front_headers=intersect(front_headers,names(data))
		data=data[,c(front_headers,setdiff(names(data), front_headers))]
	}
	if(!is.null(back_headers)){
		back_headers=intersect(back_headers,names(data))
		data=data[,c(setdiff(names(data), back_headers), back_headers)]
	}
	if(!is.null(ch_header_from)){
		names(data)= change_values(names(data), ch_header_from, ch_header_to)
	}
	data
}

define_parallel_fun<-function(nCores=1){
	if (nCores > 1) {
			library("parallel")
		if (!is.loaded("mc_fork", PACKAGE = "parallel")) {
			stop("Please load first parallel package or set parameter nCores to 1...")
		}else {
			myApply <<- function(X, FUN) {
				parallel::mclapply(X, FUN, mc.cores = nCores)
			}
		}
	}else {
		myApply <<- lapply
	}
}

