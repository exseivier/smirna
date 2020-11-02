#
#
#	Extracts the ORFs sequences and length
#	this script is part of pre-processing task
#	to run TargetScan_v7.0.
#
#	Author: Javier Montalvo-Arredondo.
#	Contact: buitrejma@gmail.com
#	Calzada Antonio Narro, 1923, San Buenavista,
#	Saltillo Coahuila, México.
#	Biotechnology Division, Bioinformatics Area.
#
library("BSgenome")
library("Biostrings")
library("GenomicFeatures")

extract_ORF_set <- function(BSgenome, gff_file, filter_tag = NULL, organism = NULL) {

	gff <- makeTxDbFromGFF(file = gff_file,
			       format = "gff3",
			       organism = organism)
	
	cds <- cdsBy(gff,
		     by = "tx")
	
	cds_names <- cds@unlistData@elementMetadata@listData$cds_name
	cds_data_table <- select(gff,
				 keys = cds_names,
				 columns = c("CDSID",
					     "CDSNAME",
					     "GENEID",
					     "TXNAME",
					     "TXID"),
				 keytype = "CDSNAME")
	cds_names <- cds_data_table[,c("GENEID", "TXID", "TXNAME")]
	cds_names <- cds_names[!duplicated(cds_names),]
	cds_names <- paste(cds_names$GENEID,
			   cds_names$TXID,
			   cds_names$TXNAME,
			   sep="_")

	names(cds) <- cds_names

	if (!is.null(filter_tag)) {
		cds <- cds[grepl(filter_tag, seqnames(cds))]
	}

	cds_seqs <- extractTranscriptSeqs(BSgenome,
					  transcripts = cds)

	cds_seqs <- cds_seqs[width(cds_seqs) > 0]

	cds_seqs	
}

calculate_ORF_length <- function(DNAStringSet_obj, organismID = "1001") {
	
	ORF_length <- data.frame(gene_name = names(DNAStringSet_obj),
				 orgID = rep(organismID, length(DNAStringSet_obj)),
				 length = width(DNAStringSet_obj))
	ORF_length

}


write_ORF_set_information <- function(DNAStringSet_obj, ORF_length, organismID = "1001", output = "out") {

#	writeXStringSet(DNAStringSet_obj,
#			filepath = paste(output, "seq", sep="."),
#			format = "fasta")

	#
	# Wirte sequences as a table.
	#
	names <- names(DNAStringSet_obj)
	orgID <- rep(organismID, length(names))
	seqs <- as.character(DNAStringSet_obj)
	names(seqs) <- NULL
	data_table <- data.frame(name = names,
				 orgID = orgID,
				 seq = seqs)

	write.table(data_table,
		    file = paste(output, "seqs.txt", sep="."),
		    quote = FALSE,
		    row.names = FALSE,
		    col.names = FALSE,
		    sep = "\t")


	write.table(ORF_length,
		    file = paste(output, "len", sep="."),
		    quote = FALSE,
		    row.names = FALSE,
		    sep = "\t")
	
	print("Objects were successfully stored in the files")

}
