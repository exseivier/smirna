#################################################################
#	extract_three_UTR.R
#	extracts UTR sequences from gff3 and fasta files
#	
#	Author: Javier Montalvo-Arredondo.
#	Version: 0.0
#	Contact: buitrejma@gmail.com
#	Departamento de Biotecnología.
#	Universidad Autónoma Agraria Antonio Narro.
#	UAAAN.
#	Calzada Antonio Narro 1923, CPXXXXX, San Buenavisat
#	Saltillo Coahuila, México.
#
#################################################################

library(Biostrings)
library(rtracklayer)
library(GenomicFeatures)

# Imports gff3 file.
extract_UTR_sequences <- function(gff_file,
				  BSgenome,
				  filter_tag = NULL,
				  UTR_length_thresshold) {

#	library(package = BSgenome_package_name, lib.loc = BSgenome_library_location)

	gff <- makeTxDbFromGFF(file = gff_file,
			       format = "gff3")
	
	threeUTR <- threeUTRsByTranscript(gff)
	keys <- threeUTR@unlistData@elementMetadata@listData$exon_name
	print(head(keys))
	
	gff_data_table <- select(gff,
				 keys = keys,
				 columns=c("EXONNAME",
					   "EXONID",
					   "GENEID",
					   "TXID",
					   "TXNAME"),
				 keytype="EXONNAME")

	print(head(gff_data_table))

	# Namming UTR by gene name
	threeUTR_names <- paste(gff_data_table[!duplicated(gff_data_table$TXID),"GENEID"],
				gff_data_table[!duplicated(gff_data_table$TXID), "TXID"],
				gff_data_table[!duplicated(gff_data_table$TXID), "TXNAME"],
				sep = "_")

	print(head(threeUTR_names))

	stopifnot(length(threeUTR_names) == length(names(threeUTR)))

	print("threeUTR names length was checked and it passed!")
	names(threeUTR) <- threeUTR_names

	if (!is.null(filter_tag)) {
		threeUTR <- threeUTR[grepl(filter_tag, seqnames(threeUTR))]
	}

	print("threeUTR data was filtered")
	
	threeUTR_seqs <- extractTranscriptSeqs(x = BSgenome,
					       transcripts = threeUTR)

	print("Three UTR sequences were extracted")
	threeUTR_seqs <- threeUTR_seqs[width(threeUTR_seqs) >= UTR_length_thresshold]

	threeUTRsData <- list(seqs = threeUTR_seqs, gff = threeUTR)
	threeUTRsData
}

export_GRangeList_DNAStringSet_to_files <- function(DNAStringSet_obj = NULL,
				       GRangesList_obj = NULL,
				       output_gff = "output.gff3",
				       output_fasta = "output.seqs.txt") {
	# CLAUSULAS
	condition1 <- !is.null(DNAStringSet_obj)
	condition2 <- !is.null(GRangesList_obj)
	condition3 <- condition1 & condition2
	stopifnot(condition3 == TRUE)

	print("Exporting GRangesList object to file...")
	GRangesList_obj <- unlist(GRangesList_obj)
	export(object = GRangesList_obj,
		    con = output_gff,
		    format = "gff3")

	print("Exporting DNAStringSet object to file...")
	writeXStringSet(x = DNAStringSet_obj,
			filepath = output_fasta,
			format = "fasta")



	print("Both objects were successfully exported.")
}

export_UTR_seqs_as_a_table <- function(DNAStringSet_obj, output = "out", organismID = "1001") {
	#
	# Exporting sequences as a table.
	#
	names <- names(DNAStringSet_obj)
	orgID <- rep(organismID, length(names))
	seqs <- as.character(DNAStringSet_obj)
	data_table <- data.frame(name = names, orgID = orgID, seq = seqs)
	write.table(data_table,
		    file = paste(output, "seqs.tab", sep="."),
		    quote = FALSE,
		    row.names = FALSE,
		    col.names = FALSE,
		    sep = "\t")
}

collapse_sequences <- function(DNAStringSet_obj, biggest_to_lowest = TRUE) {
	common_gene_names <- unique(
				    sort(
					 gsub("_[[:digit:]]*_[[:alpha:]]*.*",
					      "",
					      names(DNAStringSet_obj))
					 )
				    )
	
	collapsed_sequences <- DNAStringSet()
	pb <- txtProgressBar(min=0, max=length(common_gene_names), style=3)
	i <- 1
	for (gene_name in common_gene_names) {
		
		gene_name <- gsub("\\]",
				  '\\\\\\]',
				  gsub("\\[",
				       '\\\\\\[',
				       gene_name)
				  )
		gene_name <- paste("^", gene_name, sep="")
		tmp_DNAStringSet_obj <- DNAStringSet_obj[grepl(gene_name,
							       names(DNAStringSet_obj))]

		tmp_DNAStringSet_obj <- tmp_DNAStringSet_obj[
							     order(width(tmp_DNAStringSet_obj),
								   decreasing = biggest_to_lowest)
							     ]

		if (!length(tmp_DNAStringSet_obj) >= 1) {
			cat("\n")
			print(tmp_DNAStringSet_obj)
			print(gene_name)
			stop("Execution error!")
		}

		collapsed_sequences <- c(collapsed_sequences, tmp_DNAStringSet_obj[1])
		i <- i + 1
		setTxtProgressBar(pb, i)
	}
	collapsed_sequences
}
