
library("Biostrings")
library("BSgenome")
library("GenomicFeatures")

extract_sequences <- function(BSgenome,
			      gff_file,
			      filter_tag = "chr",
			      UTR_length_thresshold,
			      output = "out",
			      organismID = "0001") {

	#
	#	Extracts the 3' UTR and CDS sequences based on
	#	annotations from a GFF3 file. Requires a environment
	#	variable which stores a BSgenome object with the reference
	#	genome. Also it requires the gff3 file name, and a UTR length
	#	thresshold.
	#

	gff <- makeTxDbFromGFF(file = gff_file,
			       format = "gff3")

	cat("GFF file loaded...\n")
	print(gff)

	cat("Extracting 3UTR and CDS coordinates from GFF file...\n")
	threeUTR <- threeUTRsByTranscript(gff)
	cds <- cdsBy(gff,
		     by = "tx")

	cat("Tracking gene names...\n")
	keys <- threeUTR@unlistData@elementMetadata@listData$exon_name
	cds_names <- cds@unlistData@elementMetadata@listData$cds_name

	
	gff_data_table <- select(gff,
				 keys = keys,
				 columns=c("EXONNAME",
					   "EXONID",
					   "GENEID",
					   "TXID",
					   "TXNAME"),
				 keytype="EXONNAME")


	# Namming UTR by gene name
	threeUTR_names <- paste(gff_data_table[!duplicated(gff_data_table$TXID),"GENEID"],
				gff_data_table[!duplicated(gff_data_table$TXID), "TXID"],
				gff_data_table[!duplicated(gff_data_table$TXID), "TXNAME"],
				sep = "_")

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



	stopifnot(length(threeUTR_names) == length(names(threeUTR)))
	stopifnot(length(cds_names) == length(names(cds)))
	names(threeUTR) <- threeUTR_names
	names(cds) <- cds_names


	if (!is.null(filter_tag)) {
		cat("filter_tag was not null.\nFiltering 3UTR and CDS annotations...\n")
		threeUTR <- threeUTR[grepl(filter_tag, seqnames(threeUTR))]
		cds <- cds[grepl(filter_tag, seqnames(cds))]
	}

	cat("Extracting sequences form BSgenome...\n")
	threeUTR_seqs <- extractTranscriptSeqs(x = BSgenome,
					       transcripts = threeUTR)
	cds_seqs <- extractTranscriptSeqs(BSgenome,
					  transcripts = cds)

	threeUTR_seqs <- threeUTR_seqs[width(threeUTR_seqs) >= UTR_length_thresshold]
	cds_seqs <- cds_seqs[width(cds_seqs) > 0]

	cat("Collapsing 3UTR sequences...\n")
	threeUTR_seqs_collapsed <- collapse_sequences(threeUTR_seqs)
	cat("Colapsing CDS sequences...\n")
	cds_seqs_collapsed <- collapse_sequences(cds_seqs)

	cat("Returning data...\n")
	# Data returned
	list(utrs_seqs = threeUTR_seqs,
	     utrs_seqs_collapsed = threeUTR_seqs_collapsed,
	     utrs_gff = threeUTR,
	     cds_seqs = cds_seqs,
	     cds_seqs_collapsed = cds_seqs_collapsed,
	     cds_gff = cds)

}

prepare_mirna_files_to_targetscan <- function(mature_mirnas_fasta_file,
					      output = "out",
					      organismID = "0001") {
	#
	#	Takes the fasta file of mature miRNA
	#	sequences and transforms it into
	#	the two miRNA family and mature miRNA
	#	files in TargetScan v7.0 format.
	#

	#	Reads fasta file
	seqs <- readRNAStringSet(mature_mirnas_fasta_file,
				 format = "fasta")

	#	Extracting Xenopus laevis miRNAs
	seqs <- seqs[grepl("^xla-", names(seqs))]

	#	Subsetting sequences
	subseqs <- subseq(x = seqs,
			  start = 2,
			  width = 7)

	#	Getting and stripping miRNA names
	mir_names <- names(seqs)
	mir_names <- gsub(" [[:alpha:]]*[[:digit:]]*.*",
			  "",
			  mir_names)

	#	Creating organism ID vector
	orgID <- rep(organismID,
		     length(mir_names))

	#	Getting sequences as a character vector.
	mature_seqs <- as.character(seqs)
	names(mature_seqs) <- NULL
	seed_seqs <- as.character(subseqs)
	names(seed_seqs) <- NULL

	#	Creating data frames
	mirna_table_to_ts_file <- data.frame(name = seed_seqs,
					     seq = seed_seqs,
					     orgID = orgID)

	mirna_table_to_ts_file <- mirna_table_to_ts_file[!duplicated(mirna_table_to_ts_file),]

	mirna_table_to_tsCTS_file <- data.frame(seed_name = seed_seqs,
					  orgID = orgID,
					  name = mir_names,
					  seq = mature_seqs)

	#	Writing data frames to files
	write.table(x = mirna_table_to_ts_file,
		    file = paste(output, "ints", sep="."),
		    row.names = FALSE,
		    col.names = FALSE,
		    quote = FALSE,
		    sep = "\t")

	write.table(x = mirna_table_to_tsCTS_file,
		    file = paste(output, "intsCTS", sep="."),
		    row.names = FALSE,
		    col.names = FALSE,
		    quote = FALSE,
		    sep = "\t")

}

collapse_sequences <- function(DNAStringSet_obj,
			       biggest_to_lowest = TRUE) {

	#
	#	Requires a DNAStringSet object with the sequences to collapse
	#	This function finds all isoforms of transcrips for every gene
	#	and takes only the biggest (nt) isoform if biggest_to_lowest = TRUE
	#	otherwise, it takes the lowest (nt) isoform of the transcripts.
	#

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
	cat("\n")
	collapsed_sequences
}

write_sequences_as_a_table <- function(seqs, output = "out", organismID = "0001") {
	
	#
	#	Requires a DNAStringSet object (seqs) with the sequences
	#	to be stored in text file as a table.
	#	This function takes DNAStringSet object with the sequences
	#	stripps it down and all the information is allocated in
	#	a dataframe with the name of sequence, the organism ID and
	#	the sequence all separated by "tabs"
	#

	names <- names(seqs)
	orgID <- rep(organismID, length(names))
	seqs <- as.character(seqs)
	data_table <- data.frame(name = names, orgID = orgID, seq = seqs)
	write.table(data_table,
		    file = paste(output, "seqs.tab", sep="."),
		    quote = FALSE,
		    row.names = FALSE,
		    col.names = FALSE,
		    sep = "\t")

}

calculate_ORF_length <- function(DNAStringSet_obj, organismID = "0001") {
	
	ORF_length <- data.frame(gene_name = names(DNAStringSet_obj),
				 orgID = rep(organismID, length(DNAStringSet_obj)),
				 length = width(DNAStringSet_obj))
	ORF_length

}

prepare_BL_PCT_file <- function(ts_out_file) {

	#
	#	This is the trickiest part of the pipeline.
	#	It transforms the predicted targets output file
	#	after targetscan_70.pl analysis into BL_PCT file,
	#	that is the input of targetscan_70_context_score.pl
	#

	data_table <- read.delim(ts_out_file,
				 header = TRUE,
				 as.is = FALSE,
				 sep = "\t")
	
	#
	#	I don not know how to handle this.
	#	I will see how to solve it from
	#	read.delim function.
	#	The problem is that at species.ID column,
	#	it stores numeric values that begins with
	#	three zeros "000", so it is transformed
	#	to numeric value and read.delim erases the
	#	trailing zeros.
	#
	#	So I transform the values of that column with
	#	the following instruction.

	data_table$species_ID <- rep("0001", length(data_table[,1]))


	print(head(data_table))

	data_table <- data_table[,1:11]
	zeros <- rep(0, length(data_table[,1]))
	data_table$BranchLS <- zeros
	data_table$Pct <- zeros
	data_table$Conserved <- rep("", length(data_table[,1]))
	write.table(x = data_table,
		    file = "targets_BL_PCT.tab",
		    row.names = FALSE,
		    col.names = TRUE,
		    quote = FALSE,
		    sep = "\t")
}

make_a_list_from_dataframe <- function(data, factor) {
	
	#
	#	This function takes a dataframe object and
	#	transforms it into list	based on elements of a factor.
	#	Each element of the list corresponds to a dataframe
	#	with data related to one element of the factor vector.
	#

	mir_names <- unique(sort(data[, factor]))
	data_list <- list()
	for (name in mir_names) {
		data_list[[name]] <- data[data[,factor] == name,]
	}
	data_list
}
