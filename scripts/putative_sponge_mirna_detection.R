#
#	putative_sponge_mirna_detection.R
#	This script finds those RNA transcripts
#	which are significantly enriched of miRNA
#	seeds.
#
#	Author: Javier Montalvo-Arredondo.
#	Contact: buitrejma@gmail.com
#	Departamento de Bioinformática
#	Universidad Autónoma Agraria Antonio Narro (UAAAN).
#	Calzada Antonio Narro, 1923, CP XXXXX, San Buenavista,
#	Saltillo, Coahuila, México.
#	+52 1 844 3 521 904
#

library("Biostrings")

transform_mature_mirna_to_seed_sequence <- function(mature_mirnas_fasta_filename) {
	# Requires a DNAStringSet object with the mature miRNA sequences
	# Returns the reverse complement of the miRNA seeds sequence

	mature <- readRNAStringSet(mature_mirnas_fasta_filename,
				   format="fasta")
	mature <- mature[grepl("^xla-",
			       names(mature))]
	offset_k6mer <- subseq(mature,
			       start=3,
			       width=6)

	k6mer <- subseq(mature,
			start=2,
			width=6)

	k7mer_A1 <- subseq(mature,
			   start=1,
			   width=7)

	k7mer_m8 <- subseq(mature,
			   start=2,
			   width=7)

	k8mer <- subseq(mature,
			start=1,
			width=8)

	seeds <- list(
		      offset_k6mer = reverseComplement(offset_k6mer),
		      k6mer = reverseComplement(k6mer),
		      k7mer_A1 = reverseComplement(k7mer_A1),
		      k7mer_m8 = reverseComplement(k7mer_m8),
		      k8mer = reverseComplement(k8mer)
		      )
	seeds
}

map_mirna_seeds <- function(transcripts, seeds) {
	# Requires a DNAStringSet object with the trasncript sequences.
	# Requires a data frame with the miRNA's ID and seed sequence.
	# Calculates the number of seeds matched in transcript sequence
	# for every trasncript and for every miRNA seed.
	# It returns a data frame in TIDY format.

# CHECKPOINT TESTED... [OK]
#	print(transcripts)
#	head(seeds)

	seeds_types <- names(seeds)
	tidy_data <- data.frame()
	transcript_names <- names(transcripts)
	transcript_width <- width(transcripts)
	for (seed_type in seeds_types) {
		table_seed_type <- rep(seed_type, length(transcript_names))
		cat("\n")
		print(paste("Processing ", seed_type, "...", sep=""))

		counter <- 1
		pb <- txtProgressBar(min=counter,
				     max=length(names(seeds[[seed_type]])),
				     style=3)

		for (seed_name in names(seeds[[seed_type]])) {
			table_seed_name <- rep(seed_name,
					       length(transcript_names))
			tmp_table <- data.frame(transcript_names,
					   transcript_width,
					   table_seed_type,
					   table_seed_name)
			tmp_table <- cbind(tmp_table,
			seeds_per_transcript = vcountPattern(pattern = chartr(
								as.character(
									seeds[[seed_type]][seed_name]
									),
								old="U",
								new="T"
								),
								subject = transcripts,
								algorithm = "boyer-moore"))

			if (length(tidy_data) == 0) {
				tidy_data <- tmp_table
			}
			else {
				tidy_data <- rbind(tidy_data,
						   tmp_table)
			}
			setTxtProgressBar(pb, counter)
			counter <- counter + 1
		}
	}
# CHECKPOINT TESTED... [OK]
	tidy_data$seeds_per_1kb <- tidy_data$seeds_per_transcript * (1000 / tidy_data$transcript_width)
	tidy_data$seeds_per_1kb <- round(tidy_data$seeds_per_1kb, 0)
	cat("\n")
	print("Finished Task!")
	tidy_data
}

calculate_significance <- function(mapped_seeds) {
	# Requires a data frame with the names of transcripts, miRNA
	# and the number of matched seeds in sequence by kb or
	# transcript.
	lambda_spt <- sum(mapped_seeds$seeds_per_transcript) / length(mapped_seeds$seeds_per_transcript)
	lambda_spkb <- sum(mapped_seeds$seeds_per_1kb) / length(mapped_seeds$seeds_per_1kb)

	# Seeds per transcript (spt) p-value
	mapped_seeds$p_value_spt <- ppois(q = mapped_seeds$seeds_per_transcript,
					  lambda = lambda_spt,
					  lower.tail = F,
					  log.p = F)

	# Seeds per transcript (spt) q-value
	mapped_seeds$q_value_spt <- p.adjust(p = mapped_seeds$p_value_spt,
					     method = "BH")

	# Seeds per transcript (spt) -log(p-value)
	mapped_seeds$log_p_value_spt <- log(x = mapped_seeds$p_value_spt,
					    base = 10) * -1

	# Seeds per transcript (spt) -log(q-value)
	mapped_seeds$log_q_value_spt <- log(x = mapped_seeds$q_value_spt,
					    base = 10) * -1

#########################################################################################

	# Seeds per 1kb (spkb) p-value
	mapped_seeds$p_value_spkb <- ppois(q = mapped_seeds$seeds_per_1kb,
					   lambda = lambda_spkb,
					   lower.tail = F,
					   log.p=F)

	# Seeds per 1kb (spkb) q-value
	mapped_seeds$q_value_spkb <- p.adjust(p = mapped_seeds$p_value_spkb,
					      method = "BH")

	# Seeds per 1kb (spkb) -log(p-value)
	mapped_seeds$log_p_value_spkb <- log(x = mapped_seeds$p_value_spkb,
					     base = 10) * -1

	# Seeds per 1kb (spkb) -log(q-value)
	mapped_seeds$log_q_value_spkb <- log(x = mapped_seeds$q_value_spkb,
					     base = 10) * -1

##########################################################################################
	# Lambda parameters for both metrics
	attr(mapped_seeds, "parameters") <- list(lambda_spt = lambda_spt,
						 lambda_spkb = lambda_spkb)
	# Returns data frame
	mapped_seeds
}

topN <- function(mapped_seeds_with_signifiqueishon, thresshold) {
	# Requires a dataframe with the number of matched seeds
	# of every transcript and the p-values and q-values.
	# Returns a data frame with all transcript, that has
	# q-values lower or equal a thresshold.
	tmp <- mapped_seeds_with_signifiqueishon[mapped_seeds_with_signifiqueishon$q_value_spkb <= thresshold
					  & mapped_seeds_with_signifiqueishon$transcript_width >= 1000,]
	tmp[order(tmp$log_q_value_spkb, decreasing=T),]
}

identify_putative_sponge_mirnas <- function(filenames, selected_seed_type = NULL) {
	# This is the wrapper function. Takes the file names
	# where transcripts and seed sequences are.
	# It makes the mapping and statistical process.
	# It returns the topN table.
	# filenames argument is a character vector with the
	# names of complete path to transcript and seed sequences
	# files. The first position is for transcripts and the second
	# position is for seed sequences fila.
	

	transcripts <- readDNAStringSet(filenames[1], format="fasta")
	seeds <- transform_mature_mirna_to_seed_sequence(filenames[2])

	if ( ! is.null(selected_seed_type)) {
		tmp <- list()
		tmp[[selected_seed_type]] <- seeds[[selected_seed_type]]
		seeds <- tmp
	}

#	print(head(seeds))

	mapped_seeds <- map_mirna_seeds(transcripts, seeds)
#CHECKPOINT TESTED... [OK]
	mapped_seeds <- calculate_significance(mapped_seeds)
	mapped_seeds
	
}
