


prepare_BL_PCT_file <- function(ts_out_file) {
	data_table <- read.delim(ts_out_file,
				 header = TRUE,
				 as.is = FALSE,
				 sep = "\t")
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
