


anabel <- function(seqs) {

	#	seqs is a DNA o RNAStringSet object
	#	with the sequences.
	anagram_seqs <- RNAStringSet(
				sapply(
				       seqs, function(x){
				       			x[sample(1:length(x),
								replace=FALSE)]
				       			}
				       )
	)
	anagram_seqs

}
