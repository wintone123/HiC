# load library
genome <- library(BSgenome.Mmusculus.UCSC.mm10)

# load genome info
genome <- seqinfo(genome)
bin_size = 100000

# import data
nrow_import <- nrow(import)

# reorder
for (i in nrow(import_1)) {
	if (import_1$chr1[i] > import_1$chr2[i]) {
		temp1 <- import_1[i,1:3]
		temp2 <- import_1[i,4:6]
		import_1[i,1:3] <- temp1
		import_1[i,4:6] <- temp2
	} else if (import_1$chr1[i] == import_1$chr2[i] & import_1$start1[i] >  import_1$start2[i]) {
		temp1 <- import_1[i,1:3]
		temp2 <- import_1[i,4:6]
		import_1[i,1:3] <- temp1
		import_1[i,4:6] <- temp2
	}
}

# alignmnet
output_1 <- data.frame(chr1 = import$chr1,
					   start1 = import$start1 %/% bin_size * bin_size,
					   end1 = (import$end1 %/% bin_size + 1) * bin_size, 
					   chr2 = import$chr2,
					   start2 = import$start2 %/% bin_size * bin_size,
					   end2 = (import$end2 %/% bin_size + 1) * bin_size,
					   stringsAsFactors = FALSE)
chr_length_1 <- seqlength(genome)[output_1$chr1]
output_1[output_1$end1 >= chr_length_1, "end1"] <- chr_length_1[output_1$end1 >= chr_length_1]
chr_length_2 <- seqlength(genome)[output_1$chr2]
output_1[output_1$end2 >= chr_length_2, "end1"] <- chr_length_2[output_1$end2 >= chr_length_2]

# scoring
s = 1
n = 0
cond = TRUE
while (cond) {
	for (i in s:nrow(output_1)) {
		data <- output_1$temp[s]
		if (output_1$temp[i] == data) {
			n = n + 1
			if (i == nrow(output_1)) {
				print(paste0(data,n))
				cond = FALSE
			}
		} else {
			s = i 
			print(paste0(data,n))
			n = 0
 		    break
		}
	}
}

# output data