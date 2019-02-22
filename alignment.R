# load genome info
bin_size = 100000

# import data

# alignmnet
nrow_import <- nrow(import)
output_1 <- data.frame(chr1 = rep(NA, nrow_import),
					 start1 = rep(NA, nrow_import),
					 end1 = rep(NA, nrow_import),
					 chr2 = rep(NA, nrow_import),
					 start2 = rep(NA, nrow_import),
					 end2 = rep(NA, nrow_import),
					 stringsAsFactors = FALSE)
for (i in nrow_import) {
	# alignment part 1
	bin1_start <- import$start1[i] %/% bin_size * bin_size + 1
	bin1_end <- (import$end1[i] %/% bin_size + 1) * bin_size
	if (bin1_end >= length(genome[import$chr1[i]])) {
		bin1_end <- length(genome[import$chr1[i]])
	}

	# alignment part 2
	bin2_start <- import$start2[i] %/% bin_size * bin_size + 1
	bin2_end <- (import$end2[i] %/% bin_size + 1) * bin_size
	if (bin2_end >= length(genome[import$chr2[i]])) {
		bin2_end <- length(genome[import$chr2[i]])
	}

	# output
	output_1$chr1[i] <- chr1
	output_1$start1[i] <- bin1_start
	output_1$end1[i] <- bin1_end
	output_1$chr2[i] <- chr2
	output_1$start2[i] <- bin2_start
	output_1$end2[i] <- bin2_end
}

# data arrangement
output_1 <- unite(output_1, temp, c(1:6), sep = "_")
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