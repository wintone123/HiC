# load library
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)

# load genome info
genome <- seqinfo(BSgenome.Mmusculus.UCSC.mm10)
bin_size = 100000

# import data

# reorder
# for (i in nrow(import)) {
# 	if (import$chr1[i] > import$chr2[i]) {
# 		temp1 <- import[i,1:3]
# 		temp2 <- import[i,4:6]
# 		import[i,1:3] <- temp2
# 		import[i,4:6] <- temp1
# 	} else if (import$chr1[i] == import$chr2[i] & import$start1[i] >  import$start2[i]) {
# 		temp1 <- import[i,1:3]
# 		temp2 <- import[i,4:6]
# 		import[i,1:3] <- temp2
# 		import[i,4:6] <- temp1
# 	}
# }

# alignmnet to bin
output_1 <- data.frame(chr1 = import$chr1,
					   start1 = import$start1 %/% bin_size + 1,
					   end1 = import$end1 %/% bin_size + 1, 
					   chr2 = import$chr2,
					   start2 = import$start2 %/% bin_size + 1,
					   end2 = import$end2 %/% bin_size + 1,
					   stringsAsFactors = FALSE)
output_1 <- arrange(output_1, chr1, start1, end1, chr2, start2, end2)

# scoring
temp_df <- unite(output_1, temp, c(1:6), sep = "_")
output_2 <- data.frame(name = NA, score = NA)
s = 1
n = 0
cond = TRUE
while (cond) {
	for (i in s:nrow(temp_df)) {
		data <- temp_df$temp[s]
		if (temp_df$temp[i] == data) {
			n = n + 1
			if (i == nrow(temp_df)) {
				temp <- data.frame(name = data, score = n)
				output_2 <- rbind(output_2, temp)
				cond = FALSE
			}
		} else {
			s = i 
			temp <- data.frame(name = data, score = n)
			output_2 <- rbind(output_2, temp)
			n = 0
 		    break
		}
	}
}
output_2 <- output_2[2:nrow(output_2),]
output_2 <- separate(output_2, temp, c("chr1","start1","end1","chr2","start2","end2"), sep = "_")
output_2$start1 <- as.numeric(output_2$start1)
output_2$end1 <- as.numeric(output_2$end1)
output_2$start2 <- as.numeric(output_2$start2)
output_2$end2 <- as.numeric(output_2$end2)

# matrixing data
chr_xy <- sort(unique(c(output_2$chr1, output_2$chr2)))
xy_list <- vector()
for (i in chr_xy) {
	length <- seqlength(genome)[i]
	bin_num <- length %/% bin_size + 1
	xy_list <- append(xy_list, paste(i, 1:bin_num, sep = "_"))
}
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- y_list
rownames(mat) <- x_list
for (i in nrow(output_2)) {
	chr1 <- output_2$chr1[i]
	start1 <- paste(chr1, output_2$start1[i], sep = "_")
	end1 <- paste(chr, output_2$end1[i], sep = "_")
	chr2 <- output_2$chr2[i]
	start2 <- paste(chr2, output_2$start2[i], sep = "_")
	end2 <- paste(chr2, output_2$end2[i], sep = "_")
	if (start1 == end1) {
		if (start2 == end2) {
			mat[start1,start2] = output_2$score[i] + mat[start1,start2]
			mat[start2,start1] = output_2$score[i] + mat[start2,start1]
		} else {
			mat[start1,start2] = output_2$score[i] + mat[start1,start2]
			mat[start1,end2] = output_2$score[i] + mat[start1,end2]
			mat[start2,start1] = output_2$score[i] + mat[start2,start1]
			mat[end2,start1] = output_2$score[i] + mat[end2,start1]
		}
	} else {
		if (start2 == end2) {
			mat[start1,start2] = output_2$score[i] + mat[start1,start2]
			mat[end1,start2] = output_2$score[i] + mat[end1,start2]
			mat[start2,start1] = output_2$score[i] + mat[start2,start1]
			mat[start2,end1] = output_2$score[i] + mat[start2,end1]
		} else {
			mat[start1,start2] = output_2$score[i] + mat[start1,start2]
			mat[start1,end2] = output_2$score[i] + mat[start1,end2]
			mat[end1,start2] = output_2$score[i] + mat[end1,start2]
			mat[end1,end2] = output_2$score[i] + mat[end1,end2]
			mat[start2,start1] = output_2$score[i] + mat[start2,start1]
			mat[end2,start1] = output_2$score[i] + mat[end2,start1]
			mat[start2,end1] = output_2$score[i] + mat[start2,end1]
			mat[end2,end1] = output_2$score[i] + mat[end2,end1]
		}
	}
}

# output data