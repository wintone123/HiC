# load library
library(tidyverse)

# load genome info
genome_length <- 107043718 # hg18 chr14
bin_size = 100000 

# import data
path = "/mnt/c/HiC/test1"
import <- read.delim(file.path(path, "chr14.csv"), sep = ",", row.names = 1, header = TRUE)

# binning
print("----------binning......----------")
output_1 <- data.frame(chr1 = import$chr1,
					   start1 = import$position1 %/% bin_size + 1,
					   chr2 = import$chr2,
					   start2 = import$position2 %/% bin_size + 1,
					   stringsAsFactors = FALSE)
output_1 <- arrange(output_1, start1, start2)
output_1 <- filter(output_1, start1 != start2)

# scoring
print("----------scoring......----------")
temp_df <- unite(output_1, temp, c(1:4), sep = "_")
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
				print(paste0(data,"---",n))
				temp <- data.frame(name = data, score = n)
				output_2 <- rbind(output_2, temp)
				cond = FALSE
			}
		} else {
			s = i 
			print(paste0(data,"---",n))
			temp <- data.frame(name = data, score = n)
			output_2 <- rbind(output_2, temp)
			n = 0
 		    break
		}
	}
}
output_2 <- output_2[2:nrow(output_2),]
output_2 <- separate(output_2, name, c("chr1","start1","chr2","start2"), sep = "_")
output_2$start1 <- as.numeric(output_2$start1)
output_2$start2 <- as.numeric(output_2$start2)

# output data
print("----------writing csv......----------")
write.csv(output_2, file.path(path, "output_100k.csv"))
q()

# matrixing 
print("----------matrixing......----------")
chr_xy <- sort(unique(c(output_2$chr1, output_2$chr2)))
xy_list <- vector()
for (i in chr_xy) {
	length <- genome_length
	bin_num <- length %/% bin_size + 1
	xy_list <- append(xy_list, paste(i, 61:bin_num, sep = "_"))
}
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list
output_2 <- unite(output_2, chr1, c(1:2), sep = "_")
output_2 <- unite(output_2, chr2, c(2:3), sep = "_")
for (i in 1:nrow(output_2)) {
	chr1 <- output_2$chr1[i]
	chr2 <- output_2$chr2[i]
	mat[chr1,chr2] <- output_2$score[i] + mat[chr1,chr2]
	mat[chr2,chr1] <- output_2$score[i] + mat[chr2,chr1]

}

# output data
print("----------writing csv......----------")
write.csv(mat, file.path(path, "output2.csv"))