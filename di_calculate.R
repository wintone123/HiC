# load library
library(tidyr)
library(dplyr)

# load info
path <- "/mnt/c/HiC/test2"
imput_csv <- "chr14_1.csv"
output_csv <- "chr14_40k.csv"
bin_size <- 100000
extra_size <- 400000
chosen_area <- c(20000000,30000000)

# load file
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

# choose area
bin_start <- chosen_area[1] %/% 40000 + 1
bin_end <- chosen_area[2] %/% 40000 +1
extra_bin <- extra_size %/% bin_size
upstream <- bin_end + extra_bin
downstream <- bin_start - extra_bin

# binning 
print("----------binning......----------")
bins_temp <- data.frame(start1 = import$start1 %/% bin_size + 1,
				   start2 = import$start2 %/% bin_size + 1,
				   stringsAsFactors = FALSE)

# change position & filting
print("----------exchanging & filtering......----------")
bins <- data.frame(start1 = ifelse(bins_temp$start1 < bins_temp$start2, bins_temp$start1, bins_temp$start2),
				   start2 = ifelse(bins_temp$start1 < bins_temp$start2, bins_temp$start2, bins_temp$start1))
bins <- arrange(bins, start1, start2)
bins <- filter(bins, start1 != start2)
bins <- filter(bins, start1 >= downstream & start1 <= upstream)

# scoring
print("----------scoring......----------")
scores <- data.frame(position = NA, score = NA)
s = 1
n = 0
cond = TRUE
while (cond) {
	for (i in s:nrow(bins)) {
		data <- bins$start1[s]
		if (bins$start1[i] == data) {
			n = n + 1
			if (i == nrow(bins)) {
				temp <- data.frame(position = data, score = n)
				scores <- rbind(scores, temp)
				cond = FALSE
			}
		} else {
			s = i 
			temp <- data.frame(position = data, score = n)
			scores <- rbind(scores, temp)
			n = 0
 		    break
		}
	}
}
scores <- scores[2:nrow(scores),]

# calculating
print("----------calculating......----------")
output <- data.frame(position = c(bin_start: bin_end),
                     hits = rep(0, length(c(bin_start: bin_end))))
for (i in 1:nrow(output)) {
    bin <- c(bin_start: bin_end)[i]
    A <- sum(scores[scores$position %in% c((bin-extra_bin):(bin-1)),2])
    B <- sum(scores[scores$position %in% c((bin+1):(bin+extra_bin)),2])
    # E <- (A + B) / 2
    # output$hits[i] <- ifelse(B>A, 1, -1)*(((A-E)^2+(B-E)^2)/E)
	output$hits[i] <- ifelse(B>A, 1, -1)*(A-B)^2/(A+B)
}

# output data
print("----------writing csv......----------")
write.csv(output, file.path(path, output_csv))