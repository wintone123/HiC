# load library
library(tidyr)
library(dplyr)

# load info
path <- "/mnt/c/HiC/test2"
imput_csv <- "chr14_1.csv"
output_csv <- "chr14_40k_DI.csv"
output_scores <- "chr14_40k_score.csv"
bin_size <- 20000
extra_size <- 100000
chosen_area <- c(65000000,75000000)

# load file
import <- read.delim(file.path(path, imput_csv), header = TRUE)

# choose area
bin_start <- chosen_area[1] %/% bin_size + 1
bin_end <- chosen_area[2] %/% bin_size +1
extra_bin <- extra_size %/% bin_size
upstream <- bin_end + extra_bin
downstream <- bin_start - extra_bin

# binning 
print("----------binning......----------")
bins <- data.frame(start1 = import$start1 %/% bin_size + 1,
			       start2 = import$start2 %/% bin_size + 1)

# filting
print("---------filtering......---------")
bins <- filter(bins, start1 != start2)
bins_fil <- filter(bins, start1 >= downstream & start1 <= upstream & start2 >= downstream & start2 <= upstream)
bins_fil <- rbind(bins_fil, filter(bins, start1 >= downstream & start1 <= upstream & start2 < downstream & start2 > upstream))
bins_fil <- rbind(bins_fil, filter(bins, start1 < downstream & start1 > upstream & start2 >= downstream & start2 <= upstream))
bins_fil <- arrange(bins_fil, start1, start2)

# scoring
print("----------scoring......----------")
scores <- data.frame(position = c(downstream:upstream), score = 0)
for (i in 1:nrow(bins_fil)) {
	scores[scores$position == bins$start1[i],"score"] = scores[scores$position == bins$start1[i],"score"] + 1
	scores[scores$position == bins$start2[i],"score"] = scores[scores$position == bins$start2[i],"score"] + 1

}

# output data
print("--------writing csv......--------")
write.csv(scores, file.path(path, output_scores))

# calculating
print("--------calculating......--------")
output <- data.frame(position = c(bin_start: bin_end),
                     score = rep(0, length(c(bin_start: bin_end))))
for (i in 1:nrow(output)) {
    bin <- c(bin_start: bin_end)[i]
    A <- sum(scores[scores$position %in% c((bin-extra_bin):(bin)),2])
    B <- sum(scores[scores$position %in% c((bin):(bin+extra_bin)),2])
    # E <- (A + B) / 2
    # output$hits[i] <- ifelse(B>A, 1, -1)*(((A-E)^2+(B-E)^2)/E)
	output$score[i] <- ifelse(B>A, 1, -1)*(A-B)^2/(A+B)
}

# output data
print("--------writing csv......--------")
write.csv(output, file.path(path, output_csv))