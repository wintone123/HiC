# load library
library(tidyr)
library(dplyr)

# load info
path <- "/mnt/c/HiC/test2"
imput_csv <- "chr14_1.csv"
output_csv <- "chr14_40k.csv"
bin_size <- 40000
chosen_area <- c(40000000,50000000)

# load file
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

# choose area
bin_start <- chosen_area[1] %/% 40000 + 1
bin_end <- chosen_area[2] %/% 40000 +1
upstream <- bin_end + 50
downstream <- bin_start - 50

# binning & filting
print("----------binning & filtering......----------")
bins <- data.frame(start1 = import$start1 %/% bin_size + 1,
				   start2 = import$start2 %/% bin_size + 1,
				   stringsAsFactors = FALSE)
bins <- filter(bins, start1 >= downstream & start1 <= upstream & start2 >= downstream & start2 <= upstream)

# change position
print("----------changing position......----------")
for (i in 1:nrow(bins)) {
    if (bins$start1[i] > bins$start2[i]) {
        temp <- bins$start1[i]
        bins$start1[i] <- bins$start2[i]
        bins$start2[i] <- temp
    }
}
bins <- arrange(bins, start1, start2)

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
				print(paste0(data,"---",n))
				temp <- data.frame(position = data, score = n)
				scores <- rbind(scores, temp)
				cond = FALSE
			}
		} else {
			s = i 
			print(paste0(data,"---",n))
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
                     hits = rep(NA, length(c(bin_start: bin_end))))
for (i in nrow(output)) {
    bin <- c(bin_start: bin_end)[i]
    A_temp <- filter(scores, position >= bin-50 & position <= bin-1)
    A <- sum(A_temp$score)
    B_temp <- filter(scores, position >= bin+1 & position <= bin+50)
    B <- sum(B_temp$score)
    print(paste0(A,"--",B))
    E <- (A + B) / 2
    output$hits[i] <- (B-A)/(ifelse(B>A, B-A, A-B))*(((A-E)^2)/E+((B-E)^2)/E)
}

# output data
print("----------writing csv......----------")
write.csv(output, file.path(path, output_csv))