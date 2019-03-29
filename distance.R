# library
library(dplyr)
library(tidyr)

# load info
path <- "/mnt/c/HiC/test4"
imput_csv <- "chr1.csv"
output_name <- "chr1_100k"
bin_size <- 100000

# load file
cat("--------------loading csv--------------", "\n")
import <- read.delim(file.path(path, imput_csv), header = TRUE)

# binning
cat("----------------binning----------------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)

# filtering
cat("---------------filtering---------------", "\n")
bins <- filter(bins, start1 != start2)
bins_fil <- data.frame(chr1 = bins$chr1,
                       start1 = ifelse(bins$start1 <= bins$start2, bins$start1, bins$start2),
                       chr2 = bins$chr2,
                       start2 = ifelse(bins$start1 <= bins$start2, bins$start2, bins$start1),
                       stringsAsFactors = FALSE)
bins_fil <- arrange(bins_fil, chr1, start1, chr2, start2)
# cat(nrow(bins_fil), "\n")

# scoring
cat("----------------scoring----------------", "\n")
temp_df <- unite(bins_fil, temp, c(1:4), sep = "_")
scores <- data.frame(name = NA, score = NA)
s <- 1
n <- 0
j <- 0
k <- 0
cond = TRUE
while (cond) {
	for (i in s:nrow(temp_df)) {
		data <- temp_df$temp[s]
		if (temp_df$temp[i] == data) {
			n <- n + 1
			if (i == nrow(temp_df)) {
                # cat(paste0(data,"--",n))
				scores <- rbind(scores, data.frame(name = data, score = n))
				cond = FALSE
			}
		} else {
			s <- i 
            # cat(paste0(data,"--",n))
			scores <- rbind(scores, data.frame(name = data, score = n))
			n <- 0
 		    break
		}
        # progress bar
        if ((j / nrow(temp_df) * 100) %/% 1 == k) { 
            cat("[", paste(rep("#", (34*k/100) %/% 1),collapse = ""), 
            paste(rep("_", 34-(34*k/100) %/% 1),collapse = ""), k, "%]","\r", sep = "")
            k <- k + 1
        }
        j <- j + 1
    }
}
scores <- scores[2:nrow(scores),]
scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")
scores$start1 <- as.numeric(scores$start1)
scores$start2 <- as.numeric(scores$start2)
scores <- arrange(scores, score)

# scoring
cat("----------------scoring----------------", "\n")
for (i in 1:nrow(scores)) {
    dis <- scores$start2[i] - scores$start1[i]
    if (i == 1) {
        output <- data.frame(distance = dis, score = scores$score[i])
    } else {
        if (dis %in% output$distance) {
            output[output$distance == dis, "score"] <- output[output$distance == dis, "score"] + scores$score[i]
        } else {
            output <- rbind(output, data.frame(distance = dis, score = scores$score[i]))
        }
    }
}
output <- arrange(output, score)
output$percentage <- output$score / sum(output$score)

# output
cat("--------------writing csv--------------", "\n")
write.csv(output, paste0(path, "/", output_name, "_distance.csv"))