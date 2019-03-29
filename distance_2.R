# library
library(dplyr)
library(tidyr)

# function
add_scores <- function(temp_df) {
    scores <- data.frame(distance = NA, score = NA)
    s <- 1
    n <- 0
    j <- 0
    k <- 0
    cond = TRUE
    while (cond) {
        for (i in s:nrow(temp_df)) {
            data <- temp_df$distance[s]
            if (temp_df$distance[i] == data) {
                n <- n + 1
                if (i == nrow(temp_df)) {
                    # cat(paste0(data,"--",n))
                    scores <- rbind(scores, data.frame(distance = data, score = n))
                    cond = FALSE
                }
            } else {
                s <- i 
                # cat(paste0(data,"--",n))
                scores <- rbind(scores, data.frame(distance = data, score = n))
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
    scores <- arrange(scores, distance)
    return(scores)
}

# load info
path <- "/mnt/c/HiC/test4"
imput_csv <- "test_fil.csv"
output_name <- "chr1-X_100k"
min_size <- 100000
chosen_chr <- c(1:19,"X")

# load file
cat("--------------loading csv--------------", "\n")
import <- readr::read_delim(file.path(path, imput_csv), delim = "\t")

# filtering
cat("---------------filtering---------------", "\n")
import_fil <- data.frame(chr1 = as.character(import$chr1), 
                         chr2 = as.character(import$chr2),
                         distance = abs(import$start2 - import$start1) %/% min_size,
                         stringsAsFactors = FALSE)
import_fil <- filter(import_fil, chr1 == chr2 & distance >= 1)
# cat(nrow(import_fil), "\n")

# loop
cat("--------------loop start---------------", "\n")
scores <-data.frame(distance = NA, score = NA, percentage = NA, chr = NA)
for (chr in chosen_chr) {
    cat("-----------processing chr", chr, "------------", "\n", sep = "")
    
    # filtering
    temp_df <- filter(import_fil, chr1 == chr & chr2 == chr)
    temp_df <- arrange(temp_df, distance)
    # cat(nrow(temp_df), "\n")
    
    # scoring
    scores_temp <- add_scores(temp_df)
    scores_temp$percentage <- scores_temp$score / sum(scores_temp$score)
    scores_temp$chr <- rep(chr, nrow(scores_temp))
    scores <- rbind(scores, scores_temp)
}
scores <- scores[2:nrow(scores),]

# output
cat("--------------writing csv--------------", "\n")
write.csv(scores, paste0(path, "/", output_name, "_distance.csv"))