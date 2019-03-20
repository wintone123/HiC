# load library
library(dplyr)
library(tidyr)
ulimit::memory_limit(4096)

# load function
add_scores <- function(temp_df) {
    temp_df <- unite(temp_df, temp, c(1:4), sep = "_")
    scores <- data.frame(name = NA, score = NA)
    s = 1
    n = 0
    j <- 0
    k <- 0
    cond = TRUE
    while (cond) {
        for (l in s:nrow(temp_df)) {
            data <- temp_df$temp[s]
            if (temp_df$temp[l] == data) {
                n <- n + 1
                if (l == nrow(temp_df)) {
                    # cat(paste0(data,"--",n))
                    temp <- data.frame(name = data, score = n)
                    scores <- rbind(scores, temp)
                    cond = FALSE
                }
            } else {
                s <- l 
                # cat(paste0(data,"--",n))
                temp <- data.frame(name = data, score = n)
                scores <- rbind(scores, temp)
                n <- 0
                break
            }
            # progress bar
            if ((j / nrow(temp_df) * 100) %/% 1 == k) { 
                # cat("--------------------", k, "%--------------------", "\r", sep = "")
                cat("[", paste(rep("#", (37*k/100) %/% 1),collapse = ""), 
                paste(rep("_", 37-(37*k/100) %/% 1),collapse = ""), "]","\r", sep = "")
                k <- k + 1
            }
            j <- j + 1
        }
    }
    scores <- scores[2:nrow(scores),]
    scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")
    return(scores)
}

# load info
path <- "/mnt/c/HiC/test3"
imput_csv <- "test_fil.csv"
output_name <-"chr1_10k"
bin_size <- 200000
chosen_chr <- c(1)

# load file
cat("--------------loading csv--------------", "\n")
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)
# cat(object.size(import), "\n")

# binning
cat("----------------binning----------------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)
# cat(object.size(bins), "\n")

# loop
cat("---------------loop start--------------", "\n")
bins <- filter(bins, chr1 == chr2 & start1 != start2)
# cat(object.size(bins), "\n")
loops <- data.frame(length = NA, score = NA)
for (i in 1:length(chosen_chr)) {
    cat(paste0("-----------processing chr",chosen_chr[i], "-------------", "\n"))

    # filtering
    bins_fil <- filter(bins, chr1 == chosen_chr[i] & chr2 == chosen_chr[i])
    bins_fil <- arrange(bins_fil, chr1, start1, chr2, start2)
    # cat(object.size(bins_fil), "\n")

    # scoring
    scores <- add_scores(bins_fil)
    # cat(object.size(scores), "\n")
        
    # loops 
    for (j in 1:nrow(scores)) {
        l <- abs(as.numeric(scores$start2[j]) - as.numeric(scores$start1[j]))
        if (l %in% loops$length) {
            loops[loops$length == l, "score"] <- loops[loops$length == l, "score"] + scores$score[j]
        } else {
            loops <- rbind(loops, data.frame(length = l, score = scores$score[j]))
            if (i == 1 & j == 1) {
                loops <- loops[2:nrow(loops),]
            }
        }
    }
}
loops <- arrange(loops, length)
cat("----------------loop end---------------", "\n")

# output loops
cat("--------------writing csv--------------", "\n")
write.csv(loops, paste0(path, "/", output_name, "_loops.csv"))