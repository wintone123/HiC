# load library
library(dplyr)
library(tidyr)

# load info
path <- "/mnt/c/HiC/test3"
imput_csv <- "test_fil.csv"
output_name <-"chr1_2k_hub"
bin_size <- 2000
chosen_chr <- c(1)

# load file
cat("-------------loading file--------------", "\n")
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

# binning & filtering
cat("----------binning & filtering----------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)
bins <- filter(bins, chr1 %in% chosen_chr & chr2 %in% chosen_chr &
               chr1 == chr2 & start1 != start2)
bins <- arrange(bins, chr1, chr2, start1,start2)
bins <- unite(bins, chr1, c("chr1","start1"), sep = "_")
bins <- unite(bins, chr2, c("chr2","start2"), sep = "_")

# scoring
cat("----------------scoring----------------", "\n")
# k <- 0
# for (i in 1:nrow(bins)) {
#     p1 <- bins$chr1[i]
#     p2 <- bins$chr2[i]
#     if (i == 1) {
#         scores <- data.frame(position = p1, score = 1)
#         scores <- rbind(scores, data.frame(position = p2, score = 1))
#     } else {
#         if (p1 %in% scores$position) {
#             scores[scores$position == p1, "score"] = scores[scores$position == p1, "score"] + 1    
#         } else {
#             scores <- rbind(scores, data.frame(position = p1, score = 1))
#         }
#         if (p2 %in% scores$position) {
#             scores[scores$position == p2, "score"] = scores[scores$position == p2, "score"] + 1
#         } else {
#             scores <- rbind(scores, data.frame(position = p2, score = 1))
#         }
#     }
#     # progress bar
#     if ((i / nrow(bins) * 100) %/% 1 == k) { 
#         # cat("--------------------", k, "%--------------------", "\r", sep = "")
#         cat("[", paste(rep("#", (37*k/100) %/% 1), collapse = ""), 
#         paste(rep("_", 37-(37*k/100) %/% 1),collapse = ""), "]","\r", sep = "")
#         k <- k + 1
#     }
# }

p_list <- unique(c(bins$chr1, bins$chr2))
# cat(length(p_list), "\n")
temp <- c(bins$chr1, bins$chr2)
scores <- data.frame(position = p_list, score = rep(0, length(p_list)))
k <- 0
for (i in 1:length(p_list)) {
    p <- p_list[i]           
    num <- sum((temp == p) == TRUE)
    temp <- temp[temp != p]
    scores[scores$position == p, "score"] <- scores[scores$position == p, "score"] + num
    # progress bar
    if ((i / length(p_list) * 100) %/% 1 == k) { 
        # cat("--------------------", k, "%--------------------", "\r", sep = "")
        cat("[", paste(rep("#", (37*k/100) %/% 1), collapse = ""), 
        paste(rep("_", 37-(37*k/100) %/% 1),collapse = ""), "]","\r", sep = "")
        k <- k + 1
    }
}

scores <- separate(scores, position, c("chr", "position"), sep = "_")
scores <- arrange(scores, chr, desc(score))

# output data
cat("--------------writing csv--------------", "\n")
write.csv(scores, paste0(path, "/", output_name, ".csv"))