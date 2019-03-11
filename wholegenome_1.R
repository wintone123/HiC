# load library
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(BSgenome.Hsapiens.UCSC.hg38)

#load function
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

max_legend <- function(a) {
    cond <- TRUE
    n <- 0
    while (cond) {
        n <- n + 1
        b = a %/% 2
        if (b == 1) {
            cond <- FALSE
            return(n + 1)
        } else {
            a <- b
        }
    }
}

add_scores <- function(temp_df) {
    temp_df <- unite(temp_df, temp, c(1:4), sep = "_")
    scores <- data.frame(name = NA, score = NA)
    s = 1
    n = 0
    cond = TRUE
    while (cond) {
        for (l in s:nrow(temp_df)) {
            data <- temp_df$temp[s]
            if (temp_df$temp[l] == data) {
                n <- n + 1
                if (l == nrow(temp_df)) {
                    # print(paste0(data,"--",n))
                    temp <- data.frame(name = data, score = n)
                    scores <- rbind(scores, temp)
                    cond = FALSE
                }
            } else {
                s <- l 
                # print(paste0(data,"--",n))
                temp <- data.frame(name = data, score = n)
                scores <- rbind(scores, temp)
                n <- 0
                break
            }
        }
    }
    scores <- scores[2:nrow(scores),]
    scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")
    return(scores)
}

# load info
path <- "/mnt/c/HiC/test3"
imput_csv <- "test_fil.csv"
output_name <-"chr1-X_5m"
bin_size <- 5000000
chosen_chr <- c(1:22,"X")
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
print("----------loading csv......----------")
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

# load genome info
genome <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
seqlength_list <- seqlengths(genome)[paste0("chr", chosen_chr)]
bins_list <- seqlength_list %/% bin_size + 1

# binning
print("------------binning......------------")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)

# making matrix
print("---------making matrix......---------")
xy_list <- vector()
for (i in 1:length(chosen_chr)) {
    xy_list <- append(xy_list, paste0(chosen_chr[i], "_", c(1:bins_list[i])))
}
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list

# loop
print("-------------loop start-------------")
bins <- filter(bins, chr1 != chr2 | start1 != start2)
for (i in 1:length(chosen_chr)) {
    for (j in 1:length(chosen_chr)) {
        print(paste0("--------proceesing chr",chosen_chr[i], "-chr", chosen_chr[j], "--------"))

        # filtering
        bins_fil <- filter(bins, chr1 == chosen_chr[i] & start1 >= 1 & start1 <= bins_list[i] &
                            chr2 == chosen_chr[j] & start2 >= 1 & start2 <= bins_list[j])
        bins_fil <- arrange(bins_fil, chr1, start1, chr2, start2)

        # scoring
        scores <- add_scores(bins_fil)

        # matrixing
        scores <- unite(scores, chr1, c("chr1","start1"), sep = "_")
        scores <- unite(scores, chr2, c("chr2","start2"), sep = "_")
        for (k in 1:nrow(scores)) {
            chr1 <- scores$chr1[k]
            chr2 <- scores$chr2[k]
            mat[chr1, chr2] <- scores$score[k] + mat[chr1, chr2]
            mat[chr2, chr1] <- scores$score[k] + mat[chr2, chr1]
        }
    }
}
print("--------------loop end--------------")

# heatmapping
print("----------heatmapping......----------")
breaks <- c(0, 2^(0:max_legend(max(mat))))
# breaks <- seq(min(mat),max(mat),length.out = 256)
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
                filename = paste0(path, "/", output_name, ".png"))
# print("----------writing pdf......----------")
# save_pheatmap_pdf(pic, paste0(path, "/", output_name, ".pdf"))