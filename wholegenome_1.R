# load library
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
ulimit::memory_limit(4096)

# load function
# save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
#    stopifnot(!missing(x))
#    stopifnot(!missing(filename))
#    pdf(filename, width=width, height=height)
#    grid::grid.newpage()
#    grid::grid.draw(x$gtable)
#    dev.off()
# }

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
                    scores <- rbind(scores, data.frame(name = data, score = n))
                    cond = FALSE
                }
            } else {
                s <- l 
                # cat(paste0(data,"--",n)
                scores <- rbind(scores, data.frame(name = data, score = n))
                n <- 0
                break
            }
            # progress bar
            if ((j / nrow(temp_df) * 100) %/% 1 == k) { 
                cat("[", paste(rep("#", (33*k/100) %/% 1),collapse = ""), 
                paste(rep("_", 33-(33*k/100) %/% 1),collapse = ""), k, "%]","\r", sep = "")
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
path <- "/mnt/c/HiC/test4"
imput_csv <- "test_fil.csv"
output_name <-"chr1-X_1m"
bin_size <- 1000000
chosen_chr <- c(1:19,"X")
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("--------------loading csv--------------", "\n")
import <- readr::read_delim(file.path(path, imput_csv), delim = "\t")
# cat(object.size(import), "\n")

# load genome info
genome <- GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
seqlength_list <- GenomeInfoDb::seqlengths(genome)[paste0("chr", chosen_chr)]
bins_list <- seqlength_list %/% bin_size + 1

# binning
cat("----------------binning----------------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)
# cat(object.size(bins), "\n")

#filtering
bins <- filter(bins, chr1 != chr2)

# making matrix
cat("-------------making matrix-------------", "\n")
xy_list <- vector()
for (i in 1:length(chosen_chr)) {
    xy_list <- append(xy_list, paste0(chosen_chr[i], "_", c(1:bins_list[i])))
}
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list

# loop
cat("---------------loop start--------------", "\n")
# cat(object.size(bins), "\n")
for (i in 1:length(chosen_chr)) {
    for (j in 1:length(chosen_chr)) {
        chrA <- chosen_chr[i]
        chrB <- chosen_chr[j]
        if (chrA != chrB) {
            cat("---------processing chr", chrA, "-chr", chrB, "----------", "\n", sep = "")

            # filtering
            bins_fil <- filter(bins, chr1 == chrA & start1 >= 1 & start1 <= bins_list[i] &
                            chr2 == chrB & start2 >= 1 & start2 <= bins_list[j])
            bins_fil <- arrange(bins_fil, chr1, start1, chr2, start2)
            # cat(object.size(bins_fil))

            # scoring
            scores <- add_scores(bins_fil)
            # cat(object.size(scores), "\n")

            # matrixing
            scores <- unite(scores, chr1, c("chr1","start1"), sep = "_")
            scores <- unite(scores, chr2, c("chr2","start2"), sep = "_")
            for (k in 1:nrow(scores)) {
                chr1 <- scores$chr1[k]
                chr2 <- scores$chr2[k]
                mat[chr1, chr2] <- scores$score[k] + mat[chr1, chr2]
                mat[chr2, chr1] <- scores$score[k] + mat[chr2, chr1]
            }
            # cat(object.size(mat), "\n")
        }
    }
}
cat("----------------loop end---------------", "\n")

# output matrix
# cat("--------------writing csv--------------", "\n")
# write.csv(mat, paste0(path, "/", output_name, "_mat.csv"))

# heatmapping
cat("--------------heatmapping--------------", "\n")
breaks <- c(0, 2^(0:max_legend(max(mat))))
# breaks <- seq(min(mat),max(mat),length.out = 256)
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
                filename = paste0(path, "/", output_name, ".png"))
# cat("--------------writing pdf--------------", "\n")
# save_pheatmap_pdf(pic, paste0(path, "/", output_name, ".pdf"))