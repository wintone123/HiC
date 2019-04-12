# load library
library(dplyr)
library(tidyr)
library(parallel)
library(pheatmap)
library(RColorBrewer)
ulimit::memory_limit(4096)

# load function
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
                    scores <- rbind(scores, data.frame(name = data, score = n))
                    cond = FALSE
                }
            } else {
                s <- l 
                scores <- rbind(scores, data.frame(name = data, score = n))
                n <- 0
                break
            }

        }
    }
    scores <- scores[2:nrow(scores),]
    scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")
    return(scores)
}

matrxing <- function(chrB) {
    mat_temp <- mat

    # arrange
    bins_fil <- filter(bins_temp, chr2 == chrB)
    bins_fil <- arrange(bins_fil, chr1, start1, chr2, start2)
    bins_fil <- unite(bins_fil, temp, c(1:4), sep = "_")

    # scoring
    scores <- data.frame(name = NA, score = NA)
    s = 1
    n = 0
    cond = TRUE
    while (cond) {
        for (l in s:nrow(bins_fil)) {
            data <- bins_fil$temp[s]
            if (bins_fil$temp[l] == data) {
                n <- n + 1
                if (l == nrow(bins_fil)) {
                    scores <- rbind(scores, data.frame(name = data, score = n))
                    cond = FALSE
                }
            } else {
                s <- l 
                scores <- rbind(scores, data.frame(name = data, score = n))
                n <- 0
                break
            }
        }
    }   
    scores <- scores[2:nrow(scores),]
    scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")

    # matrixing
    scores <- unite(scores, chr1, c("chr1","start1"), sep = "_")
    scores <- unite(scores, chr2, c("chr2","start2"), sep = "_")
    for (k in 1:nrow(scores)) {
        chr1 <- scores$chr1[k]
        chr2 <- scores$chr2[k]
        mat_temp[chr1, chr2] <- scores$score[k] + mat_temp[chr1, chr2]
        mat_temp[chr2, chr1] <- scores$score[k] + mat_temp[chr2, chr1]
    }

    return(mat_temp)
}

# load info
path <- "/mnt/c/HiC/test4"
imput_csv <- "test_fil.csv"
output_name <- "chr1-X_1m"
bin_size <- 1000000
chosen_chr <- c(1:19,"X")
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("-------------loading csv-------------", "\n")
import <- readr::read_delim(file.path(path, imput_csv), delim = "\t")
# cat(object.size(import), "\n")

# load genome info
genome <- GenomeInfoDb::seqinfo(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
seqlength_list <- GenomeInfoDb::seqlengths(genome)[paste0("chr", chosen_chr)]
bins_list <- seqlength_list %/% bin_size + 1

# binning
cat("---------------binning---------------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)
# cat(object.size(bins), "\n")

# filtering
bins <- filter(bins, chr1 != chr2 | start1 != start2)

# making matrix
cat("------------making matrix------------", "\n")
xy_list <- vector()
for (i in 1:length(chosen_chr)) {
    xy_list <- append(xy_list, paste0(chosen_chr[i], "_", c(1:bins_list[i])))
}
mat <- matrix(0, length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list

# loop
cat("-------------loop start--------------", "\n")
mat_out <- mat
for (i in 1:length(chosen_chr)) {
    chrA <- chosen_chr[i]
    cat("-----------processing chr",chosen_chr[i], "-----------", "\n", sep = "")
    bins_temp <- filter(bins, chr1 == chrA)

    # 2 thread processing
    cl <- makeForkCluster(2)
    # clusterEvalQ(cl, c("dplyr", "tidyr"))
    clusterExport(cl, c("bins_temp", "bins_list", "chrA", "mat"))
    mat_list <- parLapply(cl, chosen_chr, matrxing)
    mat_temp <- Reduce("+", mat_list)
    stopCluster(cl)
    mat_out <- mat_out + mat_temp
}
cat("---------------loop end--------------", "\n")

# output matrix
# cat("-------------writing csv-------------", "\n")
# write.csv(mat_out, paste0(path, "/", output_name, "_mat.csv"))
# q()

# heatmapping
cat("-------------heatmapping-------------", "\n")
breaks <- c(0, 2^(0:max_legend(max(mat_out))))
# breaks <- seq(min(mat_out),max(mat_out),length.out = 256)
pic <- pheatmap(mat_out, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
                filename = paste0(path, "/", output_name, ".png"))
# cat("-------------writing pdf-------------", "\n")
# save_pheatmap_pdf(pic, paste0(path, "/", output_name, ".pdf"))