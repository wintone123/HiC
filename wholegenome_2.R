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

# load info
path <- "/mnt/c/HiC/test4"
imput_csv <- "test_fil.csv"
output_name <-"chr1-X"
chosen_chr <- c(1:19,"X")
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("--------------loading csv--------------", "\n")
import <- readr::read_delim(file.path(path, imput_csv), delim = "\t")
# cat(object.size(import), "\n")

# filtering
import <- filter(import, chr1 != chr2)
import <- arrange(import, chr1, chr2)

# making matrix
cat("-------------making matrix-------------", "\n")
xy_list <- vector()
mat <- matrix(rep(0, length(chosen_chr)*length(chosen_chr)), length(chosen_chr), length(chosen_chr))
colnames(mat) <- chosen_chr
rownames(mat) <- chosen_chr

# loop
cat("---------------loop start--------------", "\n")
# cat(object.size(bins), "\n")
for (i in 1:length(chosen_chr)) {
    for (j in 1:length(chosen_chr)) {
        a <- chosen_chr[i]
        b <- chosen_chr[j]
        if (a != b) {
            cat("---------processing chr", a, "-chr", b, "---------", "\r", sep = "")
 
            # scoring & matrixing
            mat[a, b] <- mat[a, b] + nrow(filter(import, chr1 == a & chr2 == b))
            mat[b, a] <- mat[b, a] + nrow(filter(import, chr1 == a & chr2 == b))
        }
    }
}
cat("----------------loop end---------------", "\n")

# output matrix
# cat("--------------writing csv--------------", "\n")
# write.csv(mat, paste0(path, "/", output_name, "_mat.csv"))

# heatmapping
cat("--------------heatmapping--------------", "\n")
# breaks <- c(0, 2^(0:max_legend(max(mat))))
breaks <- seq(min(mat),max(mat),length.out = 100)
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
                filename = paste0(path, "/", output_name, ".png"))
# cat("--------------writing pdf--------------", "\n")
# save_pheatmap_pdf(pic, paste0(path, "/", output_name, ".pdf"))