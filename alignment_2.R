# load library
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

#load function
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
imput_csv <- "chr1.csv"
output_name <-"chr1_1m"
bin_size <- 1000000
chosen_chr1 <- c(1)
chosen_aera1 <- c(1,195000000)
chosen_chr2 <- c(1)
chosen_aera2 <- c(1,195000000)
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
cat("--------------loading csv--------------", "\n")
import <- readr::read_delim(file.path(path, imput_csv), delim = "\t")

# judgement
if (all(chosen_chr1 == chosen_chr2)) {
    if (all(chosen_aera1 == chosen_aera2)) {
        mode <- "SCSA" # same_chr_same_area
    } else {
        mode <- "SCDA" # same_chr_diff_area
    }
} else {
    mode <- "DC" # diff_chr
}

# binning
cat("----------------binning----------------", "\n")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)

# filtering
cat("---------------filtering---------------", "\n")
chosen_bin1 <- chosen_aera1 %/% bin_size + 1
chosen_bin2 <- chosen_aera2 %/% bin_size + 1
if (mode == "SCSA") {
    bins_fil <- filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                       chr2 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2])
} else {
    bins_fil <- data.frame(chr1 = NA, start1 = NA, chr2 = NA, start2 = NA)
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                                       chr2 == chosen_chr1[1] & start2 >= chosen_bin1[1] & start2 <= chosen_bin1[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                                       chr2 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr2[1] & start1 >= chosen_bin2[1] & start1 <= chosen_bin2[2] & 
                                       chr2 == chosen_chr1[1] & start2 >= chosen_bin1[1] & start2 <= chosen_bin1[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr2[1] & start1 >= chosen_bin2[1] & start1 <= chosen_bin2[2] & 
                                       chr2 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2]))
    bins_fil <- bins_fil[2:nrow(bins_fil),]
} 
if (mode != "DC") { # remove same bins
    bins_fil <- filter(bins_fil, start1 != start2) 
    bins_fil2 <- data.frame(chr1 = bins_fil$chr1,
                            start1 = ifelse(bins_fil$start1 <= bins_fil$start2, bins_fil$start1, bins_fil$start2),
                            chr2 = bins_fil$chr2,
                            start2 = ifelse(bins_fil$start1 <= bins_fil$start2, bins_fil$start2, bins_fil$start1),
                            stringsAsFactors = FALSE)
} else {
    bins_fil2 <- filter(bins_fil, chr1 != chr2 | start1 != start2)
} 
bins_fil2 <- arrange(bins_fil2, chr1, start1, chr2, start2)
# cat(nrow(bins_fil2), "\n")

# output data
# cat("--------------writing csv--------------", "\n")
# write.csv(bins_fil2, paste0(path, "/", output_name, ".csv"))
# q()

# scoring
cat("----------------scoring----------------", "\n")
temp_df <- unite(bins_fil2, temp, c(1:4), sep = "_")
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
# scores <- filter(scores, score >= 1) # remove fewer than 1 reads
scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")

# output data
# cat("--------------writing csv--------------", "\n")
# write.csv(scores, paste0(path, "/", output_name, "_scores.csv"))

# matrixing 
cat("---------------matrixing---------------", "\n")
xy_list <- vector()
if (mode == "SCSA") {
    bin_start <- chosen_aera1[1] %/% bin_size + 1
    bin_end <- chosen_aera1[2] %/% bin_size + 1
    xy_list <- append(xy_list, paste0(chosen_chr1[1], "_", c(bin_start:bin_end)))
} else {
    bin_start_1 <- chosen_aera1[1] %/% bin_size + 1
    bin_end_1 <- chosen_aera1[2] %/% bin_size + 1
    xy_list <- paste0(chosen_chr1[1],"_",c(bin_start_1:bin_end_1))
    bin_start_2 <- chosen_aera2[1] %/% bin_size + 1
    bin_end_2 <- chosen_aera2[2] %/% bin_size + 1
    xy_list <- append(xy_list, paste0(chosen_chr2[1], "_", c(bin_start_2:bin_end_2)))
} 
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list
scores <- unite(scores, chr1, c(1:2), sep = "_")
scores <- unite(scores, chr2, c(2:3), sep = "_")
for (i in 1:nrow(scores)) {
	chr1 <- scores$chr1[i]
	chr2 <- scores$chr2[i]
	mat[chr1, chr2] <- scores$score[i] + mat[chr1, chr2]
	mat[chr2, chr1] <- scores$score[i] + mat[chr2, chr1]
}

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