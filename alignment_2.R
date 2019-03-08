# load library
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

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

# load info
path <- "/mnt/c/HiC/test1"
imput_csv <- "chr14_1.csv"
output_csv <- "chr14_40k.csv"
output_mat <- "chr14_40k_mat.csv"
output_pdf <-  "chr14_40k.pdf"
output_png <- "chr14_5k.png"
bin_size <- 5000
chosen_chr1 <- c(14)
chosen_aera1 <- c(70000000,71000000)
chosen_chr2 <- c(14)
chosen_aera2 <- c(70000000,71000000)
col <- colorRampPalette(brewer.pal(9,"YlOrRd"))

# load file
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

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
print("------------binning......------------")
bins <- data.frame(chr1 = as.character(import$chr1),
			       start1 = import$start1 %/% bin_size + 1,
			       chr2 = as.character(import$chr2),
			       start2 = import$start2 %/% bin_size + 1,
			       stringsAsFactors = FALSE)

# filtering
print("-----------filtering......-----------")
chosen_bin1 <- chosen_aera1 %/% bin_size + 1
chosen_bin2 <- chosen_aera2 %/% bin_size + 1
if (mode == "SCSA") {
    bins_fil <- filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                       chr2 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2])
} else {
    bins_fil <- data.frame(chr1 = NA, start1 = NA, chr2 = NA, start2 = NA)
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                                       chr1 == chosen_chr1[1] & start2 >= chosen_bin1[1] & start2 <= chosen_bin1[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr1[1] & start1 >= chosen_bin1[1] & start1 <= chosen_bin1[2] & 
                                       chr1 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr2[1] & start1 >= chosen_bin2[1] & start1 <= chosen_bin2[2] & 
                                       chr1 == chosen_chr1[1] & start2 >= chosen_bin1[1] & start2 <= chosen_bin1[2]))
    bins_fil <- rbind(bins_fil, filter(bins, chr1 == chosen_chr2[1] & start1 >= chosen_bin2[1] & start1 <= chosen_bin2[2] & 
                                       chr1 == chosen_chr2[1] & start2 >= chosen_bin2[1] & start2 <= chosen_bin2[2]))
    bins_fil <- bins_fil[2:nrow(bins_fil),]
} 
if (mode != "DC") { # remove same bins
    bins_fil <- filter(bins_fil, start1 != start2) 
} else {
    bins_fil <- filter(bins_fil, chr1 != chr2 & start1 != start2)
} 
if (mode != "DC") {
    bins_fil2 <- data.frame(chr1 = bins_fil$chr1,
                            start1 = ifelse(bins_fil$start1 <= bins_fil$start2, bins_fil$start1, bins_fil$start2),
                            chr2 = bins_fil$chr2,
                            start2 = ifelse(bins_fil$start1 <= bins_fil$start2, bins_fil$start2, bins_fil$start1),
                            stringsAsFactors = FALSE)
}
bins_fil2 <- arrange(bins_fil2, chr1, start1, chr2, start2)

# scoring
print("------------scoring......------------")
temp_df <- unite(bins_fil2, temp, c(1:4), sep = "_")
scores <- data.frame(name = NA, score = NA)
s = 1
n = 0
cond = TRUE
while (cond) {
	for (i in s:nrow(temp_df)) {
		data <- temp_df$temp[s]
		if (temp_df$temp[i] == data) {
			n = n + 1
			if (i == nrow(temp_df)) {
                # print(paste0(data,"--",n))
				temp <- data.frame(name = data, score = n)
				scores <- rbind(scores, temp)
				cond = FALSE
			}
		} else {
			s = i 
            # print(paste0(data,"--",n))
			temp <- data.frame(name = data, score = n)
			scores <- rbind(scores, temp)
			n = 0
 		    break
		}
	}
}
scores <- scores[2:nrow(scores),]
scores <- filter(scores, score >= 1) # remove fewer than 1 reads
scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")

# output data
print("----------writing csv......----------")
write.csv(scores, file.path(path, output_csv))

# matrixing 
print("-----------matrixing......-----------")
xy_list <- vector()
if (mode == "SCSA") {
    for (i in length(chosen_chr1)){
        bin_start <- chosen_aera1[i] %/% bin_size + 1
        bin_end <- chosen_aera1[i*2] %/% bin_size + 1
        xy_list <- append(xy_list, paste0(chosen_chr1[i],"_",c(bin_start:bin_end)))
    }
} else {
    bin_start_1 <- chosen_aera1[1] %/% bin_size + 1
    bin_end_1 <- chosen_aera1[2] %/% bin_size + 1
    xy_list <- paste0(chosen_chr1[1],"_",c(bin_start_1:bin_end_1))
    bin_start_2 <- chosen_aera2[1] %/% bin_size + 1
    bin_end_2 <- chosen_aera2[2] %/% bin_size + 1
    xy_list <- append(xy_list, paste0(chosen_chr2[1],"_",c(bin_start_2:bin_end_2)))
} 
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list
scores <- unite(scores, chr1, c(1:2), sep = "_")
scores <- unite(scores, chr2, c(2:3), sep = "_")
for (i in 1:nrow(scores)) {
	chr1 <- scores$chr1[i]
	chr2 <- scores$chr2[i]
	mat[chr1,chr2] <- scores$score[i] + mat[chr1,chr2]
	mat[chr2,chr1] <- scores$score[i] + mat[chr2,chr1]
}

# heatmapping
breaks <- c(0,2^(0:max_legend(max(mat))))
# breaks <- seq(min(mat),max(mat),length.out = 256)
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(length(breaks)), breaks = breaks, legend = FALSE, border_color = NA,
                filename = file.path(path, output_png))

# output mat
# output mat
print("----------writing csv......----------")
write.csv(mat, file.path(path, output_mat))
print("----------writing pdf......----------")
save_pheatmap_pdf(pic, file.path(path, output_pdf))