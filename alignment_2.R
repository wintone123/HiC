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

# load info
path <- "/mnt/c/HiC/test1"
imput_csv <- "chr14_1.csv"
output_csv <- "chr14_40k.csv"
output_pdf <-  "chr14_40k.pdf"
output_png <- "chr14_40k.png"
bin_size <- 100000
chosen_chr1 <- c(14)
chosen_aera1 <- c(20000000,30000000)
chosen_chr2 <- c(14)
chosen_aera2 <- c(20000000,30000000)
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

# filtering
print("----------filtering......----------")
filt_temp <- data.frame(chr1 = NA, start1 = NA, chr2 = NA, start2 = NA)
if (mode == "SCSA") {
    for (i in length(chosen_chr1)) {
        for (j in length(chosen_chr2)) {
            filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr1[i] & start1 >= chosen_aera1[i] & start1 <= chosen_aera1[i*2] & 
                                                 chr2 == chosen_chr2[j] & start2 >= chosen_aera2[j] & start2 <= chosen_aera2[j*2]))
            filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr1[j] & start1 >= chosen_aera1[j] & start1 <= chosen_aera1[j*2] & 
                                                 chr2 == chosen_chr2[i] & start2 >= chosen_aera2[i] & start2 <= chosen_aera2[i*2]))
        }
    }
} else {
    filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr1[1] & start1 >= chosen_aera1[1] & start1 <= chosen_aera1[2] & 
                                         chr2 == chosen_chr2[1] & start2 >= chosen_aera2[1] & start2 <= chosen_aera2[2]))
    filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr1[1] & start1 >= chosen_aera1[1] & start1 <= chosen_aera1[2] & 
                                         chr2 == chosen_chr1[1] & start2 >= chosen_aera1[1] & start2 <= chosen_aera1[2]))
    filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr2[1] & start1 >= chosen_aera2[1] & start1 <= chosen_aera2[2] & 
                                         chr2 == chosen_chr1[1] & start2 >= chosen_aera1[1] & start2 <= chosen_aera1[2]))
    filt_temp <- rbind(filt_temp, filter(import, chr1 == chosen_chr2[1] & start1 >= chosen_aera2[1] & start1 <= chosen_aera2[2] & 
                                         chr2 == chosen_chr2[1] & start2 >= chosen_aera2[1] & start2 <= chosen_aera2[2]))                                 
}
filt_temp <- filt_temp[2:nrow(filt_temp),]

# binning
print("----------binning......----------")
bins <- data.frame(chr1 = as.character(filt_temp$chr1),
				   start1 = filt_temp$start1 %/% bin_size + 1,
				   chr2 = as.character(filt_temp$chr2),
				   start2 = filt_temp$start2 %/% bin_size + 1,
				   stringsAsFactors = FALSE)
bins <- arrange(bins, chr1, start1, chr2, start2)
bins <- filter(bins, chr1 != chr2 | start1 != start2)

# scoring
print("----------scoring......----------")
temp_df <- unite(bins, temp, c(1:4), sep = "_")
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
				temp <- data.frame(name = data, score = n)
				scores <- rbind(scores, temp)
				cond = FALSE
			}
		} else {
			s = i 
			temp <- data.frame(name = data, score = n)
			scores <- rbind(scores, temp)
			n = 0
 		    break
		}
	}
}
scores <- scores[2:nrow(scores),]
scores <- separate(scores, name, c("chr1","start1","chr2","start2"), sep = "_")

# output data
print("----------writing csv......----------")
write.csv(scores, file.path(path, output_csv))

# matrixing 
print("----------matrixing......----------")
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
breaks <- seq(min(scores$score), max(scores$score), length.out = 256)
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(256), breaks = breaks, legend = FALSE, border_color = NA,
                filename = file.path(path, output_png))

# output data
print("----------writing pdf......----------")
save_pheatmap_pdf(pic, file.path(path, output_pdf))