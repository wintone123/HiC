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
imput_csv <- "chr14.csv"
output_csv <- "chr14_200k.csv"
output_pdf <-  "chr14_200k.pdf"
bin_size <- 200000 
chosen_chr1 <- 14
chosen_start1 <- 30000000
chosen_end1 <- 67043718
chosen_chr2 <- 14
chosen_start2 <- 30000000
chosen_end2 <- 67043718

# import data
import <- read.delim(file.path(path, imput_csv), sep = ",", row.names = 1, header = TRUE)

# judgement
if (chosen_chr1 == chosen_chr2 & chosen_start1 == chosen_start2 & chosen_end1 == chosen_end2) {
	mode <- "same"
} else {
	mode <- "diff"
} 

# filtering
if (mode == "diff") {
	filt_temp1 <- filter(import, chr1 == chosen_chr1 & position1 >= chosen_start1 & position1 <= chosen_end1)
	filter_temp1_1 <- filter(filt_temp1, chr2 == chosen_chr1 & position2 >= chosen_start1 & position2 <= chosen_end1)
	filter_temp1_2 <- filter(filt_temp1, chr2 == chosen_chr2 & position2 >= chosen_start2 & position2 <= chosen_end2)
	filt_temp2 <- filter(import, chr2 == chosen_chr2 & position2 >= chosen_start2 & position2 <= chosen_end2)
	filter_temp2_1 <- filter(filt_temp2, chr1 == chosen_chr1 & position1 >= chosen_start1 & position1 <= chosen_end1)
	filter_temp2_2 <- filter(filt_temp2, chr1 == chosen_chr2 & position1 >= chosen_start2 & position1 <= chosen_end2)
	import_2 <- rbind(filter_temp1_1, filter_temp1_2, filter_temp2_1, filter_temp2_2)
} else {
	filt_temp <- filter(import, chr1 == chosen_chr1 & position1 >= chosen_start1 & position1 <= chosen_end1)
	import_2 <- filter(filt_temp, chr2 == chosen_chr2 & position2 >= chosen_start2 & position2 <= chosen_end2)
}

# binning
print("----------binning......----------")
output_1 <- data.frame(chr1 = import_2$chr1,
					   start1 = import_2$position1 %/% bin_size + 1,
					   chr2 = import_2$chr2,
					   start2 = import_2$position2 %/% bin_size + 1,
					   stringsAsFactors = FALSE)
output_1 <- arrange(output_1, chr1, start1, chr2, start2)
output_1 <- filter(output_1, chr1 != chr2 | start1 != start2)

# scoring
print("----------scoring......----------")
temp_df <- unite(output_1, temp, c(1:4), sep = "_")
output_2 <- data.frame(name = NA, score = NA)
s = 1
n = 0
cond = TRUE
while (cond) {
	for (i in s:nrow(temp_df)) {
		data <- temp_df$temp[s]
		if (temp_df$temp[i] == data) {
			n = n + 1
			if (i == nrow(temp_df)) {
				print(paste0(data,"---",n))
				temp <- data.frame(name = data, score = n)
				output_2 <- rbind(output_2, temp)
				cond = FALSE
			}
		} else {
			s = i 
			print(paste0(data,"---",n))
			temp <- data.frame(name = data, score = n)
			output_2 <- rbind(output_2, temp)
			n = 0
 		    break
		}
	}
}
output_2 <- output_2[2:nrow(output_2),]
output_2 <- separate(output_2, name, c("chr1","start1","chr2","start2"), sep = "_")

# output data
print("----------writing csv......----------")
write.csv(output_2, file.path(path, output_csv))

# matrixing 
print("----------matrixing......----------")
xy_list <- vector()
chr1_bin_end <- chosen_end1/ bin_size + 1
chr1_bin_start <- chosen_start1 / bin_size + 1
chr2_bin_end <- chosen_end2/ bin_size + 1
chr2_bin_start <- chosen_start2 / bin_size + 1
if (mode == "diff") {
	xy_list <- paste0(chosen_chr1,"_",c(chr1_bin_start:chr1_bin_end))
	xy_list <- append(xy_list, paste0(chosen_chr2,"_",c(chr2_bin_start:chr2_bin_end)))
} else {
	xy_list <- paste0(chosen_chr1,"_",c(chr1_bin_start:chr1_bin_end))
}
mat <- matrix(rep(0, length(xy_list)*length(xy_list)), length(xy_list), length(xy_list))
colnames(mat) <- xy_list
rownames(mat) <- xy_list
output_2 <- unite(output_2, chr1, c(1:2), sep = "_")
output_2 <- unite(output_2, chr2, c(2:3), sep = "_")
for (i in 1:nrow(output_2)) {
	chr1 <- output_2$chr1[i]
	chr2 <- output_2$chr2[i]
	mat[chr1,chr2] <- output_2$score[i] + mat[chr1,chr2]
	mat[chr2,chr1] <- output_2$score[i] + mat[chr2,chr1]
}

# heatmapping
breaks <- seq(min(output_2$score), max(output_2$score), length.out = 256)
col <- colorRampPalette(rev(brewer.pal(11,"RdYlBu")))
pic <- pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,
				show_rownames = FALSE, show_colnames = FALSE, 
				col = col(256), breaks = breaks, legend = FALSE)

# output data
print("----------writing pdf......----------")
save_pheatmap_pdf(pic, file.path(path, output_pdf))