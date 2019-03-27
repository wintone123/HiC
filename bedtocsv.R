gff# import library
library(dplyr)
library(tidyr)
library(rtracklayer)

# load info
path <- "/mnt/c/Hic/test3"
import_bed <- "test.bed"
output_csv <- "test_fil.csv"
score=20

# import file
print("--------impoting bed......--------")
import <- import.bed(file.path(path, import_bed))

# filtering
print("---------filtering......---------")
import_df <- data.frame(chr = import@seqnames, start = start(import),
                        score = import$score)

# arrange and remove low quality pairs
print("---------arranging......---------")
output_df <- data.frame(chr1 = import_df[c(1:(nrow(import_df)/2))*2-1, "chr"], 
                        start1 = import_df[c(1:(nrow(import_df)/2))*2-1, "start"],
                        score1 = import_df[c(1:(nrow(import_df)/2))*2-1, "score"],
                        chr2 = import_df[c(1:(nrow(import_df)/2))*2, "chr"], 
                        start2 = import_df[c(1:(nrow(import_df)/2))*2, "start"],
                        score2 = import_df[c(1:(nrow(import_df)/2))*2, "score"],
                        stringsAsFactors = FALSE)
output_df <- filter(output_df, score1 >= score & score2 >= score)
output_df$score1 <- NULL
output_df$score2 <- NULL

# output file
print("--------writing csv......--------")
write.csv(output_df, file.path(path, output_csv))