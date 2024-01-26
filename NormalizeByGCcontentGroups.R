library("dplyr")

args <- commandArgs(trailingOnly = TRUE)

groups_file <- args[2]
print(groups_file)


#file <- list.files(".", pattern = ".regions.bed.gz$")
file <- args[1]

moscov <- read.table(gzfile(file))
colnames(moscov) <- c("CHR", "start", "end", "coverage")

eid <- as.integer(gsub(".regions.bed.gz", "", file))

groups <-  read.table(groups_file, sep="\t")
colnames(groups) <- c("CHR", "start", "end", "chrY", "subclass")

norm_factor <-  moscov %>% merge(groups, by=c("CHR", "start", "end")) %>% 
  filter(chrY == 0) %>% group_by(subclass) %>% 
  summarize(median_norm=median(coverage))

result <- moscov %>% merge(groups, by=c("CHR", "start", "end"))   %>%
  filter(CHR == "chrY") %>% 
  merge(norm_factor, by="subclass", all.x=T) %>% 
  mutate(chrY_norm=coverage/median_norm) %>% 
  summarize(chrY_norm_median=median(chrY_norm)) %>% 
  mutate(eid=eid) %>% select(eid, chrY_norm_median)


write.table(result, paste0("./out_files/", eid, ".txt"), quote = F, row.names = F, sep="\t", col.names = F)