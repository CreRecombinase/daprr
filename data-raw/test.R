library(archive)
library(daprcpp)
library(tidyverse)
library(daprcpp)
anno_f <- "/home/nwknoblauch//Dropbox/Repos/dap/torus_src/little_gwas_anno.tsv.zstd"
gwas_f <- "/home/nwknoblauch/Dropbox/Repos/dap/torus_src/little_gwas_i.tsv.zstd"
fr <- file_read(gwas_f)
gw_df <- read_tsv(fr)
# write_tsv(gw_df,"~/Downloads/dap/torus_src/gwas_z.txt.gz")
anno_df <- read_tsv(file_read(anno_f))
# write_tsv(anno_df,"~/Downloads/dap/torus_src/anno.txt.gz")
new_df <- left_join(gw_df,anno_df) %>%
    mutate(Conserved_LindbladToh_d=replace_na(Conserved_LindbladToh_d,0))
anno_m <- data.matrix(select(new_df,Conserved_LindbladToh_d))
ret <- torus(locus_id = new_df$locus,z_hat = new_df$`z-stat`,anno_mat = anno_m,names = colnames(anno_m)
      )
