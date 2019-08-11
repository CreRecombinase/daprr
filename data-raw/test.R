library(archive)
library(daprcpp)
library(tidyverse)
library(daprcpp)
gf <- system.file("gwas_z_t.txt.gz", package = "daprcpp")
af <- system.file("gwas_anno_t.txt.gz", package = "daprcpp")
gw_df <- read_delim(gf,delim = " ") %>% mutate(snp_id=1:n(),locus_id=as.integer(factor(locus)))  
anno_df <- read_tsv(af) %>% inner_join(select(gw_df,SNP,snp_id))
anno_f <- "/home/nwknoblauch//Dropbox/Repos/dap/torus_src/little_gwas_anno.tsv.zstd"
gwas_f <- "/home/nwknoblauch/Dropbox/Repos/dap/torus_src/little_gwas_i.tsv.zstd"
fr <- file_read(gwas_f)
gw_df <- read_tsv(fr)
# write_tsv(gw_df,"~/Downloads/dap/torus_src/gwas_z.txt.gz")
anno_df <- read_tsv(file_read(anno_f))
# write_tsv(anno_df,"~/Downloads/dap/torus_src/anno.txt.gz")


p <- 1000
o_anno_df <- tibble::tibble(SNP=1:p,aa=sample(0:1,p,replace=T),ab=sample(0:1,p,replace=T))
t_anno_m <- data.matrix(select(o_anno_df,-SNP)) 
dimnames(t_anno_m) <- NULL
t_anno_l <- list(annomat=t_anno_m,names=c("aa","ab"))
anno_df <- o_anno_df%>% tidyr::gather(key = "feature",value="value",-SNP) %>% dplyr::filter(value==1L)
annomat_l <- make_matrix(p,anno_df)
expect_equal(annomat_l,t_anno_l)

y <- runif(p)


new_df <- left_join(gw_df,anno_df) %>%
    mutate(Conserved_LindbladToh_d=replace_na(Conserved_LindbladToh_d,0))
anno_m <- data.matrix(select(new_df,Conserved_LindbladToh_d))
ret <- torus(locus_id = new_df$locus,z_hat = new_df$`z-stat`,anno_mat = anno_m,names = colnames(anno_m)
      )
