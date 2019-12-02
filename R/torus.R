write_anno <- function(anno_df=tibble(SNP=integer(),feature=character()),p=max(max(anno_df$SNP),1L),af=tempfile(fileext=".txt.gz")){
  if(is.null(anno_df)){
    return(write_anno(af=af))
  }
  spread_anno_l <- make_matrix(p = p,anno_df = anno_df)
  spread_anno_df <- tibble::as_tibble(magrittr::set_colnames(spread_anno_l$annomat,paste0(spread_anno_l$names,"_d"))) %>%
    mutate(SNP=1:n()) %>% dplyr::select(SNP,dplyr::everything()) %>% filter_at(.vars = vars(-SNP),any_vars(. != 0))
  readr::write_tsv(spread_anno_df,path=af)
  return(af)
}




#' Title
#'
#' @param X 
#' @param gw_df 
#' @param alpha 
#' @param lambda 
#' @param EM_thresh 
#' @param prior_init 
#' @param parallel 
#'
#' @return
#' @export
#'
#' @examples
torus_glmnet <- function(X,z,region_id,nfolds=10,alpha=0.95,lambda=NULL,EM_thresh=0.05,prior_init=1e-3,...){
  snpnum <- length(z)
  uf_id <- sample(unique(region_id))
  n_reg <- length(uf_id)
  uid <- ((1:n_reg) %% nfolds)+1
  names(uid) <- uf_id
  folds <- uid[region_id]
  stopifnot(snpnum==nrow(X),
            snpnum==length(region_id))
  
  prior <- rep(prior_init,snpnum)
  BF <- daprcpp:::zhat2BF(z)
  splt <- daprcpp:::new_splitter(as.integer(factor(gw_df$region_id)))
  old_lik <- -9999
  new_lik <- 0
  it <- 1
  while(abs(new_lik-old_lik)>EM_thresh){
    ES_o <- daprcpp:::Esteps(BF,prior,splt)
    new_lik <- attr(ES_o,"lik")
    ym <- cbind(1-ES_o,ES_o)
    mg <- glmnet::cv.glmnet(x = X,y=ym,family="binomial",alpha=alpha,lambda=lambda,foldid = folds,...)
    a_prior <- map(mg$lambda,
                   ~c(glmnet::predict.cv.glmnet(object = mg,newx = X,type="response",s=.x))
                   )
    a_ES <- map(a_prior,~daprcpp:::Esteps(BF,.x,splt))
    a_lik <- map_dbl(a_ES,~attr(.x,"lik"))
    prior <- c(glmnet::predict.cv.glmnet(object = mg,newx = X,type="response",s="lambda.min"))
    ES <- daprcpp:::Esteps(BF,prior,splt)
    old_lik <- new_lik
    new_lik <- attr(ES,"lik")
    stopifnot(new_lik>old_lik)
    if(it>20){
      break
    }
    it <- it+1
    cat(it,"\n")
  }
 return(mg)
}



fps <- function(betav,X,BF,splt){
  prior <- 1/(1+exp(-c(as.matrix(betav[1]+X%*%betav[-1]))))
  ES_o <- daprcpp:::Esteps(BF,prior,splt)
  ym <- cbind(1-ES_o,ES_o)
  mg <- glmnet::cv.glmnet(x = X,
                          y=ym,
                          family="binomial",
                          alpha=0.95,
                          parallel = TRUE,
                          lambda=10^(-seq(from=5,to=10,length.out=100)))
  return(c(as.matrix(coef(mg))))
}

#' Title
#'
#' @param s 
#' @param lik 
#'
#' @return
#' @export
#'
#' @examples
parse_torus_s <- function(s,lik) {
    df <- read.table(file = textConnection(s),skip=1,header=F,stringsAsFactors = F)
    colnames(df) <- c("term", "estimate", "low", "high")

    df <- dplyr::mutate(df,term=stringr::str_replace(term,pattern = "\\.[0-9]+$",replacement = ""),
                        sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
    # lik <- scan(lik_file,what=numeric())
    # file.remove(lik_file)mutate
    df <- tidyr::nest(df,data = tidyr::everything()) %>% dplyr::mutate(lik=lik)
}

write_gwas <- function(gwas_df,gf=tempfile(fileext=".txt.gz")){
  
  dplyr::select(gwas_df,SNP,region_id,`z-stat`) %>% write_tsv(path=gf)
  return(gf)
}


gen_torus_cmd <- function(gf,af,torus_p=character(0)){
  torus_d <- fs::file_temp()
  lik_file <- fs::file_temp()
  if(length(torus_p)>0){
    p_f <- fs::path(torus_d,torus_p,ext="prior")
    stopifnot(!fs::dir_exists(torus_d))
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est",
      "-dump_prior",
      torus_d)
  } else{
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est"
    )
  }
  return(res_args)
}

##' Run torus with torus-formatted data
##'
##'
##'
##' @param gf gwas summary stats file
##' @param af annotation file
##' @param torus_p regions at which to estimate priors
##' @return a torus result object
##' @author
run_torus_cmd <- function(gf,af,torus_p=character(0)){

  torus_path <- system.file("dap-master/torus_src/torus",package = "daprcpp")
  stopifnot(file.exists(torus_path),torus_path!="")
  fo <- fs::file_info(torus_path)
  stopifnot((fo$permissions & "u+x") == "u+x")
  torus_d <- fs::file_temp()
  lik_file <- fs::file_temp()
  if(length(torus_p)>0){
    p_f <- fs::path(torus_d,torus_p,ext="prior")
    stopifnot(!fs::dir_exists(torus_d))
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est",
      "-dump_prior",
      torus_d)
  } else{
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est"
    )
  }
  res <- processx::run(torus_path,args = res_args,echo_cmd = TRUE)
  df <- read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F)
  colnames(df) <- c("term", "estimate", "low", "high")

  df <- dplyr::mutate(df,term=stringr::str_replace(term,pattern = "\\.[0-9]+$",replacement = ""),
                      sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  lik <- scan(lik_file,what=numeric())
  file.remove(lik_file)
  df <- tidyr::nest(df,data = tidyr::everything()) %>% dplyr::mutate(lik=lik)
  if(length(torus_p)>0){
    stopifnot(all(fs::file_exists(p_f)))
    prior_l <- map(torus_p,function(x){
      fp <- as.character(fs::path(torus_d,x,ext="prior"))
      suppressMessages(
        ret <- vroom::vroom(file = fp,delim = "  ",trim_ws = T,col_names = c("SNP","prior"),col_types = cols("SNP"="c","prior"="d")) %>% mutate(region_id=x)
      )
      return(ret)
    })
    fs::file_delete(p_f)

    names(prior_l) <- torus_p
    return(list(df=df,priors=prior_l))
  }else{
    return(list(df=df))
  }
}



test_torus_cmd <- function(gf,af,torus_p=character(0)){
  torus_path <- system.file("dap-master/torus_src/torus",package = "daprcpp")
  stopifnot(file.exists(torus_path),torus_path!="")
  fo <- fs::file_info(torus_path)
  stopifnot((fo$permissions & "u+x") == "u+x")
  torus_d <- fs::file_temp()
  lik_file <- fs::file_temp()
  if(length(torus_p)>0){
    p_f <- fs::path(torus_d,torus_p,ext="prior")
    stopifnot(!fs::dir_exists(torus_d))
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est",
      "-dump_prior",
      torus_d)
  } else{
    res_args <- c(
      "-d",
      fs::path_expand(gf),
      "-annot",
      fs::path_expand(af),
      "--load_zval",
      "-lik",
      lik_file,
      "-est"
    )
  }
  res <- dap_torus(res_args)
  df <- read.table(file = textConnection(nret),skip=1,header=F,stringsAsFactors = F)
  colnames(df) <- c("term", "estimate", "low", "high")
  
  df <- dplyr::mutate(df,term=stringr::str_replace(term,pattern = "\\.[0-9]+$",replacement = ""),
                      sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  lik <- scan(lik_file,what=numeric())
  file.remove(lik_file)
  df <- tidyr::nest(df,data = tidyr::everything()) %>% dplyr::mutate(lik=lik)
  if(length(torus_p)>0){
    stopifnot(all(fs::file_exists(p_f)))
    prior_l <- map(torus_p,function(x){
      fp <- as.character(fs::path(torus_d,x,ext="prior"))
      suppressMessages(
        ret <- vroom::vroom(file = fp,delim = "  ",trim_ws = T,col_names = c("SNP","prior"),col_types = cols("SNP"="i","prior"="d")) %>% mutate(region_id=x)
      )
      return(ret)
    })
    fs::file_delete(p_f)
    names(prior_l) <- torus_p
    return(list(df=df,priors=prior_l))
  }else{
    return(list(df=df))
  }
}





coef.torus <- function(fit){
  ret <- tidyr::unnest(fit$df)
  retvec <- ret$estimate
  names(retvec) <- ret$term
  return(retvec)
}




#' Perform a 'forward selection'
#'
#' @param f A function that fits the underlying model, taking `X` (and anything passed to `...`) as argument.
#'  `f` should return a list with a `lik` element, which will be used for the forward selection.
#' @param params a list over which to perform the forward selection.

#' @param combo_fun a function for generating input to `f` by combining elements of `params`.  `combo_fun` should
#' accept subsets of `params` that are of length 0(for the NULL model) or greater (for subsequent cases). At each step of 
#' the forward selection, each remaining element (not in the current model) of `params` will be added to the model
#' using `combo_fun`,the results of which will be passed to `f`
#' @param extract_terms a function for extracting the relevant terms from the model, so the params[extract_terms]
#' 
#' @param 
#' @return
#' @export
#'
#' @examples
forward_select_fun <- function(f,params,combo_fun,extract_terms,steps=1L,ret_all=FALSE){
  
  rest_terms <- params
  term_selection <- list(character())
  all_results <- list()
  
  db_fun <- function(ts){
    init_d <- combo_fun(ts)
    f(init_d)
  }
  
  for(i in seq_len(steps+1L)){
    all_fit <- purrr::map(term_selection,db_fun)
    all_results[[i]] <- all_fit
    lik_vec <- map_dbl(all_fit,~.x$df$lik)
    best_fit <- all_fit[[which.max(lik_vec)]]
    best_terms <- extract_terms(best_fit)
    stopifnot(all(best_terms %in% params))
    rest_terms <- rest_terms[! rest_terms %in% best_terms]
    term_selection <- purrr::map(rest_terms,~c(.x,best_terms))
  }
  if(ret_all){
    return (all_results)
  }else{
    return(best_fit) 
  }
}



#' Get locus-level FDR from torus
#'
#' @param gf gwas summary statistics
#' @param af annotation file
#'
#' @return torus result object
#' @export
#'
#' @examples
torus_fdr <- function(gf,af){
  torus_path <- system.file("dap-master/torus_src/torus",package = "daprcpp")
  stopifnot(file.exists(torus_path),torus_path!="")
  stopifnot(file.exists(gf),
            file.exists(af))
  fo <- fs::file_info(torus_path)
  stopifnot((fo$permissions & "u+x") == "u+x")
  qtl_file <- fs::file_temp()
  res_args <- c(
    "-d",
    fs::path_expand(gf),
    "-annot",
    fs::path_expand(af),
    "--load_zval",
    "-qtl",
    qtl_file
  )
  res <- processx::run(torus_path,args = res_args,echo_cmd = TRUE,echo = TRUE)
  fdr_res <- readr::read_tsv(qtl_file,col_names = c("rej","region_id","fdr","decision"))
  fs::file_delete(qtl_file)
  return(fdr_res)
}




fs_torus <- function(gwas_df,full_anno_df,steps=1L,p_cutoff=1,torus_p=character(0)){
  
  params <- unique(full_anno_df$feature)
  p <- nrow(gwas_df)
  gf <- write_gwas(gwas_df = gwas_df)
  torus_f <- partial(run_torus_cmd,gf=gf)
  
  combo_fun <- function(params){
    write_anno(anno_df = dplyr::filter(full_anno_df,feature %in% params),p = p)  
  }
  et_fun <- function(x){
    nc <- names(coef.torus(x))
    return(nc[nc!="Intercept"])
  }
  all_ret <- forward_select_fun(f = torus_f,params=params,combo_fun = combo_fun,extract_terms = et_fun,steps = steps)
  
  final_terms <- unnest(all_ret) %>% 
    filter(term!="Intercept",p<p_cutoff) %>% pull(term)
  
  taf <- combo_fun(final_terms)
  final_ret <- run_torus_cmd(gf=gf,af = taf,torus_p = torus_p)
  return(final_ret)
}
  

forward_op_torus_cmd <- function(gf,anno_df,f_feat_df,term_df,i,p,prior=NA_integer_,verbose=F){
  fn <- slice(term_df,i)
  miss_df <- anti_join(fn,unnest(f_feat_df))
  if(nrow(miss_df)==0){
    return(NULL)
  }
  tdf <- filter(anno_df,feature %in% c(term_list[i],f_feat))
  af <- write_anno(anno_df)
    tret_f <- daprcpp::run_torus_cmd(gf,af,torus_p=prior)
  fs::file_delete(af)
  return(lik_fun(tret_f))
}



