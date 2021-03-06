#' Title
#'
#' @param anno_df 
#' @param p 
#' @param af 
#'
#' @return
#' @export
#'
#' @examples
write_anno <- function(anno_df=tibble(SNP=integer(),feature=character()),p=max(max(anno_df$SNP),1L),af=tempfile(fileext=".txt.gz")){
  if(is.null(anno_df)){
    return(write_anno(af=af))
  }
  spread_anno_l <- make_matrix(p = p,anno_df = anno_df)
  spread_anno_df <- tibble::as_tibble(magrittr::set_colnames(spread_anno_l$annomat,paste0(spread_anno_l$names,"_d"))) %>%
    dplyr::mutate(SNP=1:dplyr::n()) %>% 
    dplyr::select(SNP,dplyr::everything()) %>%
    dplyr::filter_at(.vars = dplyr::vars(-SNP),dplyr::any_vars(. != 0))
  readr::write_tsv(spread_anno_df,path=af)
  return(af)
  
}


                        
                           


#' Title
#'
#' @param gwas_df 
#' @param gf 
#'
#' @return
#' @export
#'
#' @examples
write_gwas <- function(gwas_df,gf=tempfile(fileext=".txt.gz")){
  
  dplyr::select(gwas_df,SNP,region_id,`z-stat`) %>% write_tsv(path=gf)
  return(gf)
}






#' Title
#'
#' @param gf 
#' @param af 
#' @param torus_p 
#'
#' @return
#' @export
#'
#' @examples
run_torus_cmd <- function(gf,af,torus_p=character(0),l1=NA_real_,l2=NA_real_,torus_path=system("which torus")){


  stopifnot(file.exists(gf),
            file.exists(af))
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
      lik_file
    )
  }
  if(!is.na(l1)){
    res_args <- c(res_args,"-l1_lambda",l1)
  }
  if(!is.na(l2)){
    res_args <- c(res_args,"-l2_lambda",l2)
  }
  res <- processx::run(torus_path,args = res_args,echo_cmd = TRUE,echo = TRUE)
  df <- read.table(file = textConnection(res$stdout),skip=1,header=F,stringsAsFactors = F)
  colnames(df) <- c("term", "estimate", "low", "high")
  
  df <- dplyr::mutate(df,term=stringr::str_replace(term,pattern = "\\.[0-9]+$",replacement = ""),
                      sd=(low-estimate)/(-1.96),z=estimate/sd,p=pnorm(abs(z),lower.tail = FALSE))
  lik <- scan(lik_file,what=numeric())
  file.remove(lik_file)
  df <- tidyr::nest(df) %>% dplyr::mutate(lik=lik)
  if(length(torus_p)>0){
      stopifnot(all(fs::file_exists(p_f)))
      prior_l <- purrr::map(torus_p,function(x){
          fp <- as.character(fs::path(torus_d,x,ext="prior"))
          suppressMessages(
              ret <- vroom::vroom(file = fp,delim = "  ",trim_ws = T,col_names = c("SNP","prior"),col_types = cols("SNP"="i","prior"="d")) %>% dplyr::mutate(region_id=x)
          )
          return(ret)
      })
      fs::file_delete(p_f)
      names(prior_l) <- torus_p
      ret <- list(df=df,priors=prior_l)
  }else{
      ret <- list(df=df)
  }
  return(ret)
}




#' Title
#'
#' @param gf 
#' @param af 
#'
#' @return
#' @export
#'
#' @examples
torus_fdr <- function(gf,af,torus_path=system("which torus")){
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


#' Title
#'
#' @param fit 
#'
#' @return
#' @export
#'
#' @examples
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
forward_select_fun <- function(f,params,combo_fun,extract_terms,steps=1L,parallel=parallel,init_params=character()){
  

  best_terms <- init_params
  rest_terms <- params[!params %in% best_terms]

  term_selection <- purrr::map(rest_terms,~c(.x,best_terms))
#  all_results <- list()
  
  db_fun <- function(ts){
    init_d <- combo_fun(ts)
    f(init_d)
  }
  
  for(i in seq_len(steps)){
    if(parallel){
      all_fit <- furrr::future_map(term_selection,db_fun)
    }else{
      all_fit <- purrr::map(term_selection,db_fun)
    }
    #    all_results[[i]] <- all_fit
    lik_vec <- map_dbl(all_fit,~.x$df$lik)
    if(length(lik_vec)==0){
      break
    }
    best_fit <- all_fit[[which.max(lik_vec)]]
    best_terms <- extract_terms(best_fit)
    stopifnot(all(best_terms %in% params))
    rest_terms <- rest_terms[! rest_terms %in% best_terms]
    term_selection <- purrr::map(rest_terms,~c(.x,best_terms))
  }
  # if(ret_all){
  #   return (all_results)
  # }else{
    return(best_fit) 
#  }
}


#' Title
#'
#' @param gwas_df 
#' @param full_anno_df 
#' @param steps 
#' @param p_cutoff 
#' @param torus_p 
#'
#' @return
#' @export
#'
#' @examples
fs_torus <- function(gf,p,full_anno_df,steps=1L,p_cutoff=1,torus_p=character(0),parallel=FALSE,init_terms=character()){
  
  params <- unique(full_anno_df$feature)
  torus_f <- purrr::partial(run_torus_cmd,gf=gf)
  stopifnot(all(init_terms %in% full_anno_df$feature))
  
  combo_fun <- function(params){
    write_anno(anno_df = dplyr::filter(full_anno_df,feature %in% params),p = p)  
  }
  et_fun <- function(x){
    nc <- names(coef.torus(x))
    return(nc[nc!="Intercept"])
  }
  all_ret <- forward_select_fun(f = torus_f,params=params,combo_fun = combo_fun,extract_terms = et_fun,steps = steps,parallel=parallel,init_params = init_terms)
  
  final_terms <- unnest(all_ret$df) %>% 
    filter(term!="Intercept",p<p_cutoff) %>% pull(term)
  
  taf <- combo_fun(final_terms)
  final_ret <- run_torus_cmd(gf=gf,af = taf,torus_p = torus_p)

  return(final_ret)
}
  

#' Title
#'
#' @param gf 
#' @param anno_df 
#' @param f_feat_df 
#' @param term_df 
#' @param i 
#' @param p 
#' @param prior 
#' @param verbose 
#'
#' @return
#' @export
#'
#' @examples
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



