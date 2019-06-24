#include <RcppGSL.h>
#include "classdef.hpp"
#include <gsl/gsl_multifit.h>
#include <cmath>
#include <variant>

using	loc_it = Rcpp::IntegerVector::iterator;
using namespace Rcpp;


std::pair<RcppGSL::matrix<int>,Rcpp::StringVector> gen_anno_mat(const size_t p,Rcpp::DataFrame anno_df){
        
    Rcpp::IntegerVector anno_row_id = anno_df["SNP"];
    SEXP tr = anno_df["feature"];
    
    auto t = TYPEOF(tr);
    int k=0;
    const size_t nr = anno_row_id.size();

    Rcpp::IntegerVector anno_col_id;
    Rcpp::StringVector names;
    if(t ==INTSXP){
        anno_col_id = anno_df["feature"];
        Rcpp::Nullable<Rcpp::StringVector> n_names =  anno_col_id.attr("levels");
        if(n_names.isNull()){
            Rcpp::stop("Features must be named");
        }
        names=n_names;
        k=names.size();
    } else{
        if(t!=STRSXP){
            Rcpp::stop("column `feature` in anno_df  must be int(factor) or character");
        }
        std::unordered_map<std::string,int> u_names;
        Rcpp::StringVector feat_v = anno_df["feature"];
        anno_col_id = Rcpp::IntegerVector(nr);

        std::string buffer ;
        bool chk;
        int i=0;
        for(auto fv: feat_v ){
            buffer = fv;
            auto u_name_i = u_names.find(buffer);
            if(u_name_i==u_names.end()){
                auto mp = u_names.insert({buffer,++k});
                u_name_i = mp.first;
            }
            anno_col_id(i++)=u_name_i->second;
        }
        u_names.size();
        std::vector<std::string> tnames(k);
        for(auto [mf,ti] :u_names ){
            tnames[ti-1]=mf;
        }
        names = Rcpp::wrap(tnames);
    }
    
    
    RcppGSL::matrix<int> anno_mat(p,std::max(k,1));
    // auto dp = anno_mat.data;
    // for(int i=0; i<p*k; i++){}
    //   
    std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));
    
    for(int i=0; i<nr; i++){
        anno_mat(anno_row_id(i)-1,anno_col_id(i)-1)=1;
    }
    if(nr==0){
      names=anno_df.attr("feature");
    }
    return(std::make_pair(anno_mat,names));
  
};


//[[Rcpp::export]]
Rcpp::List make_matrix(const size_t p,Rcpp::DataFrame anno_df){
  auto [mat,names] = gen_anno_mat(p,anno_df);
  using namespace Rcpp;
  return(List::create(_["annomat"]=mat,_["names"]=names));
}


bool verbose=true;

// [[Rcpp::export]]
Rcpp::List torus_df(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,Rcpp::DataFrame anno_df,const bool prior=false,const bool do_verbose=false){

        
    const size_t p = locus_id.size();


    if(!do_verbose){
        verbose=false;
    }else{
        verbose=true;
    }
    double EM_thresh = 0.05;
    double init_pi1 = 1e-3;
    int print_avg = 0;
    auto [anno_mat,names] = gen_anno_mat(p,anno_df);
    auto split = make_splitter(locus_id.begin(),locus_id.end());
    Result_obj res(locus_id.size(),anno_mat.ncol()+1);
    controller con(split,res);
    con.EM_thresh = EM_thresh;
    con.init_pi1 = init_pi1;
    con.print_avg = print_avg;
    con.load_data_R(z_hat);
    con.load_annotations_R(anno_mat,Rcpp::as<std::vector<std::string>>(names));
    try{
      auto result = con.estimate();
      using namespace Rcpp;
      if(!prior){
        return Rcpp::List::create(_["est"]=result);
      }else{

        return (Rcpp::List::create(_["prior"]=res.prior,_["est"]=result));           
      }
    }catch(std::exception e){
      Rcpp::Rcerr<<e.what()<<std::endl;
      Rcpp::stop("Caught exception from torus!");
    }catch(int e){
      Rcpp::Rcerr<<e<<std::endl;
      Rcpp::stop("Caught exception from GSL!");

    }

}


// [[Rcpp::export]]
Rcpp::List torus(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,RcppGSL::matrix<int> anno_mat,Rcpp::StringVector names,const bool prior=false,const bool do_verbose=false){


    if(!do_verbose){
        verbose=false;
    }else{
        verbose=true;
    }
    double EM_thresh = 0.05;
    double init_pi1 = 1e-3;
    int print_avg = 0;
    auto split = make_splitter(locus_id.begin(),locus_id.end());
    Result_obj res(locus_id.size(),anno_mat.ncol()+1);
    controller con(split,res);
    con.EM_thresh = EM_thresh;
    con.init_pi1 = init_pi1;
    con.print_avg = print_avg;
    con.load_data_R(z_hat);
    con.load_annotations_R(anno_mat,Rcpp::as<std::vector<std::string>>(names));
    auto result = con.estimate();
    using namespace Rcpp;
    if(!prior){
        return Rcpp::List::create(_["est"]=result);
    }else{

        return (Rcpp::List::create(_["prior"]=res.prior,_["est"]=result));           
    }
}

// Rcpp::NumericVector logistic_cat_pred(RcppGSL::vector<double> beta,  // Vector of parameters length = 1 + Sum_k(C_k - 1)
//                          RcppGSL::matrix<int> X,  //Matrix Nobs x K
//                          RcppGSL::vector<int> nlev,
//                          RcppGSL::vector<double> yhat
//                                         )
// {

//   for(int i = 0; i < X->size1; ++i) {
//     double Xbetai=beta->data[0];
//     int iParm=1;
//     for(int k = 0; k < X->size2; ++k) {
//       if(gsl_matrix_int_get(X,i,k)>0)
        
//         Xbetai+=beta->data[gsl_matrix_int_get(X,i,k)-1+iParm];
//       iParm+=nlev->data[k]-1;
//     }
//     yhat->data[i]=1/(1 + std::exp(-Xbetai));
//   }
//   return(Rcpp::wrap(yhat));
// }
