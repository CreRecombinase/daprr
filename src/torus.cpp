#define EIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS
#include "classdef.h"
#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <gsl/gsl_multifit.h>
#include <cmath>
//[[Rcpp::depends(RcppEigen)]]
using	loc_it = Rcpp::IntegerVector::iterator;












// [[Rcpp::export]]
Rcpp::List torus(Rcpp::IntegerVector locus_id, Rcpp::NumericVector z_hat,Eigen::Map<Eigen::SparseMatrix<double> > anno_mat,Rcpp::StringVector names,const bool prior=false){


    double EM_thresh = 0.05;
    double init_pi1 = 1e-3;
    int print_avg = 0;
    auto split = make_splitter(locus_id.begin(),locus_id.end());
    Result_obj res(locus_id.size(),anno_mat.ncol()+1);
    controller con(split,res,anno_mat);
    con.EM_thresh = EM_thresh;
    con.init_pi1 = init_pi1;
    con.print_avg = print_avg;
    con.load_data_R(z_hat);
    con.load_annotations_R(Rcpp::as<std::vector<std::string>>(names));
    auto result = con.estimate();
    if(!prior){
        return result;
    }else{
        using namespace Rcpp;
        return (Rcpp::List::create(_["prior"]=res.r_prior,_["est"]=result));           
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
