#include "classdef.hpp"

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

        int i=0;
	std::string buffer;
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

    std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));

    for(int i=0; i<nr; i++){
        anno_mat(anno_row_id(i)-1,anno_col_id(i)-1)=1;
    }
    if(nr==0){
      names=anno_df.attr("feature");
    }
    return(std::make_pair(anno_mat,names));
};
