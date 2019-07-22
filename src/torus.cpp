#include <RcppGSL.h>
#include <RcppEigen.h>
#include <algorithm>
#include <cmath>

using namespace Eigen;


using	loc_it = Rcpp::IntegerVector::iterator;
using namespace Rcpp;
bool verbose=true;

template<typename T>
std::string sexp_f(const T tr){

  if  (tr == 0){ return "NILSXP";}    /* nil = NULL */
  if  (tr == 1){ return "SYMSXP";}    /* symbols */
  if  (tr == 2){ return "LISTSXP";}    /* lists of dotted pairs */
  if  (tr == 3){ return "CLOSXP";}    /* closures */
  if  (tr == 4){ return "ENVSXP";}    /* environments */
  if  (tr == 5){ return "PROMSXP";}    /* promises: [un]evaluated closure arguments */
  if  (tr == 6){ return "LANGSXP";}    /* language constructs (special lists) */
  if  (tr == 7){ return "SPECIALSXP";}    /* special forms */
  if  (tr == 8){ return "BUILTINSXP";}    /* builtin non-special forms */
  if  (tr == 9){ return "CHARSXP";}    /* "scalar" string type (internal only)*/
  if  (tr == 10){ return "LGLSXP";}   /* logical vectors */
  if  (tr == 13){ return "INTSXP";}   /* integer vectors */
  if  (tr == 14){ return "REALSXP";}   /* real variables */
  if  (tr == 15){ return "CPLXSXP";}   /* complex variables */
  if  (tr == 16){ return "STRSXP";}   /* string vectors */
  if  (tr == 17){ return "DOTSXP";}   /* dot-dot-dot object */
  if  (tr == 18){ return "ANYSXP";}   /* make "any" args work */
  if  (tr == 19){ return "VECSXP";}   /* generic vectors */
  if  (tr == 20){ return "EXPRSXP";}   /* expressions vectors */
  if  (tr == 21){ return "BCODESXP";}   /* byte code */
  if  (tr == 22){ return "EXTPTRSXP";}   /* external pointer */
  if  (tr == 23){ return "WEAKREFSXP";}   /* weak reference */
  if  (tr == 24){ return "RAWSXP";}   /* raw bytes */
  if  (tr == 25){ return "S4SXP";}
  if  (tr == 30){ return "NEWSXP";}   /* fresh node creaed in new page */
  if  (tr == 31){ return "FREESXP";}   /* node released by GC */
  if  (tr == 99){ return "FUNSXP";}
  return "Don't Know";
};



class factor_col{
public:
  Rcpp::IntegerVector vec;
  Rcpp::StringVector names;
  factor_col(SEXP tr,const bool need_names=false){
    auto t = TYPEOF(tr);
    int k=0;
    if(t ==INTSXP){
      vec = tr;
      Rcpp::Nullable<Rcpp::StringVector> n_names =  vec.attr("levels");
      if(need_names){
	if(n_names.isNull()){
	  Rcpp::stop("Features must be named");
	}else{
	  names=n_names;
	  k=names.size();
	}
      }
    } else{
      if(t!=STRSXP){
	Rcpp::Rcerr<<"t is of type"<<sexp_f(t)<<std::endl;
	Rcpp::stop("column  must be int(factor) or character");
      }
      std::unordered_map<std::string,int> u_names;
      Rcpp::StringVector feat_v = tr;
      vec = Rcpp::IntegerVector(feat_v.size());
      int i=0;
      std::string buffer;
      for(auto fv: feat_v ){
	buffer = fv;
	auto u_name_i = u_names.find(buffer);
	if(u_name_i==u_names.end()){
	  auto mp = u_names.insert({buffer,++k});
	  u_name_i = mp.first;
	}
	vec(i++)=u_name_i->second;
      }
      u_names.size();
      std::vector<std::string> tnames(k);
      
      for(auto mp :u_names ){
        auto mf = mp.first;
        auto ti = mp.second;
        tnames[ti-1]=mf;
      }
      names = Rcpp::wrap(tnames);
    }
  }
};

template<typename T>
class ProxyTrip {
  
  Rcpp::IntegerVector::iterator rvec;
  Rcpp::IntegerVector::iterator cvec;
  using trip_t =   Eigen::Triplet<T,typename Eigen::SparseMatrix<T>::StorageIndex>;
  trip_t ret;

public:
  ProxyTrip(Rcpp::IntegerVector::iterator rvec_, Rcpp::IntegerVector::iterator cvec_):rvec(rvec_),cvec(cvec_){
  }
  trip_t* operator->(){
    ret=trip_t((*rvec)-1,(*cvec)-1,1);

    return(&ret);
  }
  ProxyTrip& operator++(){     // prefix
    rvec++;
    cvec++;
    return(*this);
  }

  ProxyTrip& operator--(){  // prefix
    rvec--;
    cvec--;
    return *this;

  }

  ProxyTrip operator++(int){ // postfix

    ProxyTrip temp(rvec,cvec);
    // Use prefix operator to increment this digit
    ++(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }

  bool operator ==(const ProxyTrip<T> other)const {
    return( (rvec==other.rvec)&&(cvec==other.cvec));
  }
  bool operator !=(const ProxyTrip<T> other)const {
    return !((rvec==other.rvec)&&(cvec==other.cvec)) ;
  }

  ProxyTrip operator--(int){ // postfix

    ProxyTrip temp(rvec,cvec);
    // Use prefix operator to increment this digit
    --(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }
 
};


class SparseDF{
  factor_col anno_row_id;
  factor_col anno_col_id;
  const size_t p;
  const int k;
public:
  SparseDF(Rcpp::DataFrame anno_df,const size_t p_,const std::string row_name="SNP",const std::string col_name = "feature"):
  anno_row_id(anno_df[row_name],false),
  anno_col_id(anno_df[col_name],true),
  p(p_),
  k(anno_col_id.names.size()){

  }

  RcppGSL::matrix<int> getMat(){
    RcppGSL::matrix<int> anno_mat(p,std::max(k,1));

    //std::fill(anno_mat.begin(),anno_mat.end(),0);
        std::memset(anno_mat->data,0,anno_mat->size1*anno_mat->size2*sizeof(int));
    const size_t nr = anno_row_id.vec.size();
    for(int i=0; i<nr; i++){
      anno_mat(anno_row_id.vec(i)-1,anno_col_id.vec(i)-1)=1;
    }
    return(anno_mat);
  }
  Rcpp::StringVector names(){
    return anno_col_id.names;
  }
};




//[[Rcpp::export]]
Rcpp::List make_matrix(const size_t p,Rcpp::DataFrame anno_df){


  SparseDF spdf( anno_df, p,"SNP","feature");

  using namespace Rcpp;
  return(List::create(_["annomat"]=spdf.getMat(),_["names"]=spdf.names()));
}


