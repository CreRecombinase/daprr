#pragma once
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <iostream>
#include <fstream>
#include <sstream>
#include "gsl/span"
#include <RcppEigen.h>
 #include <stdio.h>

 #ifdef NDEBUG
 # define assertr(EX)
 #else
 # define assertr(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))
 #endif





inline  void __assert(const char *msg, const char *file, const int line) {
  char buffer[100];
  snprintf(buffer, 100, "Assert Failure: %s at %s line #%d", msg, file, line);
  throw Rcpp::exception(buffer);
}

template<typename T,typename IT=Rcpp::IntegerVector::iterator>
class ProxyTrip {

  IT rvec;
  IT cvec;
  using trip_t =   Eigen::Triplet<T,typename Eigen::SparseMatrix<T>::StorageIndex>;
  trip_t ret;

public:
  ProxyTrip(IT rvec_, IT cvec_):rvec(rvec_),cvec(cvec_){
  }
  trip_t* operator->(){
    ret=trip_t((*rvec),(*cvec),1);

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

  bool operator ==(const ProxyTrip<T,IT> other)const {
    return( (rvec==other.rvec)&&(cvec==other.cvec));
  }
  bool operator !=(const ProxyTrip<T,IT> other)const {
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



template<typename T,typename IT=Rcpp::IntegerVector::iterator>
class Zip {
  IT rvec;
  IT cvec;
  T  dvec;
  using ret_tt = typename T::value_type;
  using trip_t =   Eigen::Triplet<ret_tt,typename Eigen::SparseMatrix<ret_tt>::StorageIndex>;
  trip_t ret;

public:
  Zip(IT rvec_, IT cvec_,T dvec_):rvec(rvec_),cvec(cvec_),dvec(dvec_){
  }
  trip_t* operator->(){
    ret=trip_t((*rvec),(*cvec),(*dvec));

    return(&ret);
  }
  Zip& operator++(){     // prefix
    rvec++;
    cvec++;
    dvec++;
    return(*this);
  }

  Zip& operator--(){  // prefix
    rvec--;
    cvec--;
    dvec--;
    return *this;

  }

  Zip operator++(int){ // postfix

    Zip temp(rvec,cvec,dvec);
    // Use prefix operator to increment this digit
    ++(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }

  bool operator ==(const Zip<T,IT> other)const {
    return( (rvec==other.rvec)&&(cvec==other.cvec)&&(other.dvec==dvec));
  }
  bool operator !=(const Zip<T,IT> other)const {
    return !((rvec==other.rvec)&&(cvec==other.cvec)&&(dvec==other.dvec)) ;
  }

  Zip operator--(int) { // postfix

    Zip temp(rvec, cvec, dvec);
    // Use prefix operator to increment this digit
    --(*this); // apply operator

    // return temporary result
    return temp; // return saved state
  }
};


class Writer{


  virtual void writeData(std::string field,double data) = 0;
  virtual void writeData(std::string field,int data) = 0;
  virtual void writeData(std::string field,std::string data) = 0;

  virtual double readData_d(std::string field,double data) = 0;
  virtual void readData_i(std::string field,int data) = 0;
  virtual void readData_s(std::string field,std::string data) = 0;

};



class FileWriter{


  virtual void writeData(std::string field,double data) = 0;
  virtual void writeData(std::string field,int data) = 0;
  virtual void writeData(std::string field,std::string data) = 0;

  virtual double readData_d(std::string field,double data) = 0;
  virtual void readData_i(std::string field,int data) = 0;
  virtual void readData_s(std::string field,std::string data) = 0;

};






class AnnotationParser {

public:
  using Xd_pt = ProxyTrip<int, gsl::span<int>::iterator>;
  using Xc_pt = Zip<gsl::span<double>::iterator, gsl::span<int>::iterator>;
  virtual Xd_pt Xd_begin() = 0;
  virtual Xd_pt Xd_end() = 0;
  virtual Xc_pt Xc_begin() = 0;
  virtual Xc_pt Xc_end() = 0;
  virtual std::vector<std::string> name_vec_d() const = 0;
  virtual std::vector<std::string> name_vec_c() const = 0;
  virtual int num_discrete() const = 0;
  virtual int num_continuous() const = 0;
};

class	DataAnnotationParser : public AnnotationParser{
  const size_t p;
  int kc; // number of continuous cov`ariate
  int kd; // number of discrete covariate
  std::vector<std::string> dvar_name_vec;
  std::vector<std::string> cvar_name_vec;
  gsl::span<int> row_d;
  gsl::span<int> col_d;
  gsl::span<int> row_c;
  gsl::span<int> col_c;
  gsl::span<double> data_c;
public:
  DataAnnotationParser(const size_t p_):
    p(p_){
  }
  void load_d(const gsl::span<int> row,const gsl::span<int> col,std::vector<std::string> names){
    dvar_name_vec = names;
    kd=names.size();
    const size_t anno_p = row.size();
    if (anno_p != col.size()) {
      throw std::invalid_argument("row.size()!=col.size() in load_d");
    }
    row_d = row;
    col_d = col;
  }
  void load_c(gsl::span<int> row,gsl::span<int> col,gsl::span<double> data,std::vector<std::string> names){
    cvar_name_vec = names;
    kc=names.size();
    const size_t anno_p = row.size();
    if(anno_p!=col.size() || anno_p != data.size()){
      throw std::invalid_argument( "row.size()!=col.size()!=data.size() in load_c" );
    }
    row_c=row;
    col_c=col;
    data_c=data;
  }
  std::vector<std::string> name_vec_d() const{
    return dvar_name_vec;
  }
  std::vector<std::string> name_vec_c() const{
    return cvar_name_vec;
  }
  Xd_pt Xd_begin() {
    Xd_pt ret(row_d.begin(), col_d.begin());
    return (ret);
  }
  Xd_pt Xd_end() {
    Xd_pt ret(row_d.end(), col_d.end());
    return (ret);
  }
  Xc_pt Xc_begin(){
    Xc_pt ret(row_c.begin(),col_c.begin(),data_c.begin());
    return(ret);
  }
  Xc_pt Xc_end(){
    Xc_pt ret(row_c.end(),col_c.end(),data_c.end());
    return(ret);
  }
  int	num_discrete() const{
    return kd;
  }
  int num_continuous() const { return kc; }
};

class FileAnnotationParser : public AnnotationParser{
  const size_t p;
  int kc; // number of continuous cov`ariate
  int kd; // number of discrete covariate
  std::vector<int> row_d;
  std::vector<int> col_d;
  std::vector<int> row_c;
  std::vector<int> col_c;
  std::vector<double> data_c;
  std::vector<std::string> cvar_name_vec;
  std::vector<std::string> dvar_name_vec;
  gsl::span<int> rows_d;
  gsl::span<int> cols_d;
  std::string snp;
public:
  FileAnnotationParser(const std::map<std::string,int> &snp_hash,const char* annot_file):
    p(snp_hash.size()),
    kc(0),
    kd(0){
    std::ifstream afile(annot_file, std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream inf;
    std::string line;
    std::istringstream ins;
    int count_d;
    int count_c;
    enum class cat { disc, cont };
    std::map<int, cat> col2cat;
    std::map<int, int> col2cpos;
    std::map<int, int> col2dpos;

    inf.push(boost::iostreams::gzip_decompressor());
    inf.push(afile);
    int col_count = -1;
    while(std::getline(inf,line)){
      ins.clear();
	  ins.str(line);
	  std::string token;
	  while(ins>>token){
	    if(col_count==-1){
	      if(token == "SNP" || token == "snp"){
		col_count = 0;
		continue;
	      }else{
		break;
	      }
	    }
	    std::string cat = token.substr(token.size()-2, 2);
	  // continuous
	    if(cat == "_c" || cat =="_C"){
	      col2cat[col_count] = cat::cont;
	      col2cpos[col_count] = kc;
	      std::string name = token.substr(0,token.size()-2);
	      cvar_name_vec.push_back(name);
	      kc++;
	    } // discrete/categorical
	    else{
	      col2cat[col_count] = cat::disc;
	      col2dpos[col_count] = kd;
	      std::string name = token.substr(0,token.size()-2);
	      dvar_name_vec.push_back(name);
	    kd++;
	  }
	  col_count++;
	}
	if(col_count!=-1)
	  break;
      }
      if(col_count==-1){
	throw std::length_error("col_count is -1");
      }
      while(getline(inf,line)){
	ins.clear();
	ins.str(line);
	snp.clear();
	ins>>snp;
	auto snp_h  = snp_hash.find(snp);
	if(snp_h != snp_hash.end()){
	  int col_ct=0;
	  double val;
	  while(ins>>val){
	    if(col2cat[col_ct]==cat::cont){
	      row_c.push_back(snp_h->second);
	      col_c.push_back(col2cpos[col_ct]);
	      data_c.push_back(val);
	    }else{
              if (val != 0) {
                row_d.push_back(snp_h->second);
                col_d.push_back(col2dpos[col_ct]);
              }
            }
            col_ct++;
          }
        }
      }
      rows_d=gsl::span<int>(row_d);
      cols_d=gsl::span<int>(col_d);

  }
  Xd_pt Xd_begin(){
    Xd_pt r(rows_d.begin(),cols_d.begin());
    return(r);
  }
  Xc_pt Xc_begin(){
    gsl::span<int> rows(row_c);
    gsl::span<int> cols(col_c);
    gsl::span<double> datas(data_c);
    Xc_pt r(rows.begin(),cols.begin(),datas.begin());
    return(r);
  }
  Xd_pt Xd_end(){
    Xd_pt r(rows_d.end(),cols_d.end());
    return(r);
  }
  Xc_pt Xc_end(){
    gsl::span<int> rows(row_c);
    gsl::span<int> cols(col_c);
    gsl::span<double> datas(data_c);
    Xc_pt r(rows.end(),cols.end(),datas.end());
    return(r);
  }
  int	num_discrete() const{
    return kd;
  }
  int num_continuous() const { return kc; }
    std::vector<std::string> name_vec_d() const{
    return dvar_name_vec;
  }
  std::vector<std::string> name_vec_c() const{
    return cvar_name_vec;
  }

};

class ZscoreParser {

public:
  using linetype_s = std::tuple<std::string, std::string, double>;
  using linetype_i = std::tuple<int, int, double>;
  virtual linetype_s* getline_s() = 0;
  virtual linetype_i* getline_i() = 0;
};


class DataZscoreParser : public ZscoreParser {
  gsl::span<std::string> SNP;
  gsl::span<std::string> locus;
  gsl::span<double> z;
  linetype_i line_tupi;
  linetype_s line_tups;
  const size_t p;
  int i;
  std::unordered_set<std::string> loc_map;

public:
  DataZscoreParser(gsl::span<std::string> SNP_, gsl::span<std::string> locus_,
                   gsl::span<double> z_)
    : SNP(SNP_), locus(locus_), z(z_),p(SNP.size()),i(0) {}

  linetype_s *getline_s() {
    if (i < p) {
      loc_map.insert(locus[i]);
      i++;
      line_tups = std::make_tuple(SNP[i], locus[i], z[i]);
      return (&line_tups);
    } else {
      return (nullptr);
    }
  }
  linetype_i *getline_i() {
    if (i < p) {
      loc_map.insert(locus[i]);
      i++;
      line_tupi = {i, loc_map.size(), z[i]};
      return (&line_tupi);
    } else {
      return (nullptr);
    }
  }
};

class FileZscoreParser: public ZscoreParser {

  const char *filename;
  std::ifstream dfile;
  boost::iostreams::filtering_istream insd;
  std::string line;
  std::istringstream ins;

  linetype_i line_tupi;
  linetype_s line_tups;
  std::string &snp_id;
  std::string &loc_id;
  std::unordered_set<std::string> loc_map;
  double &z_val;
  int i;
public:
  FileZscoreParser(const char *filename_)
      : filename(filename_),
        dfile(filename, std::ios_base::in | std::ios_base::binary), line_tupi(),
	line_tups(),
        snp_id(std::get<0>(line_tups)), loc_id(std::get<1>(line_tups)),
        z_val(std::get<2>(line_tups)),i(0) {
    insd.push(boost::iostreams::gzip_decompressor());
    insd.push(dfile);
    ins.clear();
    ins.str(line);
    std::getline(insd, line);
    if (line != "SNP locus z-val") {
      throw std::length_error(line);
    }
  }
  linetype_s *getline_s() {
    if (std::getline(insd, line)) {
      ins.clear();
      ins.str(line);
      ins >> std::get<0>(line_tups) >> std::get<1>(line_tups) >>
          std::get<2>(line_tups);
      loc_map.insert(std::get<1>(line_tups));
      i++;
      return (&line_tups);
    } else {
      return(nullptr);
    }
  }
  linetype_i *getline_i() {

    if (std::getline(insd, line)) {
      ins.clear();
      ins.str(line);
      ins >> std::get<0>(line_tups) >> std::get<1>(line_tups) >>
          std::get<2>(line_tups);
      loc_map.insert(std::get<1>(line_tups));
      line_tupi = {i, loc_map.size(), z_val};
      i++;
      return (&line_tupi);
    } else {
      return (nullptr);
    }
  }
};
