#include "dap_classdef.hpp"
#include <fstream>
#include <sstream>

namespace torus {


void show_banner(){

  fprintf(stderr, "\n\nTORUS: QTL Discovery Utilizing Genomic Annotations\n\n");
  fprintf(stderr, "Usage: torus -est|qtl|dump_prior  -d input_data.gz [-smap snp_map.gz] [-gmap gene_map.gz] [-annot annotation_file.gz] [--load_bf | --load_zval]\n\n");
  


}



controller cont_wrapper(const std::string gwas_file, const std::string annotation_file){
  int force_logistic = 0;
  int prob_annot = 0;

  int data_format = 1;

  double EM_thresh = 0.05;
  double dist_bin_size = -1;
  double l1_lambda = 0;
  double l2_lambda = 0;
  int find_egene = 0;
  int est = 0;
  double alpha = 0.05;

  double init_pi1 = 1e-3;

  controller con;
  con.EM_thresh = EM_thresh;

  if(dist_bin_size > 0){
    con.dist_bin_size = dist_bin_size;
  }

  if(force_logistic){
    con.force_logistic = 1;
  }

  if(prob_annot){
    con.single_fuzzy_annot = 1;
  }

  if(l1_lambda!=0){
    con.l1_lambda = l1_lambda;
    con.force_logistic = 1;
  }

  if(l2_lambda!=0){
    con.l2_lambda = l2_lambda;
    con.force_logistic =1;
  }
  int print_avg = 0;


  con.init_pi1 = init_pi1;
  con.print_avg = print_avg;

  con.load_data_zscore(gwas_file.c_str());


  con.load_annotation(annotation_file.c_str());
  fprintf(stderr,"Initializing ... \n");

  if(est==0 && find_egene==0){
    est = 1;
  }
  return con;
}




}
