#include "classdef.h"
#include <fstream>
#include <sstream>


using namespace torus;

void show_banner(){

  fprintf(stderr, "\n\nTORUS: QTL Discovery Utilizing Genomic Annotations\n\n");
  fprintf(stderr, "Usage: torus -est|qtl|dump_prior  -d input_data.gz [-smap snp_map.gz] [-gmap gene_map.gz] [-annot annotation_file.gz] [--load_bf | --load_zval]\n\n");
  


}



int main(int argc, char **argv){
  
  // creating the grid
  
  //olist.push_back(0.1);
  //phlist.push_back(0.05);

  std::string data_file;
  std::string gmap_file;
  std::string smap_file;
  std::string annot_file;
  std::string prior_dir;
  std::string lik_file;
  std::string output_pip;

  int csize = -1;
  int gsize = -1;
  int nthread = 1;
  int print_avg = 0;
  
  int fastqtl_use_dtss = 1;






  std::string init_file;
  std::string qtl_file;

  int find_egene = 0;
  int est = 0;
  double alpha = 0.05;

  double init_pi1 = 1e-3;

  std::string ci_file;

  int force_logistic = 0;
  int prob_annot = 0;
  
  int data_format = 1;

  double EM_thresh = 0.05;
  double dist_bin_size = -1;
  double l1_lambda = 0;
  double l2_lambda = 0;
  
  for(int i=1;i<argc;i++){
    
    if(strcmp(argv[i], "-d")==0 || strcmp(argv[i], "-data")==0){
      data_file=argv[++i];
      continue;
    }

    if(strcmp(argv[i], "-lik")==0){
      lik_file=argv[++i];
      continue;
    }



    if(strcmp(argv[i], "-gmap")==0){
      gmap_file=argv[++i];
      continue;
    }

    
    if(strcmp(argv[i], "-smap")==0 ){
      smap_file=argv[++i];
      continue;
    }


    if(strcmp(argv[i], "-annot")==0 ){
      annot_file=argv[++i];
      continue;
    }
    
    if(strcmp(argv[i], "-t")==0 || strcmp(argv[i], "-thresh")==0){
      EM_thresh = atof(argv[++i]);
      continue;
    }
    
    if(strcmp(argv[i], "-l1_lambda")==0){
      l1_lambda = atof(argv[++i]);
      continue;
    }

    if(strcmp(argv[i], "-l2_lambda")==0){
      l2_lambda = atof(argv[++i]);
      continue;
    }
    

    if(strcmp(argv[i], "--single_fuzzy_annot")==0){
      prob_annot = 1;
      continue;
    }

    
    if(strcmp(argv[i], "-init_pi1")==0){
      init_pi1 = atof(argv[++i]);
      continue;
    }


    if(strcmp(argv[i], "--force_logistic")==0){
      force_logistic = 1;
      continue;
    }
    
    if(strcmp(argv[i], "--load_bf")==0 || strcmp(argv[i], "--bf")==0 ){
      data_format = 2;
      continue;
    }
    
   
    if(strcmp(argv[i], "--load_zval")==0 || strcmp(argv[i], "--zval")==0){
      data_format = 3;
      continue;
    }
       
    
    if(strcmp(argv[i], "--load_fastqtl")==0 || strcmp(argv[i], "--fastqtl")==0){
      data_format = 4;
      continue;
    }

    if(strcmp(argv[i], "--load_matrixeqtl")==0 || strcmp(argv[i], "--matrixeqtl")==0){
      data_format = 1;
      continue;
    }

    if(strcmp(argv[i], "--no_dtss") == 0) {
      fastqtl_use_dtss = 0;
      continue;
    }
    
    



    if(strcmp(argv[i], "-dist_bin_size") == 0){
      dist_bin_size = atof(argv[++i]);
      continue;
    }

    
    if(strcmp(argv[i], "-est")==0){
      est = 1;
      continue;
    }
    
    
    if(strcmp(argv[i], "-egene")==0 || strcmp(argv[i], "-qtl")==0 ){
      qtl_file=argv[++i];
      find_egene = 1;
      continue;
    }

    if(strcmp(argv[i], "-dump_prior")==0){
      prior_dir=argv[++i];
      continue;
    }

    if(strcmp(argv[i], "-dump_pip")==0){
      output_pip=argv[++i];
      continue;
    }


    if(strcmp(argv[i], "-h")==0 || strcmp(argv[i], "-help")==0 ){
      show_banner();
      continue;
    }


    if(strcmp(argv[i], "--print_avg")==0 ){
      print_avg = 1;
      continue;
    }


    if(strcmp(argv[i], "-alpha")==0){
      alpha = atof(argv[++i]);
      continue;
    }


    fprintf(stderr, "Error: undefined option %s\n", argv[i]);
    show_banner();
    exit(0);

  }    



  // checking mandatory arguments
if(strlen(data_file.c_str())==0){
    fprintf(stderr,"Error: data file unspecified\n");
    show_banner();
    exit(0);
  }

    
  // a global variable 
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
  
  
  con.init_pi1 = init_pi1;
  con.print_avg = print_avg;

  switch(data_format){
  case 1:
    con.load_data(data_file.c_str());
    break;
  case 2:
    con.load_data_BF(data_file.c_str());
    break;
  case 3:
    con.load_data_zscore(data_file.c_str());
    break;
  case 4:
    con.fastqtl_use_dtss = fastqtl_use_dtss;
    con.load_data_fastqtl(data_file.c_str());
    gmap_file[0]=smap_file[0] = 0;
    break;
  default:
    con.load_data(data_file.c_str());
    break;
  }
   

  con.load_map(gmap_file.c_str(), smap_file.c_str());  
  con.load_annotation(annot_file.c_str());
  fprintf(stderr,"Initializing ... \n");

  if(est==0 && find_egene==0){
    est = 1;
  }


  if(est)
    con.estimate();
  if(find_egene){
    con.find_eGene(qtl_file.c_str(),alpha);
  }
  if(prior_dir.size()){
    con.dump_prior(prior_dir.c_str());
  }
  if(output_pip.size()>0){
    fprintf(stderr,"#-dump_pip output_pip file: %s \n",output_pip.c_str());
con.dump_pip(output_pip.c_str());
  }

if(lik_file.size()>0){
    std::string filename = lik_file;
    std::ofstream ostrm(filename);
    ostrm << con.final_log10_lik<<std::endl;
  }

  return 0;
  }
