#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_odeiv2.h>
#include <dlfcn.h>
#include "h5block.h"

#define DEFAULT 1
#define FREE_ON_SUCCESS 1
#define KEEP_ON_SUCCESS 2


char* model_function(char *model_name, char *suffix){
  assert(model_name);
  char *f=malloc(sizeof(char)*(strlen(model_name)+strlen(suffix)+1));
  assert(f);
  *mempcpy(mempcpy(f,model_name,strlen(model_name)),suffix,strlen(suffix))='\0';
  printf(fprintf,"[%s] function: «%s»\n",__func__,f);
  return f;
}

void *load_or_exit(void *lib, char *name, int opt){
  void *symbol=dlsym(lib,name);
  if (symbol) {
    fprintf(stderr,"[%s] «%s» loaded successfully.\n",__func__,name);
    if (opt==FREE_ON_SUCCESS){
      free(name);
    }
  }else{
    fprintf(stderr,"[%s] %s\n",__func__,dlerror());
    exit(-2);
  } 
  return symbol;
}

gsl_odeiv2_system* load_system(char *model_name, size_t n, double *p){
  char *so=model_function(model_name,".so");
  void *lib=dlopen(so);
  void *f,*dfdy,*dfdp;
  int i;
  char *symbol_name; // symbol name in .so
  if (lib){
    symbol_name=model_function(model_name,"_vf");
    f=load_or_exit(so,symbol_name,FREE_ON_EXIT);
    symbol_name=model_function(model_name,"_jac");
    dfdy=load_or_exit(so,symbol_name,FREE_ON_EXIT);
  } else {
    fprintf(stderr,"[%s] library «%s» could not be loaded: %s\n",so,dlerror());
    exit(-1)
  }
  gsl_odeiv2_system *sys={f,dfdy,n,p};
  return sys;
}

int option_is(char *s, char *l, char *value){
  int match=FALSE;
  assert(value);
  if (s) match=(strcmp(value,s)==0);
  else if (l) match|=(strcmp(value,l)==0);
  return match;
}

int main(int argc, char *argv[]){
  int i=0;
  char *model_name=NULL, *h5file=NULL;
  double abs_tol=1e-6,rel_tol=1e-5,h=1e-3;
  while (i<argc){
    if (option_is("-m","--model",argv[i])){
      i++;
      model_name=malloc(sizeof(char)*(strlen(argv[i])+1));
      strcpy(model_so,argv[i]);
    } else if (option_is("-d","--data",argv[i])){
      i++;
      h5file=malloc(sizeof(char)*(strlen(argv[i])+1));
      strcpy(h5file,argv[i]);
    } else if (option_is("-a","--abs-tol",argv[i])){
      i++;
      abs_tol=strtod(argv[i],NULL);
    } else if (option_is("-r","--rel-tol",argv[i])){
      i++;
      rel_tol=strtod(argv[i],NULL);
    } else if (option_is("-h","--step-size",argv[i])){
      i++;
      h=strtod(argv[i],NULL);
    }else {
      fprintf(stderr,"[%s] unknown option «%s»\n",__func__,argv[i]);
      exit(1);
    }
    i++;
  }
  assert(h5file);
  assert(model_name);
  // load data related to this model
  h5block_t* h5=h5block_alloc(2);
  hid_t data=H5Fopen(h5file,H5F_ACC_RDONLY,H5P_DEFAULT);
  hid_t prior=H5Gopen2(data,"prior",H5P_DEFAULT);

  gsl_vector *p=h5_to_gsl(prior,"mu",NULL);
  size_t i,n=p->size;
  fprintf(stderr,"[%s] prior: ",__func__);
  for (i=0;i<n;i++) fprintf(stderr,"%g ",gsl_vector_get(p,i));
  fprintf(stderr,"\n");

  size_t d=???; // get from initial conditions
  
  gsl_odeiv2_system* sys=load_system(model_name, d, p->data);
  gsl_odeiv2_step *s;
  gsl_odeiv2_driver *driver=gsl_odeiv2_driver_alloc_y_new(sys,gsl_odeiv2_step_msbdf,h,abs_tol,rel_tol)
  return EXIT_SUCCESS;
}
