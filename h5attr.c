#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "ndarray.h"

/* A command line option has the pattern `"-o"` or `"--option"`. This
   function checks whether a given string `value` is the name of a
   short option `s` or long option `l`.*/
int /* returns `1` if the given string is the specified option, `0` otherwise */
option_is(const char *s, /* short option name, e.g.: "-t", can be NULL */
 const char *l, /* long option name, e.g.: "--time", can also be NULL */
 const char *value) /* a string for comaprison with s and l*/
{
  int match=0;
  assert(value);
  if (s) match=(strcmp(value,s)==0);
  if (l) match|=(strcmp(value,l)==0);
  return match;
}

int file_extension_is(const char *name, char *ext){
  char *dot=strrchr(name, '.');
  return (strcmp(dot,ext)==0);
}

int main(int argc, char **argv){
  char *dataset_name;
  char *h5file=NULL;
  char *attr_name=NULL;
  char *infile=NULL;
  char *txt=NULL;
  char *value_str=NULL;
  ndarray *attr_value=NULL;
  int i;
  for (i=1;i<argc;i++){
    if (option_is("-d","--data-set",argv[i])){
      dataset_name=strdup(argv[++i]);
    } else if (option_is("-a","--attribute",argv[i])){
      attr_name=strdup(argv[++i]);
    } else if (option_is("--string-value","-s",argv[i])){
      value_str=strdup(argv[++i]);
    } else if (file_extension_is(argv[i],".h5")) {
      h5file=strdup(argv[i]);
    } else if (file_extension_is(argv[i],".txt")){
      txt=strdup(argv[i]);
    } else { //assume that it is a binary file
      infile=strdup(argv[i]);
    }
  }
  assert(h5file);
  assert(attr_name);
  assert(dataset_name);

  if (value_str){
    attr_value=ndarray_from_string(value_str);
  } else if (txt){
    attr_value=ndarray_from_text_file(txt);
  } else if (infile){
    attr_value=ndarray_from_binary_file(infile);
  } else {
    abort();
  }
  assert(attr_value && attr_value->rank==1 && attr_value->size[0]>0);
  hid_t h5f=H5Fopen(h5file,H5F_ACC_RDWR,H5P_DEFAULT);
  printf("[%s] writing attribute of dataset «%s»:\n",__func__, dataset_name);
  ndarray_print(attr_value,attr_name);
  ndarray_to_h5attr(attr_value,h5f,dataset_name,attr_name);
  return EXIT_SUCCESS;
}

