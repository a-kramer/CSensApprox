#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "h5block.h"
#include "ndarray.h"

/* A command line option has the pattern `"-o"` or `"--option"`. This
   function checks whether a given string `value` is the name of a
   short option `s` or long option `l`.*/
int /* returns `1` if the given string is the specified option, `0` otherwise */
option_is(const char *s, /* short option name, e.g.: "-t", can be NULL */
 const char *l, /* long option name, e.g.: "--time", can also be NULL */
 const char *value) /* a string for comaprison with s and l*/
{
  int match=FALSE;
  assert(value);
  if (s) match=(strcmp(value,s)==0);
  if (l) match|=(strcmp(value,l)==0);
  return match;
}

int file_extension_is(const char *name, char *ext){
  char *dot=strrchr(name, '.');
  return (strcmp(dot,ext)==0)
}

int main(int argc, char **argv){
  char *dataset_name;
  char *h5file=NULL;
  char *attr=NULL;
  char *infile=NULL;
  char *txt=NULL;
  ndarray *attr_value=NULL;
  for (i=1;i<argc;i++){
    if (option_is("-d","--data-set",argv[i])){
      dataset_name=strdup(argv[++i]);
    } else if (option_is("-a","--attribute",argv[i])){
      attr=strdup(argv[++i]);
    } else if (file_extension_is(argv[i],".h5")) {
      h5file=strdup(argv[i]);
    } else if (file_extension_is(argv[i],".txt")){
      txt=strdup(argv[i]);
    } else { //assume that it is a binary file
      infile=strdup(argv[i]);
    }
  }
  assert(h5file);
  assert(attr);
  assert(txt != infile);

  if (txt){
    attr_value=ndarray_from_text(txt);
  } else if (infile){
    attr_value=ndarray_from_binary(infile);
  } else {
    abort();
  }
  
  return EXIT_SUCCESS
}

