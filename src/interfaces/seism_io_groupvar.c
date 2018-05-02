/**
@section LICENSE
Copyright (c) 2013-2017, Regents of the University of California
All rights reserved.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdint.h>
#include"seism_io_globals.h"
//#include"seism_consts.h"
#include"seism_io_consts.h"
#ifndef _CALCREC_ONLY
  #ifndef _EXEC_SEISM_IO
    #include "mpi.h"
    #define _ADIOS
  #else   // _EXEC_SEISM_IO
    #define MPI_Comm int
    #define MPI_Barrier(comm) 0;
  #endif  // end _EXEC_SEISM_IO
  #ifdef _ADIOS
//    #include "adios.h"
//    #include "seism_io_adios.h"
  #endif

//#ifndef SEISM_CONSTS_H
//#define SEISM_CONSTS_H
// constants from seism_io.h
const int  DEBUG = 5;
const char const *opt[] = {"group","var"};
const int n_opt = 2;
//#endif

// constants
const char DEFAULT_TRANSPORT_METHOD[]="MPI";
const char EOS='\0';
const char intChar[]="integer";

Group *groups = NULL;
Group *lastGroup = NULL;

int ngroups = 0;

int groupID;

/* creates a new group
gID   : group ID integer. The index of the group
gname : group name
method: group method - ADIOS, MPI, HDF5, NETCDF, PNETCDF, etc.
*/
struct Group* newGroup(int gID, char *gname, char *method){
  int i;
  struct Group *g = (struct Group*)malloc(sizeof(struct Group));
  //g->name = (char*)malloc(sizeof(char)*(strlen(gname)+1));
  if(g==NULL){ printf("ERROR! Cannot allocate new Group!\n"); return NULL;}
  g->id = gID;
  for(i=0;i<strlen(gname);i++) g->name[i] = gname[i];
  strncpy(&(g->name[strlen(gname)]),&EOS,1);
  //strcpy(g->name, gname);
  if(!method) method = (char*)DEFAULT_TRANSPORT_METHOD;
  if(method){
    for(i=0;i<strlen(method);i++) g->method[i] = method[i];
    strncpy(&(g->method[strlen(method)]),&EOS,1);
  }
  //if(method) strcpy(g->method, method);
  else{
    if(DEBUG>0)
      printf("WARNING! Group %s transport method is not specified. Default (%s) is used.\n",g->name,DEFAULT_TRANSPORT_METHOD);
    strcpy(g->method,DEFAULT_TRANSPORT_METHOD);
  }
  if(gID==0) groups = g;
  if(lastGroup) lastGroup->next = g;
  lastGroup = g;
  g->next = NULL;
  g->firstVar = NULL;
  if(DEBUG>0)
    printf("new group created: name=%s method=%s\n",g->name,g->method);
  return g;
}

/* creates a new variable in the var-space
vname       : variable name - string
vtype       : variable type - string, e.g. integer, real*4
vdims       : variable dimensions - string. Local, per processor, e.g.
              "nxt,nyt,nzt" "1,3" "1,nxt"
globalDims  : global dimensions - string. the total dimensions for
              which this variable will exist in the file. e.g. "NX,NY,NZ"
              with all the writers, vdims must match globalDims
offsets     : offsets per processor - string. each processor writes
              a different portion of the same variable. offsets are in
              terms of indices, not bytes. e.g. "offsetx,offsety,offsetz"
              then offsetx can be defined with seism_define_var with a
              value coords[0]*nxt
attrKey     : some attribute key we want to save for the variable - string
              e.g. "info"
attrVal     : some attribute value related to the key, e.g. "dummyVar"
*/
struct Var* newVar(char *vname,char *vtype,char *vdims,
        char *globalDims,char *offsets, char *attrKey, char *attrVal){
  struct Var *tmpV;
  struct Var *v = (struct Var*)malloc(sizeof(struct Var));
  strcpy(v->name, vname);
  strcpy(v->type, vtype);
  v->val = NULL;
  if(vdims){
    v->dims = (char*)malloc(sizeof(char)*strlen(vdims));
    strcpy(v->dims, vdims);
  }
  else v->dims = NULL;
  if(globalDims){
    v->globalDims = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
    strcpy(v->globalDims, globalDims);
  }
  else v->globalDims = NULL;
  if(offsets){
    v->offsets = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
    strcpy(v->offsets, offsets);
  }
  else v->offsets = (char*)NULL;
  if(attrKey && attrVal){
    v->attrKey = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
    v->attrVal = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
    strcpy(v->attrKey, attrKey);
    strcpy(v->attrVal, attrVal);
  }
  else{ v->attrKey = NULL; v->attrVal = NULL;}
  v->next = NULL;
  v->group = lastGroup;
  if(!(lastGroup->firstVar)) lastGroup->firstVar = v;
  else{
    for(tmpV=lastGroup->firstVar;tmpV->next;tmpV=tmpV->next);
    tmpV->next = v;
  }
  v->globalDimsNumeric = NULL;
  if(DEBUG>0){
    printf("new var created: name=%s type=%s ",v->name,v->type);
    if(v->dims) printf("dims=%s ",v->dims);
    if(v->globalDims && v->offsets)
      printf("\tglobal-dims=%s offsets=%s",v->globalDims,v->offsets);
    if(v->attrKey && v->attrVal)
      printf("\tattr-key=%s attr-val=%s",v->attrKey,v->attrVal);
    printf("\n");
  }
  return v;
}

/* prints some information about the group object. mostly for
  debugging purposes.
*/
void _printGroup(struct Group *g){
  printf("%d) Group\n",mpi_rank);
  printf("%d) Group: name=%s method=%s\n\t%d) first-var=",mpi_rank,g->name,g->method,mpi_rank);
  if(g->firstVar) printf("%s next-group=",g->firstVar->name);
  else printf("NULL next-group=");
  if(g->next) printf("%s\n",g->next->name);
  else printf("NULL\n");
return;
}

/* prints some information about the variable object. mostly for
  debugging purposes.
*/
void _printVar(struct Var *v){
  printf("%d) Variable: name=%s type=%s dims=",mpi_rank,v->name,v->type);
  if(v->dims) printf("%s ",v->dims);
  else printf("<> ");
  printf("g-dims=");
  if(v->globalDims) printf("%s ",v->globalDims);
  else printf("<> ");
  printf("\tgroup=%s next-var=",v->group->name);
  if(v->next) printf("%s ",v->next->name);
  else printf("<> ");
  printf("val=");
  if(v->val) printf("%d\n",*(int*)(v->val));
  else printf("<>\n");
return;
}

/* findVarVal searches the variable object with name vname, and
  returns its value pointer
vname : variable name
TODO: modify so that returns pointer to variable, instead of its value
*/
void *_findVarVal(char *vname){
  struct Group *g;
  struct Var *v;
  void *val;
  val = NULL;
  if(DEBUG>1) printf("%d) Search for variable %s\n",mpi_rank,vname);
  if(DEBUG>1) printf("\t%d) Initial fast search...\n",mpi_rank);
  // most probably it is defined later and hence has its own
  // group with the same name.
  for(g=groups;g;g=g->next){
    if(DEBUG>2) _printGroup(g);
    if(!strncmp(g->name,vname,strlen(g->name))){
      val = (void*)(g->firstVar->val);
      if(DEBUG>1) printf("\t\t%d) FOUND!\n",mpi_rank);
      break;
    }
  }
  if(val) return val;
  // else, i.e. there is group with the same name, or this variable
  //  was defined in the initialization from the text parameters file
  else{
    if(DEBUG>1) printf("\n\t%d) Could not found. Now search all variables...\n",mpi_rank);
    // just linear search over groups and variables
    // TODO: efficient variable search can be done by creating variables
    //  in, say, binary search tree, or in alphabetical order, etc.
    for(g=groups;g;g=g->next){
      if(DEBUG>2) _printGroup(g);
      for(v=g->firstVar;v;v=v->next){
        if(DEBUG>2) _printVar(v);
        if(!strncmp(v->name,vname,strlen(v->name))){
          val = (v->val);
          if(DEBUG>1) printf("\t\t\t%d) FOUND!\n",mpi_rank);
          break;
        }
      }
      if(val) break;
    }
  }
  if(val) return val;
  else{
    printf("ERROR! %d) Value of the variable %s is not defined!\n",mpi_rank,vname);
    return NULL;
  }
return NULL;
}

/* split function - pretty standard string manipulation. searches
  line for the characters in delimiters, and returns indices in line
  where there is a delimiter, and delimTypes in the same order as
  indices. len is the number of delimiters found (length of inds and
  delimTypes)
line        : string, the actual line to be split
delim       : array of delimiter characters - e.g. " ,\n\t="
inds        : returned integer array of the indices in line where
              there is a delimiter
delimTypes  : list of delimiters such that line[inds[j]] = delimTypes[j]
len         : number of delimiters found - len=strlen(delimTypes)
*/
void _split(char *line, char *delim, int *inds, char *delimTypes, int *len){
  int i,j;
  int cnt = 0;
  for(i=0;i<strlen(line);i++){
    for(j=0;j<strlen(delim);j++)
      if(line[i] == delim[j]){
        inds[cnt] = i;
        delimTypes[cnt] = delim[j];
        cnt++;
      }
  }
  *len = cnt;
return;
}

/* resolves string. String comes from the dimensions. Hence it can be
  sth like: "nxt+4" , "2*PX-4" , "12/k2*3-a4"
  The function recognizes variable names, searches for the values and
  carries out the math and returns the dimension in integer
*/
int _resolveString(char *str){
  int len, i;
  int nums[MAX_NUMS_DELIMS];
  int numInd = 0;
  int tmpLen;
  int inds[MAX_NUMS_DELIMS];
  int lastInd = 0;
  char delimTypes[MAX_NUMS_DELIMS];
  char expr[MAX_CHAR_LEN];
  int res;

  if(DEBUG>1) printf("%d) Resolve string: \"%s\"\n",mpi_rank,str);
  _split(str, "*/-+", inds, delimTypes, &len);
  if(DEBUG>2) printf("%d) Resolve string deliminited len: %d\n",mpi_rank,len);
  // for each word delimited by the above math characters...
  for(i=0;i<len;i++){
    tmpLen = inds[i]-lastInd;
    strncpy(expr,&str[lastInd],tmpLen);
    strncpy(&expr[tmpLen],&EOS,1);
    // if expression/word starts with a digit, it has to be number
    // if is it sth like 5abc, we may get segmentation fault here!
    if((expr[0]>='0' && expr[0]<='9')){
      nums[numInd++] = atoi(expr);
    }
    // otherwise it is a variable name, find its value
    else{
      nums[numInd++] = *(int*)_findVarVal(expr);
    }
    lastInd = inds[i]+1;
  }
  // this is the last expression/word
  tmpLen = strlen(str)-lastInd;
  strncpy(expr,&str[lastInd],tmpLen);
  strncpy(&expr[tmpLen],&EOS,1);
  if(expr[0]>='0' && expr[0]<='9')
    nums[numInd++] = atoi(expr);
  else nums[numInd++] = *(int*)_findVarVal(expr);
  if(numInd!=len+1){
    printf("ERROR! %d) dim format is wrong: %s\n",mpi_rank,str);
    return -1;
  }
  // carry out math, first * and /, then later + and -
  for(i=0;i<len;i++){
    if(delimTypes[i]=='*'){
      nums[i] = nums[i]*nums[i+1];
      nums[i+1] = 0;
      delimTypes[i]='+';
    }
    else if(delimTypes[i]=='/'){
      nums[i] = nums[i]/nums[i+1];
      nums[i+1] = 0;
      delimTypes[i]='+';
    }
  }
  res = nums[0];
  for(i=0;i<len;i++){
    if(delimTypes[i]=='+'){
      res += nums[i+1];
    }
    else if(delimTypes[i]=='-'){
      res -= nums[i+1];
    }
  }
  if(DEBUG>1) printf("%d) resolved=%d\n",mpi_rank,res);

return res;
}

/* returns global dimensions and the group size
gname   : group name
retGdims: returned global dimensions of the last variable of the group
            - array of integers
ndim    : number of dimensions in returned global dimensions
gsize   : group size in bytes
*/
void _getGlobalDims(char *gname, int **retGdims, int *ndim, uint64_t *gsize){
  struct Group *g;
  struct Var *v;
  char *dims;
  char *numStr;
  size_t num;
  size_t tmpSize=0;
  size_t bsize;
  int lastInd, i;
  int gdims[MAX_DIM];
  dims = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
  // there is no group at all!
  if(!groups) printf("%d) there is no group!\n",mpi_rank);
  *gsize = 0;
  // search for the group
  for(g=groups;g;g=g->next){
    if(DEBUG>2) _printGroup(g);
    if(strlen(gname)!=strlen(g->name)) continue;
    // found the group!
    if(!strncmp(gname,g->name,strlen(gname))){
      if(DEBUG>1) printf("%d) group found: gname=%s\n",mpi_rank,gname);
        //_printGroup(g);
      // traverse each variable
      for(v=g->firstVar;v;v=v->next){
        // dimensions are returned only for the last variable in a group
        *ndim = 0;
        if(DEBUG>1)
          _printVar(v);
        if(!strncmp(v->type,"integer",strlen(v->type)))
          bsize = sizeof(int);
        else if(!strncmp(v->type,"integer*8",strlen(v->type)))
          bsize = 8;
        else if(!strncmp(v->type,"integer*4",strlen(v->type)))
          bsize = 4;
        else if(!strncmp(v->type,"real",strlen(v->type)))
          bsize = sizeof(float);
        else if(!strncmp(v->type,"real*8",strlen(v->type)))
          bsize = 8;
        else if(!strncmp(v->type,"real*4",strlen(v->type)))
          bsize = 4;
        else if(!strncmp(v->type,"float",strlen(v->type)))
          bsize = sizeof(float);
        else if(!strncmp(v->type,"double",strlen(v->type)))
          bsize = sizeof(double);
        else if(!strncmp(v->type,"char",strlen(v->type)))
          bsize = sizeof(char);
        if(DEBUG>1) printf("\t%d) bsize=%ld\n",mpi_rank,bsize);
        // if there is already numerical dimension info, save it
        if(v->globalDimsNumeric){
          if(DEBUG>1) printf("\t%d) globalDimsNumeric\n",mpi_rank);
          if(DEBUG>2) printf("\t%d) globalDimsNumeric: %s\n",mpi_rank,v->globalDimsNumeric);
          strcpy(dims,v->globalDimsNumeric);
          if(DEBUG>2) printf("\t%d) globalDimsNumeric is copied to dims\n",mpi_rank);
          tmpSize = 1;
          lastInd = 0;
          // for each comma-separated dimension, convert to int
          for(i=0;i<strlen(dims);i++){
            if(dims[i]==','){
              numStr = (char*)malloc(sizeof(char)*(i-lastInd+1));
              strncpy(numStr,&dims[lastInd],i-lastInd);
              strncpy(&numStr[i-lastInd],&EOS,1);
              num = atoi(numStr);
              free(numStr);
              lastInd = i+1;
              (*retGdims)[*ndim] = num;
              *ndim += 1;
              tmpSize *= num;
            }
          }
          // this is the last dimension
          numStr = (char*)malloc(sizeof(char)*(i-lastInd+1));
          strncpy(numStr,&dims[lastInd],i-lastInd);
          strncpy(&numStr[i-lastInd],&EOS,1);
          num = atoi(numStr);
          free(numStr);
          lastInd = i+1;
          (*retGdims)[*ndim] = num;
          *ndim += 1;
          tmpSize *= num;
          tmpSize *= bsize;
          *gsize += tmpSize;
          if(DEBUG>1) printf("\t%d) tmpSize=%ld\n",mpi_rank,tmpSize);
          continue;
        }
        // else if global dimensions exist but in string,
        //  that may not be numerical
        if(v->globalDims){
          if(DEBUG>1) printf("\t%d) globaldims\n",mpi_rank);
          if(DEBUG>1) printf("\t%d) globaldims is=%s\n",mpi_rank,v->globalDims);
          strcpy(dims,v->globalDims);
        }
        else if(v->dims){
          if(DEBUG>1) printf("\t%d) dims\n",mpi_rank);
          if(DEBUG>1) printf("\t%d) dims is=%s\n",mpi_rank,v->dims);
          strcpy(dims,v->dims);
        }
        else{
          if(DEBUG>1) printf("\t%d) dims is set to 1\n",mpi_rank);
          strcpy(dims,"1");
        }
        // we copied global dims, or dims, or 1 (if no dims info
        //  is given) of the variable to string dims
        tmpSize = 1;
        lastInd = 0;
        // for each comma-separated dimension, resolve it
        for(i=0;i<strlen(dims);i++){
          if(dims[i]==','){
            numStr = (char*)malloc(sizeof(char)*(i-lastInd+1));
            strncpy(numStr,&dims[lastInd],i-lastInd);
            strncpy(&numStr[i-lastInd],&EOS,1);
            // resolves string and returns integer
            num = _resolveString(numStr);
            free(numStr);
            lastInd = i+1;
            //printf("%d) just before gdims\n",mpi_rank);
            (*retGdims)[*ndim] = num;
            //printf("%d) just after gdims\n",mpi_rank);
            *ndim += 1;
            tmpSize *= num;
          }
        }
        // the last dimension
        numStr = (char*)malloc(sizeof(char)*(i-lastInd+1));
        strncpy(numStr,&dims[lastInd],i-lastInd);
        strncpy(&numStr[i-lastInd],&EOS,1);
        num = _resolveString(numStr);
        free(numStr);
        lastInd = i+1;
        //printf("%d) just before gdims, after loop of dims\n",mpi_rank);
        (*retGdims)[*ndim] = num;
        //printf("%d) just after gdims, after loop of dims\n",mpi_rank);
        *ndim += 1;
        tmpSize *= num;
        tmpSize *= bsize;
        *gsize += tmpSize;
        if(DEBUG>1)
          printf("\t%d) tmpSize=%ld\n",mpi_rank,tmpSize);
      }
      break;
    }
  }
  free(dims);
  //*retGdims = gdims;
  //if(mpi_rank==1) printf("returned global dims=");
  //if(mpi_rank==1) for(i=0;i<*ndim;i++) printf("%d ",(*retGdims)[i]);
  //if(mpi_rank==1) printf("\n");

return;
}

/* processes the line: finds keywords, saves new group and variable
  objects
line    : given line
optID   : index of the option (first word in a line: group/var)
*/
void _processLine(char *line, int optID){

  //struct Group* g;
  //struct Var* v;

  char *name, *method;
  char *vname, *vtype, *vdims, *vglobalDims, *voffsets, *vattrKey, *vattrVal;

  int len, i, j, ind;
  int inds[MAX_NUMS_DELIMS];
  char delimTypes[MAX_NUMS_DELIMS];

  // first split the line with delimiters
  _split(line, " ,=\t\n", inds, delimTypes, &len);
  if(DEBUG>2)
    for(i=0;i<len;i++){
      printf("%d) ind:%d, delimType:%c\n",i,inds[i],delimTypes[i]);
    }

  switch(optID){
    case 0: // if GROUP
            name = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
            method = NULL;
            for(i=0;i<len;i++){
              // if =, find the keyword "name," "method"
              if(delimTypes[i] == '='){
                if(!strncmp(&line[inds[i-1]+1],"name",4)){
                  strncpy(name,&line[inds[i]+1],inds[i+1]-inds[i]-1);
                  strncpy(&name[inds[i+1]-inds[i]-1],&EOS,1);
                  //printf("name=%s\n",name);
                }
                else if(!strncmp(&line[inds[i-1]+1],"method",6)){
                  method = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  strncpy(method,&line[inds[i]+1],inds[i+1]-inds[i]-1);
                  strncpy(&method[inds[i+1]-inds[i]-1],&EOS,1);
                }
              }
            }
            //printf("New group: %s\n",name);
            //g = newGroup(++groupID,name,method);
            newGroup(++groupID,name,method);
            //printf("Group is created!\n");
            free(name);
            if(method) free(method);
            break;
    case 1: // if VAR
            name = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
            vtype = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
            vdims = NULL;
            vglobalDims = NULL;
            voffsets = NULL;
            vattrKey = NULL;
            vattrVal = NULL;
            for(i=0;i<len;i++){
              if(delimTypes[i] == '='){
                if(!strncmp(&line[inds[i-1]+1],"name",4)){
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(name,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&name[inds[j]-inds[i]-1],&EOS,1);
                  //printf("name=%s\n",name);
                  i = j-1;
                }
                else if(!strncmp(&line[inds[i-1]+1],"type",4)){
                  strncpy(vtype,&line[inds[i]+1],inds[i+1]-inds[i]-1);
                  strncpy(&vtype[inds[i+1]-inds[i]-1],&EOS,1);
                  //printf("vtype=%s\n",vtype);
                }
                else if(!strncmp(&line[inds[i-1]+1],"dims",4)){
                  vdims = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(vdims,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&vdims[inds[j]-inds[i]-1],&EOS,1);
                  //printf("vdims=%s, j=%d-%d, inds[j]=%d-%d\n",vdims,j,i,inds[j],inds[i]);
                  i = j-1;
                }
                else if(!strncmp(&line[inds[i-1]+1],"global-dims",11)){
                  vglobalDims = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(vglobalDims,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&vglobalDims[inds[j]-inds[i]-1],&EOS,1);
                  i = j-1;
                }
                else if(!strncmp(&line[inds[i-1]+1],"offsets",7)){
                  voffsets = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(voffsets,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&voffsets[inds[j]-inds[i]-1],&EOS,1);
                  i = j-1;
                }
                else if(!strncmp(&line[inds[i-1]+1],"attr-key",8)){
                  vattrKey = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(vattrKey,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&vattrKey[inds[j]-inds[i]-1],&EOS,1);
                  i = j-1;
                }
                else if(!strncmp(&line[inds[i-1]+1],"attr-val",8)){
                  vattrVal = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                  for(j=i+1;j<len;j++){
                    if(delimTypes[j]==',') continue;
                    else break;
                  }
                  strncpy(vattrVal,&line[inds[i]+1],inds[j]-inds[i]-1);
                  strncpy(&vattrVal[inds[j]-inds[i]-1],&EOS,1);
                  i = j-1;
                }
              }
            }
            /*_readVals(line, "type=", "", names, &len);
            strcpy(vtype,names[0]);
            _readVals(line, "dims=", "", names, &len);
            strcpy(vdims,names[0]);
            _readVals(line, "name=", ",", names, &len);
            for(i=0;i<len;i++) printf("vname=%s\t",names[i]);
            */
            //printf("vtype=%s, vdims=%s\n",vtype,vdims);
            ind = 0;
            vname = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
            // the name value can be comma-separated names,
            //  then for each name we need to create a new Var
            for(j=0;j<strlen(name);j++){
              if(name[j]==','){
                strncpy(&vname[ind],&EOS,1);
                //printf("New var: name=%s, type=%s, dims=%s\n",vname,vtype,vdims);
                newVar(vname, vtype, vdims, vglobalDims, voffsets,
                      vattrKey, vattrVal);
                free(vname);
                vname = (char*)malloc(sizeof(char)*MAX_CHAR_LEN);
                ind = 0;
              }
              else{
                vname[ind++]=name[j];
              }
            }
            strncpy(&vname[ind],&EOS,1);
            //printf("New var: name=%s, type=%s, dims=%s\n",vname,vtype,vdims);
            newVar(vname, vtype, vdims, vglobalDims, voffsets,
                  vattrKey, vattrVal);
            /*for(j=0;j<MAX_CHAR_LEN;j++){
              name[j]='\0';
              vname[j]='\0';
              vtype[j]='\0';
              vdims[j]='\0';
            }*/
            free(name);
            free(vname);
            free(vtype);
            if(vdims) free(vdims);
            if(vglobalDims) free(vglobalDims);
            if(voffsets) free(voffsets);
            if(vattrKey) free(vattrKey);
            if(vattrVal) free(vattrVal);
            break;
  }

return;
}

/* reads in SEISM parameter input text file and populates group
  -variable space with objects.
   MUST BE CALLED BY ONLY 1 PROCESSOR!
fname : input file name
*/
void read_SEISM_params(char *fname){

  FILE *fin;
  // reads in line by line
  char line[MAX_CHAR_LEN];

  int i, lineStartInd;

  // keeps track of number of groups too
  groupID = -1;

  fin = fopen(fname, "r");
  if(fin == NULL){
    printf("FILE Cannot be found: %s\n",fname);
    return;
  }
  // For each line...
  while(fgets(line, sizeof(line), fin)!=NULL){
    lineStartInd = -1;
    // on each line, pass over initial space characters
    // if there is nothing but spaces and newline,
    //  just skip that line below.
    for(i=0;i<strlen(line);i++){
      //printf("(%c,%d)",line[i],(int)line[i]);
      if(line[i]==' ') continue;
      else if(line[i]!='\n'){
        lineStartInd = i;
        break;
      }
      else break;
    }
    // if line is empty (spaces+newline) skip to next line
    if(lineStartInd == -1){ if(DEBUG>1) printf("-- empty line read\n"); continue;}
    // if line starts with #, skip it since it is a comment line
    if(line[lineStartInd]=='#'){
      if(DEBUG>1)
        printf("Comment: >%s",line);
      continue;    // comments
    }
    if(DEBUG>1)
      printf("Line read: >%s",line);
    // if comes here, the line is valid and should be processed
    // check the first words of the line for each option,
    // if found process, else continue. If the option is unknown,
    // skip the line
    for(i=0;i<n_opt;i++)
      if(!strncmp(&line[lineStartInd],opt[i],strlen(opt[i])))
        _processLine(&line[lineStartInd], i);
  }

return;
}


/* init function to initialize the SEISM lib. Reads in text input file
  and creates the initial group and variable objects.
comm    : MPI_Comm pointer for initialization. ADIOS 1.3, 1.4 requires
            all processors to call init, not only the writers
rank    : MPI_rank of the caller (int* because fortran interface passes
            the pointer)
coords  : 3D vector of 0-based processor indices in the virtual topology
            (0,0,0),(1,0,0)...(PX-1,PY-1,PZ-1)
err     : error code
*/
void test_seism_init_(MPI_Comm *comm, int *rank, int **coords, int *err){
  if(DEBUG>0)
    printf("%d) SEISM_IO Init...\n",*rank);
  mpi_rank = *rank;
#ifdef _ADIOS
  seism_init_adios(comm, rank, coords, err);
#endif
  if(DEBUG>0)
    printf("%d) SEISM_IO Init returns\n",*rank);
return;
}

/* defines a new variable by creating/modifying Var object.
  If there is already a variable object with vname, only
  val is added to that object. If not, a new group is created
  with group name set as vname, and its firstVar (first
  variable) is set as a newly created Var object.
vname   : variable name, string, to be defined
val     : void* to the variable itself
err     : error code
vnlen   : length of the variable name
*/
void seism_define_var_(char *vname, void *val, int *err, unsigned int vnlen){

  struct Group *g;
  struct Var *v;
  // checks if in the search, vname belongs to a Var object
  int found=0;
  // C version of vname
  char *cvname;
  // just a char* to NULL
  char *nullChar;
  char *vtype;
  char *vdims;

  nullChar = NULL;
  cvname = (char*)malloc(sizeof(char)*(vnlen+1));
  strncpy(cvname,vname,(vnlen));
  strncpy(&cvname[vnlen],&EOS,1);

  if(DEBUG>0){
    char seismChar[] = "seism_";
    if(!strncmp(cvname,seismChar,6))
      printf("%d) Define variable: name=%s of len=%u\n",mpi_rank,cvname,vnlen);
    else
      printf("%d) Define variable: name=%s of len=%u val=%d\n",mpi_rank,cvname,vnlen,*(int*)val);
  }
  // linear search over all the groups and variables for vname
  // TODO: write a variable search function and call it here.
  for(g=groups;g;g=g->next){
    for(v=g->firstVar;v;v=v->next){
      //if(DEBUG>1) printf("%d) comparing v->name %s to vname %s\n",mpi_rank,v->name,vname);
      if(!strcmp(v->name,cvname)){
        v->val = val;
        found = 1;
        if(DEBUG>1)
          printf("%d) Found: group=%s var=%s\n",mpi_rank,g->name,v->name);
        break;
      }
    }
    // if found, no need to continue to search
    if(found) break;
  }
  // if no such variable, create a group and variable
  if(!found){
    if(DEBUG>1)
      printf("%d) Not found, create new group\n",mpi_rank);
    g = newGroup(++groupID,cvname,(char*)nullChar);
/*
    g = (struct Group*)malloc(sizeof(struct Group));
    strcpy(g->name,vname);
    strcpy(g->method,DEFAULT_TRANSPORT_METHOD);
    //printf("%d) new group name and method are copied.\n",mpi_rank);
    if(!lastGroup) printf("%d) lastGroup is NULL!\n",mpi_rank);
    lastGroup->next = g;
    lastGroup = g;
    g->next = NULL;
    g->id = ++groupID;
    //printf("%d) new group created\n",mpi_rank);
*/
    vtype = (char*)malloc(sizeof(char)*8);
    strcpy(vtype,intChar);
    strncpy(&vtype[7],&EOS,1);
    v = newVar(cvname,vtype,(char*)nullChar,
      (char*)nullChar,(char*)nullChar,(char*)nullChar,(char*)nullChar);
    v->val = val;
/*
    v = (struct Var*)malloc(sizeof(struct Var));
    strcpy(v->name,vname);
    //printf("%d) new var name copied\n",mpi_rank);
    v->val = val;
    //printf("%d) new var val set\n",mpi_rank);
    v->group = g;
    //printf("%d) new var group set\n",mpi_rank);
    g->firstVar = v;
    //printf("%d) new var's group's firstVar is set\n",mpi_rank);
    v->next = NULL;
    //printf("%d) new var created!\n",mpi_rank);
*/
  }

return;
}

/* file open function fortran interface
seism_f : file pointer to be returned to be caller after
            successful open
gname   : group name to open the file for
fname   : file name of the file
fmode   : file mode = "w" "r" "a" as write-read-append
comm    : MPI_Comm pointer to MPI communicator
err     : error code pointer. 0 is no error.
gnamelen: length of the group name
fnamelen: length of the file name
fmodelen: length of the file mode
*/
void test_seism_file_open_(int64_t *seism_f, char *gname,
    char *fname, char *fmode, MPI_Comm *comm, int *err,
    unsigned int gnamelen, unsigned int fnamelen, unsigned int fmodelen){
  // groupsize
  uint64_t gsize;
  // C version gname, fname, fmode. They will be copied from
  //  Fortran names
  char *cgname, *cfname, *cfmode;
  cgname = (char*)malloc(sizeof(char)*(gnamelen+1));
  cfname = (char*)malloc(sizeof(char)*(fnamelen+1));
  cfmode = (char*)malloc(sizeof(char)*(fmodelen+1));
  strncpy(cgname,gname,gnamelen);
  strncpy(&cgname[gnamelen],&EOS,1);
  strncpy(cfname,fname,fnamelen);
  strncpy(&cfname[fnamelen],&EOS,1);
  strncpy(cfmode,fmode,fmodelen);
  strncpy(&cfmode[fmodelen],&EOS,1);
  if(DEBUG>0) printf("%d) Opening file for group %s\n",mpi_rank,cgname);
#ifdef _ADIOS
  _calcGroupSize(cgname, &gsize);
  if(DEBUG>0) printf("%d) groupsize is calculated: %u\n",mpi_rank,gsize);
  if(DEBUG>0) printf("%d) calling adios file open: %s\n",mpi_rank,cfname);
  seism_file_open_adios(seism_f, cgname, cfname, cfmode, comm, err);
  if(DEBUG>0) printf("%d) adios file open returned %d\n",mpi_rank,*err);
#endif
  free(cgname);
  free(cfname);
  free(cfmode);
return;
}


/* write function fortran interface
seism_f : int64_t pointer file pointer. File has to be opened before.
vname   : variable name, that is going to be written to the file
err     : error code, integer pointer. 0 is no error.
vnlen   : length of the variable name (for fortran interface, since
            it returns the variable name as well
*/
void test_seism_write_(int64_t *seism_f, char *vname, int *err, unsigned int vnlen){
  void *val;
  // variable name in C. we will copy from fortran name.
  char *cvname;
  // allocate 1 more char to copy EOS (end of string) to the end
  cvname = (char*)malloc(sizeof(char)*(vnlen+1));
  strncpy(cvname,vname,vnlen);
  strncpy(&cvname[vnlen],&EOS,1);
  if(DEBUG>0){
    printf("%d) Writing to file variable %s\n",mpi_rank,cvname);
  }
#ifdef _ADIOS
  seism_write_adios(seism_f, cvname, err);
#endif
return;
}

/* seism_file_close closes an already opened file
seism_f : pointer to file
err     : error code
*/
void test_seism_file_close_(int64_t *seism_f, int *err){
#ifdef _ADIOS
  *err = adios_close(*seism_f);
#endif
return;
}

/* seism_finalize finalizes the library
err : error code
TODO: free all the memory used by the library! delete group-var space
*/
void test_seism_finalize_(int *err){
#ifdef _ADIOS
  *err = adios_finalize(mpi_rank);
#endif
return;
}
#endif

/*
calculates the receiver indices and number
n is the number of points per process
ind is the processor index in this dimension (0-based)
nbg,ned,nskp are the beginning,end,skp values of receivers
npbg,nped are the beginning and end indices of local receivers
nprec is the number of local receivers
*/
void _calcBgEdRec(int n, int ind, int nbg, int ned, int nskp,
    int *npbg, int *nped, int *nprec)
{

  if(nbg > n*(ind+1))     *nprec = 0;
  else if(ned < n*ind+1)  *nprec = 0;
  else{
    if(n*ind >= nbg)       *npbg = (nbg-n*ind-1)%nskp+1;
    else                   *npbg = nbg-n*ind;
    if(n*(ind+1) <= ned)   *nped = n;
    else                   *nped = ned-n*ind;
    if(*npbg<0) *npbg += nskp;
    *nped -= (*nped-*npbg)%nskp;
    *nprec = (*nped-*npbg)/nskp+1;
  }
return;
}

/*
seism_calc_rec calculates and counts the receivers' local indices
processor indices are 0-based because of the MPI calls
all other indices are 1-based (in GPU code npbgx was 0-based!)
*/
void seism_calc_rec_(int *nxt, int *nyt, int *nzt,
  int *nbgx, int *nedx, int *nskpx,
  int *nbgy, int *nedy, int *nskpy,
  int *nbgz, int *nedz, int *nskpz,
  int *npbgx, int *npedx, int *npbgy, int *npedy, int *npbgz, int *npedz,
  int *nprecx, int *nprecy, int *nprecz,
  int *coords){

#ifdef _CALCREC_ONLY
  printf("in seism_calc_rec\n");
#endif
  // ind is processor index - 0-based
  int indx = coords[0];
  int indy = coords[1];
  int indz = coords[2];

  *npbgx = 1;
  *npbgy = 1;
  *npbgz = 1;
  *npedx = (*nxt)*(indx+1);
  *npedy = (*nyt)*(indy+1);
  *npedz = (*nzt)*(indz+1);

#ifdef _CALCREC_ONLY
  printf("call for x\n");
#endif
  _calcBgEdRec(*nxt,indx,*nbgx,*nedx,*nskpx,npbgx,npedx,nprecx);
#ifdef _CALCREC_ONLY
  printf("now call for y\n");
#endif
  _calcBgEdRec(*nyt,indy,*nbgy,*nedy,*nskpy,npbgy,npedy,nprecy);
  _calcBgEdRec(*nzt,indz,*nbgz,*nedz,*nskpz,npbgz,npedz,nprecz);

  if(*nprecx==0 || *nprecy==0 || *nprecz==0){
    *nprecx = 0;
    *nprecy = 0;
    *nprecz = 0;
  }

return;
}

/* test function for receiver calculation only
*/
#ifdef _CALCREC_ONLY
int main(){
  int nxt=100;
  int nyt=120;
  int nzt=172;
  int nbgx=1;
  int nedx=522;
  int nskpx=7;
  int nbgy=133;
  int nedy=180;
  int nskpy=3;
  int nbgz=175;
  int nedz=175;
  int nskpz=2;
  int npbgx,npbgy,npbgz,npedx,npedy,npedz;
  int nprecx,nprecy,nprecz;
  int coord[]={2,1,1};

  seism_calc_rec_(&nxt,&nyt,&nzt,&nbgx,&nedx,&nskpx,
      &nbgy,&nedy,&nskpy,&nbgz,&nedz,&nskpz,
      &npbgx,&npedx,&npbgy,&npedy,&npbgz,&npedz,
      &nprecx,&nprecy,&nprecz,coord);

  printf("np: (%d,%d)\t(%d,%d)\t(%d,%d)\n",npbgx,npedx,npbgy,npedy,npbgz,npedz);
  printf("nprec: %d,%d,%d\n",nprecx,nprecy,nprecz);
return 0;
}
#endif

/* test function for init, variable definitions, file open/close etc.
*/
#ifdef _EXEC_SEISM_IO
int main(){

  size_t num;
  int err;
  int v1 = 100;
  int v2 = 1000;
  char lang[] = "Fortran";
  read_SEISM_params("./seism-params.txt");
  if(DEBUG>1)
    printf("\n\nGenerating adios configuration file...\n");
  writeAdiosXML("./adios-config.xml");
  DUMP(DEBUG);
  printf("\nAdios file is written\n");
  seism_define_var_("nxt",&v1,&err,3);
  seism_define_var_("nyt",&v1,&err,3);
  seism_define_var_("nzt",&v1,&err,3);
  printf("define NX as 1000\n");
  seism_define_var_("NX",&v2,&err,2);
  seism_define_var_("NY",&v2,&err,2);
  seism_define_var_("NZ",&v2,&err,2);
  printf("find var val for NZ\n");
  num = *(int*)_findVarVal("NZ");
  printf("NX=%ld\n",num);
  _calcGroupSize("output-vel-x",&num);
  printf("Group size vel-x=%ld\n",num);
  _calcGroupSize("globals",&num);
  printf("Group size globals=%ld\n",&num);

return 0;
}
#endif
