#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "individuos.h"

s_ind_dbl * genera_ind_dbl(int dim, double ini, double fin) {
  s_ind_dbl *ind = malloc(sizeof(s_ind_dbl));

  ind->x = malloc(sizeof(double)*dim);

  ind->ini = ini;
  ind->fin = fin;
  ind->dim = dim;

  ind->eval = DBL_MAX;

  return ind;
}

s_ind_dbl * genera_ind_dbl_rand(int dim, double ini, double fin) {
  s_ind_dbl *ind = genera_ind_dbl(dim, ini, fin);

  for(int i=0; i < dim; i++) {
    ind->x[i] = ini + (((double) rand())/RAND_MAX)*(fin-ini);
  }
//  printf("#individuo generado: \n");
//  imprime_ind_dbl(ind);

  return ind;
}

void copy_ind_dbl(s_ind_dbl *src, s_ind_dbl *dst) {
  //dst->ini = src->ini;
  //dst->fin = src->fin;
  dst->eval = src->eval;

  for(int i=0; i < src->dim; i++) {
    dst->x[i] = src->x[i];
  }
}

s_ind_dbl * clone_ind_dbl(s_ind_dbl *ind) {
  s_ind_dbl *res = genera_ind_dbl(ind->dim, ind->ini, ind->fin);
  copy_ind_dbl(ind, res);
  return res;
}

// Compara individuos almacenados en un arreglo
int cmp_eval_indp_dbl(const void * a, const void * b) {
  s_ind_dbl * ind_a = *((s_ind_dbl **) a);
  s_ind_dbl * ind_b = *((s_ind_dbl **) b);

  if (ind_a->eval > ind_b->eval ) return 1;
  if (ind_a->eval < ind_b->eval ) return -1;
  return 0;
}

void imprime_ind_dbl(s_ind_dbl *ind) {
  printf("eval: %e, [ ", ind->eval );

  for(int i=0; i < ind->dim; i++) {
    printf("%lf ", ind->x[i] );
  }
  printf("]\n");
}

void libera_ind_dbl(s_ind_dbl *ind) {
  free(ind->x);
  free(ind);
}

s_ind_dbl ** init_pob_dbl_rand(int tamPob, int dim, double ini, double fin) {
  s_ind_dbl **pob = malloc(sizeof(s_ind_dbl *)*tamPob);

  for(int i=0; i <tamPob; i++) {
    pob[i] = genera_ind_dbl_rand(dim, ini, fin);
  }

  return pob;
}

void libera_pob_dbl(s_ind_dbl **pob, int tamPob) {
  for(int i=0; i < tamPob; i++) {
    libera_ind_dbl(pob[i]);
  }
  free(pob);
}

double evalua_individuo(s_ind_dbl *ind, double (*f)(double *, int) ) {
  ind->eval = f(ind->x, ind->dim);

  return ind->eval;
}

/* devuelve el indice del individuo con la mejor evaluacion */
int evalua_poblacion(s_ind_dbl **pob, int tamPob, double(*f)(double *, int) ) {
  double bestEval = DBL_MAX;
  int bestInd = 0;

  double evalAct = DBL_MAX;

  for(int i=0; i < tamPob; i++) {
    evalAct = evalua_individuo(pob[i], f);
    if (cmp_eval_indp_dbl(&pob[i], &pob[bestInd]) < 0) {
      bestInd = i;
      bestEval = evalAct;
    }
  }
  return bestInd;
}