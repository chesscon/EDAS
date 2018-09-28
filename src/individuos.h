#ifndef INDIVIDUOS_H
#define INDIVIDUOS_H

typedef struct {
  double *x;
  int dim;

  double ini;
  double fin;

  double eval;
} s_ind_dbl;

//typedef **s_ind_dbl pob_dbl;

s_ind_dbl * genera_ind_dbl(int dim, double ini, double fin);

s_ind_dbl * genera_ind_dbl_rand(int dim, double ini, double fin);

void copy_ind_dbl(s_ind_dbl *src, s_ind_dbl *dst);

s_ind_dbl * clone_ind_dbl(s_ind_dbl *ind);

int cmp_eval_indp_dbl(const void * a, const void * b);

void imprime_ind_dbl(s_ind_dbl *ind);

void libera_ind_dbl(s_ind_dbl *ind);

s_ind_dbl ** init_pob_dbl_rand(int tamPob, int dim, double ini, double fin);

void libera_pob_dbl(s_ind_dbl **pob, int tamPob);

double evalua_individuo(s_ind_dbl *ind, double (*f)(double *, int) );

/* devuelve el indice del individuo con la mejor evaluacion */
int evalua_poblacion(s_ind_dbl **pob, int tamPob, double(*f)(double *, int) );

#endif