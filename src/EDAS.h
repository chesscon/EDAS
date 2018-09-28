#ifndef EDAS_H
#define EDAS_H

#include "individuos.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

typedef struct {
  double ini;
  double fin;
  int dim;

  int tamPob;
  double porSel;
  int maxEvals;

  double (*f)(double *, int);

  // Generador de numeros aleatorias para la gsl
  gsl_rng * r;
} s_paramsEDA;

typedef struct {
  s_ind_dbl *ind;

  int totalEvals;
  int totalGens;

} s_resEDAc;

gsl_vector * array_to_vector_gsl(double *arr, int n);
gsl_matrix * matrix_to_matrix_gsl(double **mat, int n, int m);

void calcula_medias(s_ind_dbl **pob, int tamPob, double *medias, int dim);

void calcula_varianzas(s_ind_dbl **pob, int tamPob, double *medias, 
  double *varianzas, int dim);

void calcula_covarianzas(s_ind_dbl **pob, int tamPob, double *medias, 
  double **covarianzas, int dim);

void remuestrea_ind_normal(s_ind_dbl *ind, double *medias, 
  double *varianzas, const gsl_rng * r);

void remuestrar_independientes(s_ind_dbl **pob, int tamPob, double *medias, 
  double *varianzas, const gsl_rng * r);

void remuestrea_dim_ind_normal_cond(s_ind_dbl *ind, double *medias,
  double **covarianzas, int ind_var_indep, int ind_var_dep,  const gsl_rng * r);

void remuestrar_cadena(s_ind_dbl **pob, int tamPob, double *medias, 
  double **covarianzas, const gsl_rng * r);

void remuestrar_multivariada(s_ind_dbl **pob, int tamPob, double *medias,
  double **covarianzas, const gsl_rng * r);

s_resEDAc umda_c(s_paramsEDA params, double (*f)(double *, int));



#endif