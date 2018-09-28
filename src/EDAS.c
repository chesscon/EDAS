#include <math.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>

#include "individuos.h"
#include "EDAS.h"

gsl_vector * array_to_vector_gsl(double *arr, int n) {
  gsl_vector *v = gsl_vector_alloc(n);
  for(int i=0; i < n; i++) {
    gsl_vector_set (v, i, arr[i]);
  }

  return v;
}

gsl_matrix * matrix_to_matrix_gsl(double **mat, int nr, int nc) {
  gsl_matrix *m = gsl_matrix_alloc(nr, nc);
  for(int i=0; i < nr; i++) {
    for(int j=0; j < nc; j++) {
      gsl_matrix_set(m, i, j, mat[i][j]); 
    }
  }
  return m; 
}

void calcula_medias(s_ind_dbl **pob, int tamPob, double *medias, int dim) {

  for(int i=0; i < dim; i++) {
    medias[i] = 0;
  }

  for(int i=0; i < tamPob; i++) {
    for(int j=0; j < dim; j++) {
      medias[j] += pob[i]->x[j];
    }
  }

//  printf("#\tMEDIAS: [");
  for(int i=0; i < dim; i++) {
    medias[i] /= tamPob;
//    printf(" %lf ", medias[i]);
  }
//  printf("]\n");

}

void calcula_varianzas(s_ind_dbl **pob, int tamPob, double *medias, 
  double *varianzas, int dim) {

  for(int i=0; i < dim; i++) {
    varianzas[i] = 0;
  }

  for(int i=0; i < tamPob; i++) {
    for(int j=0; j < dim; j++) {
      varianzas[j] += (pob[i]->x[j] - medias[j])*(pob[i]->x[j] - medias[j]);
    }
  }

//  printf("#\tVARIANZAS: [");
  for(int i=0; i < dim; i++) {
    varianzas[i] /= tamPob;
//    varianzas[i] += varianzas[i] < 1e-4 ? 1e-4 : 0;
//    printf(" %e ", varianzas[i]);
  }
//  printf("]\n");
}

void calcula_covarianzas(s_ind_dbl **pob, int tamPob, double *medias, 
  double **covarianzas, int dim) {

  for(int i=0; i < dim; i++) {
    for(int j=0; j <= i; j++) {
      covarianzas[i][j] = 0;
      covarianzas[j][i] = 0;
    }
  }

  for(int i=0; i < tamPob; i++) {
    for(int j=0; j < dim; j++) {
      for(int k=0; k <= j; k++ )
      covarianzas[j][k] += (pob[i]->x[j] - medias[j])*(pob[i]->x[k] - medias[k]);
    }
  }

  for(int i=0; i < dim; i++) {
    for(int j=0; j <= i; j++) {
      covarianzas[i][j] /= tamPob;
      covarianzas[j][i] = covarianzas[i][j];
    }
  }

}

void remuestrea_ind_normal(s_ind_dbl *ind, double *medias, 
  double *varianzas, const gsl_rng * r) {
  double media = 0;
  for(int i=0; i < ind->dim; i++) {
    //printf("#muestrando %d ..., media: %lf, var: %lf \n", i, medias[i], varianzas[i]);
    ind->x[i] = medias[i] + gsl_ran_gaussian(r, sqrt(varianzas[i]));
  }
}

void remuestrar_independientes(s_ind_dbl **pob, int tamPob, double *medias, 
  double *varianzas, const gsl_rng * r) {
  for(int i=0; i < tamPob; i++) {
    remuestrea_ind_normal(pob[i], medias, varianzas, r);
  }
}

void remuestrea_dim_ind_normal_cond(s_ind_dbl *ind, double *medias,
  double **covarianzas, int ind_var_indep, int ind_var_dep,  const gsl_rng * r) {

  double media_cond = medias[ind_var_dep] 
    + (covarianzas[ind_var_indep][ind_var_dep]/covarianzas[ind_var_indep][ind_var_indep])
      *(ind->x[ind_var_indep] - medias[ind_var_indep]);

  double var_con = covarianzas[ind_var_dep][ind_var_dep] 
    - (covarianzas[ind_var_indep][ind_var_dep]/covarianzas[ind_var_indep][ind_var_indep])
      *(covarianzas[ind_var_dep][ind_var_indep]);

  ind->x[ind_var_dep] = media_cond + gsl_ran_gaussian(r, sqrt(var_con));
}


void remuestrar_cadena(s_ind_dbl **pob, int tamPob, double *medias, 
  double **covarianzas, const gsl_rng * r) {

  int dim = pob[0]->dim;
  int permutacion[dim];
  for(int i=0; i < dim; i++) {
    permutacion[i] = i;
  }
  // Calcular informacion Mutua
  double IM[dim][dim];
  double rho=0;

  for(int i=0; i < dim; i++) {
    IM[i][i] = 0;
    for(int j=i+1; j < dim; j++) {
      rho = covarianzas[i][j]/(covarianzas[i][i]*covarianzas[j][j]);
      IM[i][j] = -0.5*log(1-rho);
      IM[j][i] = IM[i][j];
    }
  }

  // Obtener permutacion que maximiza informacion mutua
  // PENDIENTE

  // remuestrear poblacion
  for(int i=0; i < tamPob; i++) {

    pob[i]->x[permutacion[0]] = medias[permutacion[0]] 
      + gsl_ran_gaussian(r, sqrt(covarianzas[permutacion[0]][permutacion[0]] ));

    for (int j=1; j < dim; j++) {
      remuestrea_dim_ind_normal_cond(pob[i], medias, covarianzas, 
        permutacion[j-1], permutacion[j], r);
    }
  }
}

void remuestrar_multivariada(s_ind_dbl **pob, int tamPob, double *medias,
  double **covarianzas, const gsl_rng * r) {

  int resGSL = 0;

  // Crear vector de medias para GSL:
  gsl_vector * mu = array_to_vector_gsl(medias, pob[0]->dim);

  // Crear matriz de covacrianzas factorizada para GSL:
  gsl_matrix * L = matrix_to_matrix_gsl(covarianzas, pob[0]->dim, pob[0]->dim);
  resGSL = gsl_linalg_cholesky_decomp1(L);

  // Vector del resultado para la gsl:
  gsl_vector * result = gsl_vector_alloc(pob[0]->dim);

  for(int i=0; i < tamPob; i++) {
    resGSL = gsl_ran_multivariate_gaussian(r, mu, L, result);

    // copiar vector resultado al individuo
    for(int j=0; j < pob[i]->dim; j++ ) {
      pob[i]->x[j] = gsl_vector_get(result, j);
    }
  }

  gsl_vector_free(mu);
  gsl_vector_free(result);
  gsl_matrix_free(L);
}

s_resEDAc umda_c(s_paramsEDA params, double (*f)(double *, int)) {

  int tamPob = params.tamPob;
  int tamMuestra = (int) (params.porSel*tamPob);
  printf("#tamMuestra: %d \n", tamMuestra );

  // Generamos la poblacion inicial:
  s_ind_dbl **pob =  init_pob_dbl_rand(tamPob, params.dim, 
    params.ini, params.fin);

  int numEvals = 0;
  int indBest = evalua_poblacion(pob, tamPob, f);
  numEvals += tamPob;

  s_ind_dbl *best = clone_ind_dbl(pob[indBest]);
 
  int numIt = 1;
  double medias[params.dim];
  double varianzas[params.dim];

  printf("#Gen BEST\n");

  while (numEvals < params.maxEvals) {

    printf("%d %e \n", numIt, best->eval );

    // ordenamos la poblacion
    qsort(pob, tamPob, sizeof(s_ind_dbl *), cmp_eval_indp_dbl);
    indBest = 0;

    // Generarar modelo con los mejores individuos 
    calcula_medias(pob, tamMuestra, medias, params.dim);
    calcula_varianzas(pob, tamMuestra, medias, varianzas, params.dim);

    // Remuestrar con el modelo generado
    // Al best lo mantenemos en la poblacion
    remuestrar_independientes(pob+1, tamPob-1, medias, varianzas, params.r);

    // evaluar poblacion
    indBest = evalua_poblacion(pob+1, tamPob-1, f);
    numEvals += tamPob-1;

    if (cmp_eval_indp_dbl(&pob[indBest+1], &pob[indBest] ) < 0) {
      indBest += 1;
      copy_ind_dbl(pob[indBest], best);
    }

    numIt++;
  }

  libera_pob_dbl(pob, tamPob);

  s_resEDAc res;
  res.ind = best;
  res.totalEvals = numEvals;
  res.totalGens = numIt;

  return res;
}