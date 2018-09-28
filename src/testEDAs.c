#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "EDAS.h"
#include "individuos.h"
#include "funciones_objetivo.h"

s_paramsEDA read_params_EDA(char **argv);

int main(int argc, char **argv) {
 
  s_paramsEDA params = read_params_EDA(argv); 
  /* init generador de aleatorios */
  const gsl_rng_type * T;
  T = gsl_rng_default;
  gsl_rng_env_setup();

  // RAND STDLIB c
  unsigned int seed = time(NULL);
  printf("#SEED: %u \n", seed );
  srand(seed);

  // RAND GSL:
  params.r = gsl_rng_alloc (T);
  gsl_rng_set(params.r, seed);

/*  printf("# TEST RAND, GSL gaussian(0, 0.1): %lf \n",
    gsl_ran_gaussian(params.r, 0.1) ); */

  s_resEDAc res_umda = umda_c(params, params.f);

  printf("Mejor Individuo: \n");
  imprime_ind_dbl(res_umda.ind);

  libera_ind_dbl(res_umda.ind);

  gsl_rng_free (params.r);

  return 0;
}

s_paramsEDA read_params_EDA(char **argv) {
  s_paramsEDA params;
  printf("#PARAMS:\n");

  int numFun = atoi(argv[1]);
  params.dim = atoi(argv[2]);
  params.maxEvals = atoi(argv[3]);

  params.ini = atof(argv[4]);
  params.fin = atof(argv[5]);

  params.tamPob = atoi(argv[6]);
  params.porSel = atof(argv[7]);

  printf("#Num Fun: %d, dim: %d \n#\txi en [%lf, %lf]\n", 
    numFun, params.dim, params.ini, params.fin );
  printf("#Max Evals: %d \n", params.maxEvals);
  printf("#Tam Pob: %d, porcentaje seleccion: %lf \n", 
    params.tamPob, params.porSel );

  switch(numFun) {
    case ELLIPSOID:
      params.f = ellipsoid;
      break;
    case ZAKHAROV:
      params.f = zakharov;
      break;
    case ROSENBROCK:
      params.f = rosenbrock;
      break;
    case ACKLEY:
      params.f = ackley;
      break;
    case GRIEWANGK:
      params.f = griewangk;
      break;
    case RASTRIGIN:
      params.f = rastrigin;
      break;
    case SPHERE:
    default:
      params.f = sphere;
      break;
  }

  return params;
}