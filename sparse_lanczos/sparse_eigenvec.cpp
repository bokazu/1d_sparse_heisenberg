#include "../all.h"

using namespace std;

void ground_eigenvec(int mat_dim, int count, double err,
                     double groundstate_eigen_value, double *alpha,
                     double *beta, double **u, double *eigen_vec)
{
    double *x = new double[mat_dim];
    vec_init(mat_dim, x);

    double *eigen = new double[mat_dim];
    vec_init(mat_dim, eigen);

    double lambda = groundstate_eigen_value;
    cout << "lambda = " << lambda << endl;
    x[0] = 1.0;
    x[1] = (lambda - alpha[0]) / beta[0];
    for (int i = 2; i < count; i++)
    {
        x[i] = ((lambda - alpha[i - 1]) * x[i - 1] - beta[i - 2] * x[i - 2]) /
               beta[i - 1];
    }
    for (int i = 0; i < count; i++)
    {
        cblas_daxpy(mat_dim, x[i], u[i], 1, eigen, 1);
    }
    sdz(mat_dim, err, eigen);
    cblas_dcopy(mat_dim, eigen, 1, eigen_vec, 1);
    cout << "lambda = " << lambda << endl;
    for (int i = 0; i < mat_dim; i++)
    {
        cout << "eigen vec[" << i << "] = " << eigen_vec[i] << endl;
    }

    delete[] x;
    delete[] eigen;
}