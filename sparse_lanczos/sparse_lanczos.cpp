#include "../all.h"

using namespace std;

void sparse_lanczos(int mat_dim, int mat_elements, int *row, int *col,
                    double *mat_val, double *eigen_value, double *eigen_vec,
                    std::string S_L_Outpufile_name)  // S_L = Sparse_Lanczos
{
    std::cout
        << "/******************************************************************"
           "****************SPARSE LANCZOS "
           "METHOD********************************************************/"
        << endl;

    ofstream of_S_L_Outputfile(S_L_Outpufile_name);
    if (!(of_S_L_Outputfile))
    {
        cerr << "Could not open the file -'" << S_L_Outpufile_name << "'"
             << endl;
    }

    int count = 0;
    double eps = 1.0;
    double err = 1.0e-15;
    bool err_checker = true;

    // setting Initial vector & standabilization
    double **u = new double *[mat_dim];
    for (int i = 0; i < mat_dim; i++)
    {
        u[i] = new double[mat_dim];
        vec_init(mat_dim, u[i]);
    }

    srand(time(NULL));
    for (int i = 0; i < mat_dim; i++)
    {
        u[0][i] = (double)rand() / RAND_MAX;
    }
    sdz(mat_dim, err, u[0]);

    // setting varaiables
    double *v = new double[mat_dim];
    // diagonal elements of tridiagonal matrix
    double *alpha = new double[mat_dim];
    vec_init(mat_dim, alpha);
    // sub diagonal elements of tridiagonal matrix
    double *beta = new double[mat_dim - 1];
    vec_init(mat_dim - 1, beta);

    /*arrays of eigen value*/
    double *eigen_value_even = new double[mat_dim];
    double *eigen_value_odd = new double[mat_dim];
    vec_init(mat_dim, eigen_value_even);
    vec_init(mat_dim, eigen_value_odd);

    /*arrays of ground state's eigen vector*/
    double *eigen_vec_even = new double[mat_dim];
    double *eigen_vec_odd = new double[mat_dim];
    vec_init(mat_dim, eigen_value_even);
    vec_init(mat_dim, eigen_value_odd);

    /*diag = alpha , sub_diag = beta*/
    double *diag = new double[mat_dim];
    double *sub_diag = new double[mat_dim - 1];
    vec_init(mat_dim, diag);
    vec_init(mat_dim - 1, sub_diag);

    for (int k = 0; k < mat_dim; k++)
    {
        vec_init(mat_dim, v);
        if (err_checker)
        {
            if (k == mat_dim - 1)
            {
                // calculate v = Au0(k)
                sparse_dgemv(mat_dim, mat_elements, v, row, col, mat_val, u[k]);
                // calculate alpha
                alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
            }
            else
            {  // calculate v[i] = Au0(k)
                sparse_dgemv(mat_dim, mat_elements, v, row, col, mat_val, u[k]);
                if (k == 0)
                {
                    alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
                    cblas_daxpy(mat_dim, -alpha[k], u[k], 1, v, 1);
                    beta[k] = cblas_dnrm2(mat_dim, v, 1);
                    cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
                    cblas_dcopy(mat_dim, v, 1, u[k + 1], 1);
                    sdz(mat_dim, err, u[k + 1]);
                }
                else
                {
                    alpha[k] = cblas_ddot(mat_dim, v, 1, u[k], 1);
                    cblas_daxpy(mat_dim, -beta[k - 1], u[k - 1], 1, v, 1);
                    cblas_daxpy(mat_dim, -alpha[k], u[k], 1, v, 1);
                    beta[k] = cblas_dnrm2(mat_dim, v, 1);
                    cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
                    cblas_dcopy(mat_dim, v, 1, u[k + 1], 1);
                    sdz(mat_dim, err, u[k + 1]);
                }
            }
            // calculate eigenvalue of A(k)
            cblas_dcopy(mat_dim, alpha, 1, diag, 1);
            cblas_dcopy(mat_dim - 1, beta, 1, sub_diag, 1);
            if (k % 2 == 0)
            {
                if (k == mat_dim - 1)
                {
                    LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, diag,
                                  sub_diag, eigen_vec_even, mat_dim);
                    cblas_dcopy(mat_dim, diag, 1, eigen_value_even, 1);
                }
                else
                {
                    LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, diag, sub_diag,
                                  eigen_vec_even, k + 2);
                    cblas_dcopy(mat_dim, diag, 1, eigen_value_even, 1);
                }
            }
            else
            {
                if (k == mat_dim - 1)
                {
                    LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, diag,
                                  sub_diag, eigen_vec_odd, mat_dim);
                    cblas_dcopy(mat_dim, diag, 1, eigen_value_odd, 1);
                }
                else
                {
                    LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, diag, sub_diag,
                                  eigen_vec_odd, k + 2);
                    cblas_dcopy(mat_dim, diag, 1, eigen_value_odd, 1);
                }
            }
            // check erros of k & k+1's eigenvalue of groundstate
            if (k > 0)
            {
                eps = abs(eigen_value_even[0] - eigen_value_odd[0]);
                if (eps > err)
                {
                    err_checker = true;
                }
                else if (eps < err)
                {
                    err_checker = false;
                    std::cout << "Break At Count : " << k << endl;
                }
            }
        }
        else
        {
            break;
        }
        count = k;
    }

    if (count % 2 == 0)
    {
        cblas_dcopy(mat_dim, eigen_value_even, 1, eigen_value, 1);
    }
    else
    {
        cblas_dcopy(mat_dim, eigen_value_odd, 1, eigen_value, 1);
    }

    ground_eigenvec(mat_dim, count, err, eigen_value[0], alpha, beta, u,
                    eigen_vec);

    /*Calculate ground state of eigen vector*/

    /*release memory*/
    for (int i = 0; i < mat_dim; i++)
    {
        delete[] u[i];
    }

    delete[] u;
    delete[] v;
    delete[] alpha;
    delete[] beta;
    delete[] eigen_value_even;
    delete[] eigen_value_odd;
    delete[] eigen_vec_even;
    delete[] eigen_vec_odd;
    delete[] diag;
    delete[] sub_diag;
}

// void sparse_lanczos(int mat_dim, int mat_elements, int *row, int *col,
//                     double *mat_val, double *eigen_value, double *eigen_vec,
//                     std::string S_L_Outpufile_name)  // S_L = Sparse_Lanczos
// {
//     std::cout
//         <<
//         "/******************************************************************"
//            "****************SPARSE LANCZOS "
//            "METHOD********************************************************/"
//         << endl;

//     ofstream of_S_L_Outputfile(S_L_Outpufile_name);
//     if (!(of_S_L_Outputfile))
//     {
//         cerr << "Could not open the file -'" << S_L_Outpufile_name << "'"
//              << endl;
//     }

//     int count = 0;
//     double eps = 1.0;
//     double err = 1.0e-16;
//     bool err_checker = true;

//     // setting Initial vector & standabilization
//     double **u = new double *[2];
//     for (int i = 0; i < 2; i++)
//     {
//         u[i] = new double[mat_dim];
//     }

//     srand(time(NULL));
//     for (int i = 0; i < mat_dim; i++)
//     {
//         u[0][i] = (double)rand() / RAND_MAX;
//         u[1][i] = 0.0;
//     }
//     sdz(mat_dim, err, u[0]);

//     // setting varaiables
//     double *v = new double[mat_dim];
//     // diagonal elements of tridiagonal matrix
//     double *alpha = new double[mat_dim];
//     vec_init(mat_dim, alpha);
//     // sub diagonal elements of tridiagonal matrix
//     double *beta = new double[mat_dim - 1];
//     vec_init(mat_dim - 1, beta);

//     /*arrays of eigen value*/
//     double *eigen_value_even = new double[mat_dim];
//     double *eigen_value_odd = new double[mat_dim];
//     vec_init(mat_dim, eigen_value_even);
//     vec_init(mat_dim, eigen_value_odd);

//     /*arrays of ground state's eigen vector*/
//     double *eigen_vec_even = new double[mat_dim];
//     double *eigen_vec_odd = new double[mat_dim];
//     vec_init(mat_dim, eigen_value_even);
//     vec_init(mat_dim, eigen_value_odd);

//     /*diag = alpha , sub_diag = beta*/
//     double *diag = new double[mat_dim];
//     double *sub_diag = new double[mat_dim - 1];
//     vec_init(mat_dim, diag);
//     vec_init(mat_dim - 1, sub_diag);

//     for (int k = 0; k < mat_dim; k++)
//     {
//         vec_init(mat_dim, v);
//         if (err_checker && k == mat_dim - 1)
//         {
//             if (k == mat_dim - 1)
//             {
//                 // calculate v = Au0(k)
//                 sparse_dgemv(mat_dim, mat_elements, v, row, col, mat_val,
//                 u[1]);
//                 // calculate alpha
//                 alpha[k] = cblas_ddot(mat_dim, v, 1, u[1], 1);
//             }
//             else
//             {
//                 if (k == 0)
//                 {
//                     // calculate v[i] = Au0(k)
//                     sparse_dgemv(mat_dim, mat_elements, v, row, col, mat_val,
//                                  u[0]);
//                     alpha[k] = cblas_ddot(mat_dim, v, 1, u[0], 1);
//                     cblas_daxpy(mat_dim, -alpha[k], u[0], 1, v, 1);
//                     beta[k] = cblas_dnrm2(mat_dim, v, 1);
//                     cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
//                     cblas_dcopy(mat_dim, v, 1, u[1], 1);
//                     sdz(mat_dim, err, u[1]);
//                 }
//                 else
//                 {
//                     sparse_dgemv(mat_dim, mat_elements, v, row, col, mat_val,
//                                  u[1]);
//                     alpha[k] = cblas_ddot(mat_dim, v, 1, u[1], 1);
//                     cblas_daxpy(mat_dim, -beta[k - 1], u[0], 1, v, 1);
//                     cblas_daxpy(mat_dim, -alpha[k], u[1], 1, v, 1);
//                     beta[k] = cblas_dnrm2(mat_dim, v, 1);
//                     cblas_dscal(mat_dim, 1.0 / beta[k], v, 1);
//                     cblas_dcopy(mat_dim, u[1], 1, u[0], 1);
//                     cblas_dcopy(mat_dim, v, 1, u[1], 1);
//                     sdz(mat_dim, err, u[1]);
//                 }
//             }
//             // calculate eigenvalue of A(k)
//             cblas_dcopy(mat_dim, alpha, 1, diag, 1);
//             cblas_dcopy(mat_dim - 1, beta, 1, sub_diag, 1);
//             if (k % 2 == 0)
//             {
//                 if (k == mat_dim - 1)
//                 {
//                     LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, diag,
//                                   sub_diag, eigen_vec_even, mat_dim);
//                     cblas_dcopy(mat_dim, diag, 1, eigen_value_even, 1);
//                 }
//                 else
//                 {
//                     LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, diag,
//                     sub_diag,
//                                   eigen_vec_even, k + 2);
//                     cblas_dcopy(mat_dim, diag, 1, eigen_value_even, 1);
//                 }
//             }
//             else
//             {
//                 if (k == mat_dim - 1)
//                 {
//                     LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', mat_dim, diag,
//                                   sub_diag, eigen_vec_odd, mat_dim);
//                     cblas_dcopy(mat_dim, diag, 1, eigen_value_odd, 1);
//                 }
//                 else
//                 {
//                     LAPACKE_dstev(LAPACK_ROW_MAJOR, 'N', k + 2, diag,
//                     sub_diag,
//                                   eigen_vec_odd, k + 2);
//                     cblas_dcopy(mat_dim, diag, 1, eigen_value_odd, 1);
//                 }
//             }
//             // check erros of k & k+1's eigenvalue of groundstate
//             if (k > 0)
//             {
//                 eps = abs(eigen_value_even[0] - eigen_value_odd[0]);
//                 if (eps > err)
//                 {
//                     err_checker = true;
//                 }
//                 else if (eps < err)
//                 {
//                     err_checker = false;
//                     std::cout << "Break At Count : " << k << endl;
//                 }
//             }
//         }
//         else
//         {
//             // break;
//         }
//         count = k;
//     }

//     if (count % 2 == 0)
//     {
//         cblas_dcopy(mat_dim, eigen_value_even, 1, eigen_value, 1);
//     }
//     else
//     {
//         cblas_dcopy(mat_dim, eigen_value_odd, 1, eigen_value, 1);
//     }

//     ground_eigenvec(mat_dim, count, err, eigen_value[0], alpha, beta, u,
//                     eigen_vec);

//     /*Calculate ground state of eigen vector*/

//     /*release memory*/
//     for (int i = 0; i < 2; i++)
//     {
//         delete[] u[i];
//     }

//     delete[] u;
//     delete[] v;
//     delete[] alpha;
//     delete[] beta;
//     delete[] eigen_value_even;
//     delete[] eigen_value_odd;
//     delete[] eigen_vec_even;
//     delete[] eigen_vec_odd;
//     delete[] diag;
//     delete[] sub_diag;
// }