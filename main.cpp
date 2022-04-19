#include "all.h"

using namespace std;

int main()
{
    cout << "START 1D Heisenberg Model Calculation(Sparse)\n\n";

    /*Calculate Hamiltonian Matrix's Components*/
    /*INPUT MODEL DATA*/
    string M_H_setting_file_name = "./model_set/settingfile.txt";
    ifstream if_M_H_Settingfile(M_H_setting_file_name);
    if (!(if_M_H_Settingfile))
    {
        cerr << "Could not open the file'" << M_H_setting_file_name << "'"
             << endl;
    }

    string M_H_OutputFile_name, M_H_JsetFile_name, Boundary_Condition,
        D_L_OutputFile_name;

    int tot_site_num, precision;
    get_data(if_M_H_Settingfile, tot_site_num, M_H_OutputFile_name,
             M_H_JsetFile_name, D_L_OutputFile_name, Boundary_Condition,
             precision);

    std::cout << "/************************************************************"
                 "***************************"
              << "INPUT DATA"
              << "*************************************************************"
                 "**************************/\n";
    std::cout << "tot_site_num        = " << tot_site_num << endl;
    std::cout << "M_H_OutputFile_name = " << M_H_OutputFile_name << endl;
    std::cout << "M_H_JsetFile_name   = " << M_H_JsetFile_name << endl;
    std::cout << "D_L_OutputFile_name = " << D_L_OutputFile_name << endl;
    std::cout << "Boundary Condition  = " << Boundary_Condition << endl;
    std::cout << "precision           = " << precision << endl;
    std::cout << "/************************************************************"
                 "***************************"
              << "*************************************************************"
                 "************************************************/\n";

    /*******************************************************************************************************************************/
    int mat_dim = 1 << tot_site_num;

    int mat_elements = 0;
    mat_elements = sparse_count_mat_elements(
        mat_dim, tot_site_num, M_H_JsetFile_name, Boundary_Condition);
    cout << "mat_elements = " << mat_elements << endl;
    int *row = new int[mat_elements];
    int *col = new int[mat_elements];
    double *mat_val = new double[mat_elements];

    vec_init(mat_elements, row);
    vec_init(mat_elements, col);
    vec_init(mat_elements, mat_val);

    int coo_index = 0;  // row,col,mat_valの要素指定用のindex

    sparse_make_hamiltonian(mat_dim, tot_site_num, M_H_JsetFile_name,
                            M_H_OutputFile_name, precision, Boundary_Condition,
                            row, col, mat_val, coo_index);

    cout << "/********************************MATRIX INFO WITH COO "
            "FORMAT************************************************/"
         << endl;
    cout << setw(6) << scientific << setprecision(precision) << left << "row"
         << "  " << setw(6) << scientific << setprecision(precision) << left
         << "col"
         << "  " << setw(6) << scientific << setprecision(precision) << left
         << "mat_val" << endl;
    for (int i = 0; i < mat_elements; i++)
    {
        cout << setw(6) << scientific << setprecision(precision) << left
             << row[i] << "  " << setw(6) << scientific
             << setprecision(precision) << left << col[i] << "  " << setw(6)
             << scientific << setprecision(precision) << left << mat_val[i]
             << endl;
    }

    delete[] row;
    delete[] col;
    delete[] mat_val;
}