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
        S_L_OutputFile_name;

    int tot_site_num, precision;
    get_data(if_M_H_Settingfile, tot_site_num, M_H_OutputFile_name,
             M_H_JsetFile_name, S_L_OutputFile_name, Boundary_Condition,
             precision);

    std::cout << "/************************************************************"
                 "***************************"
              << "INPUT DATA"
              << "*************************************************************"
                 "**************************/\n";
    std::cout << "tot_site_num        = " << tot_site_num << endl;
    std::cout << "M_H_OutputFile_name = " << M_H_OutputFile_name << endl;
    std::cout << "M_H_JsetFile_name   = " << M_H_JsetFile_name << endl;
    std::cout << "D_L_OutputFile_name = " << S_L_OutputFile_name << endl;
    std::cout << "Boundary Condition  = " << Boundary_Condition << endl;
    std::cout << "precision           = " << precision << endl;
    std::cout << "/************************************************************"
                 "***************************"
              << "*************************************************************"
                 "************************************************/\n";

    /*******************************************************************************************************************************/
    int mat_dim = 1 << tot_site_num;

    /*************************Setting
     * J***************************************************************************/
    int bond_num;
    ifstream if_M_H_JsetFile(M_H_JsetFile_name);
    if (!(if_M_H_JsetFile))
    {
        cerr << "Could not open the file - '" << M_H_JsetFile_name << "'"
             << endl;
    }
    if (Boundary_Condition == "y")
    {
        bond_num = tot_site_num;
    }
    else
    {
        bond_num = tot_site_num - 1;
    }

    double *J = new double[bond_num];
    cout << "/************************************Jset*************************"
            "*****/"
         << endl;
    std::cout << "i"
              << "  "
              << "i+1"
              << ":  "
              << " J[i]      " << endl;
    for (int i = 0; i < bond_num; i++)
    {
        J[i] = 0.;
        if_M_H_JsetFile >> J[i];
        std::cout << i << "   " << i + 1 << "  :  " << J[i] << endl;
    }
    cout << "/*****************************************************************"
            "**"
            "*****/"
         << endl;
    if_M_H_JsetFile.close();
    /****************************************************************************************/
    int mat_elements = 0;
    mat_elements = sparse_count_mat_elements(
        mat_dim, tot_site_num, M_H_JsetFile_name, Boundary_Condition, J);
    cout << "mat_elements = " << mat_elements << endl;
    /*Hamiltonian with COO format*/
    int *row = new int[mat_elements];
    int *col = new int[mat_elements];
    double *mat_val = new double[mat_elements];
    vec_init(mat_elements, row);
    vec_init(mat_elements, col);
    vec_init(mat_elements, mat_val);

    /*Arrays of eigen value and vector*/
    double *eigen_value = new double[mat_dim];
    double *eigen_vec = new double[mat_dim];
    vec_init(mat_dim, eigen_value);
    vec_init(mat_dim, eigen_vec);

    int coo_index = 0;  // row,col,mat_valの要素指定用のindex

    sparse_make_hamiltonian(mat_dim, tot_site_num, M_H_OutputFile_name,
                            precision, Boundary_Condition, J, row, col, mat_val,
                            coo_index);

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

    sparse_lanczos(mat_dim, mat_elements, row, col, mat_val, eigen_value,
                   eigen_vec, S_L_OutputFile_name);

    delete[] row;
    delete[] col;
    delete[] mat_val;
    delete[] eigen_value;
    delete[] eigen_vec;
}