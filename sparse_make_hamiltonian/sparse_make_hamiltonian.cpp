/*M_H = Make_Hamiltonian*/
#include "../all.h"

using namespace std;

void sparse_make_hamiltonian(int mat_dim, int tot_site_num,
                             std::string M_H_JsetFile_name,
                             std::string M_H_OutputFile_name, int precision,
                             std::string Boundary_Condition, int *row, int *col,
                             double *mat_val, int &coo_index)
{
    int bond_num;
    double szz;

    /*jset.txtからのbondごとの相互作用情報の取得*/
    /*bond数の取得*/
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

    if (Boundary_Condition == "y")
    {
        for (int j = 0; j < mat_dim; j++)
        {
            szz = 0.;
            for (int site_num = 0; site_num < tot_site_num; site_num++)
            {
                // spm(j, site_num, tot_site_num, J, row, col, mat_val,
                // coo_index); smp(j, site_num, tot_site_num, J, row, col,
                // mat_val, coo_index); szz(j, site_num, tot_site_num, J, row,
                // col, mat_val, coo_index);
                spin_operator(j, site_num, tot_site_num, J, row, col, mat_val,
                              coo_index, szz);
            }
            row[coo_index] = j;
            col[coo_index] = j;
            mat_val[coo_index] = szz;
            std::cout << "/****************j=" << j << "**********************/"
                      << endl;
            std::cout << setw(6) << scientific << setprecision(5) << left
                      << "row[" << coo_index << "]"
                      << "  " << setw(6) << scientific << setprecision(5)
                      << left << "col[" << coo_index << "]"
                      << "  " << setw(6) << scientific << setprecision(5)
                      << left << "mat_val[" << coo_index << "]" << endl;
            std::cout << setw(6) << scientific << setprecision(5) << left
                      << row[coo_index] << "  " << setw(6) << scientific
                      << setprecision(5) << left << col[coo_index] << "  "
                      << setw(6) << scientific << setprecision(5) << left
                      << mat_val[coo_index] << endl;
            coo_index++;
        }
    }
    else if (Boundary_Condition == "n")
    {
        for (int j = 0; j < mat_dim; j++)
        {
            szz = 0.;
            for (int site_num = 0; site_num < tot_site_num - 1; site_num++)
            {
                // spm(j, site_num, tot_site_num, J, row, col, mat_val,
                // coo_index); smp(j, site_num, tot_site_num, J, row, col,
                // mat_val, coo_index); szz(j, site_num, tot_site_num, J, row,
                // col, mat_val, coo_index);
                spin_operator(j, site_num, tot_site_num, J, row, col, mat_val,
                              coo_index, szz);
            }
            row[coo_index] = j;
            col[coo_index] = j;
            mat_val[coo_index] = szz;
            coo_index++;
        }
    }
    else
    {
        cout << "ERROR : Maybe inputed other than \"y\" and \"n\" " << endl;
    }

    // /*OUTPUT HAMILTONIAN*/
    ofstream M_H_Output(M_H_OutputFile_name);

    M_H_Output.close();
    delete[] J;
}