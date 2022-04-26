/*M_H = Make_Hamiltonian*/
#include "../all.h"

using namespace std;

void sparse_make_hamiltonian(int mat_dim, int tot_site_num,
                             std::string M_H_OutputFile_name, int precision,
                             std::string Boundary_Condition, double *J,
                             int *row, int *col, double *mat_val,
                             int &coo_index)
{
    int bond_num;
    double szz;

    if (Boundary_Condition == "y")
    {
        for (int j = 0; j < mat_dim; j++)
        {
            szz = 0.;
            for (int site_num = 0; site_num < tot_site_num; site_num++)
            {
                spin_operator(j, site_num, tot_site_num, J, row, col, mat_val,
                              coo_index, szz);
            }
            row[coo_index] = j;
            col[coo_index] = j;
            mat_val[coo_index] = szz;
            // // std::cout << "/****************j=" << j <<
            // // "**********************/"
            // //           << endl;
            // // std::cout << setw(6) << scientific << setprecision(5) << left
            // //           << "row[" << coo_index << "]"
            // //           << "  " << setw(6) << scientific << setprecision(5)
            // //           << left << "col[" << coo_index << "]"
            // //           << "  " << setw(6) << scientific << setprecision(5)
            // //           << left << "mat_val[" << coo_index << "]" << endl;
            // // std::cout << setw(6) << scientific << setprecision(5) << left
            // //           << row[coo_index] << "  " << setw(6) << scientific
            // //           << setprecision(5) << left << col[coo_index] << "  "
            // //           << setw(6) << scientific << setprecision(5) << left
            // //           << mat_val[coo_index] << endl;
            if (szz != 0.0) coo_index++;
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
            // std::cout << "/****************j=" << j <<
            // "**********************/"
            //           << endl;
            // std::cout << setw(6) << scientific << setprecision(5) << left
            //           << "row[" << coo_index << "]"
            //           << "  " << setw(6) << scientific << setprecision(5)
            //           << left << "col[" << coo_index << "]"
            //           << "  " << setw(6) << scientific << setprecision(5)
            //           << left << "mat_val[" << coo_index << "]" << endl;
            // std::cout << setw(6) << scientific << setprecision(5) << left
            //           << row[coo_index] << "  " << setw(6) << scientific
            //           << setprecision(5) << left << col[coo_index] << "  "
            //           << setw(6) << scientific << setprecision(5) << left
            //           << mat_val[coo_index] << endl;
            coo_index++;
            if (szz != 0.0) coo_index++;
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