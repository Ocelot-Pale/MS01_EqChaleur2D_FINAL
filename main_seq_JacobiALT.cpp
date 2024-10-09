#include <iostream>
#include <vector>
#include <math.h>
#include <string>
#include <fstream>
#include <stdexcept>

using namespace std;
using Matrix = vector<vector<double>>;
using Vector = vector<double>;
const double PI = 3.14159;

Matrix init_matrix(int l, int c, double initial_value) {
    Matrix result(l, Vector(c, initial_value));
    return result;
}

void display(const Vector& vec) {
    cout << "<------------------------------------[";
    for (double composante : vec) {
        cout << composante << " ";
    }
    cout << "]\n";
}

void f_limite(int num_condition, Vector& vec, double delta, double shift) {
    // Apply boundary condition
    if (num_condition == 0) {
        int k = 0;
        double t = shift;
        for (double& composante : vec) {
            t += delta;
            vec[k] = 0;
            k++;
        }
    }

    if (num_condition == 1) {
        int k = 0;
        double t = shift;
        for (double& composante : vec) {
            t += delta;
            vec[k] = sin(-PI * t / 2);
            k++;
        }
    }

    if (num_condition == 2) {
        int k = 0;
        double t = shift;
        for (double& composante : vec) {
            t += delta;
            vec[k] = sin(PI * t / 2);
            k++;
        }
    }

    if (num_condition == 3) {
        int k = 0;
        double t = shift;
        for (double& composante : vec) {
            t += delta;
            vec[k] = cos(PI * t);
            k++;
        }
    }
}

void f_source(int num_condition, Matrix& mat, double delta, int Nx, int Ny) {
    if (num_condition == 0) {
        for (int ny = 0; ny < Ny; ny++) {
            for (int nx = 0; nx < Nx; nx++) {
                mat[ny][nx] = (-5 * PI * PI / 4) * sin(PI * nx * delta) * cos(PI * ny * delta);
            }
        }
    }
    return;
}

void init_CL(const Matrix& CLs, Matrix& solution, int Nx, int Ny) {
    solution[0] = CLs[0];
    solution[Ny - 1] = CLs[1];

    for (int ny = 0; ny < Ny; ny++) {
        solution[ny][0] = CLs[2][ny];
        solution[ny][Nx - 1] = CLs[3][ny];
    }
}

void write_matrix(const Matrix& m, string path) {
    ofstream file(path);
    if (!file.is_open()) {
        cerr << ("LE FICHIER " + path + " n'a pas pu être ouvert !");
        return;
    }
    int n_line = m.size();
    int n_row = m[0].size();
    for (int l = 0; l < n_line; l++) {
        for (int c = 0; c < n_row; c++) {
            file << m[l][c];
            if (c != n_row - 1) {
                file << "; ";
            }
        }
        if (l != n_line - 1) {
            file << endl;
        }
    }
    file.close();
}

void resolution_jacobi_sequentielle(const Matrix& T_SRC, Matrix& solution,
                                    int Lx, int Ly,
                                    int Nx, int Ny,
                                    double delta, int N_ITER) {
    Matrix sol_temp = solution;
    for (int iter = 0; iter < N_ITER; iter++) {
        for (int nx = 1; nx < Nx - 1; nx++) {
            for (int ny = 1; ny < Ny - 1; ny++) {
                sol_temp[ny][nx] = ((solution[ny][nx + 1] + sol_temp[ny][nx - 1]) +
                                    (solution[ny + 1][nx] + sol_temp[ny - 1][nx]) -
                                    delta * delta * T_SRC[ny][nx]) / 4;
            }
        }
        swap(sol_temp, solution);
    }
}

int main(int argc, char** argv) {
    double Lx = 1.0;
    double Ly = 1.0;
    double delta = 0.01;
    int N_ITER = 50000;
    int NxGlobal = static_cast<int>(Lx / delta);
    int NyGlobal = static_cast<int>(Ly / delta);

    Matrix CL_MATRIX;

    Vector CL_gauche(NyGlobal, 0);
    f_limite(0, CL_gauche, delta, 0);

    Vector CL_droite(NyGlobal, 0);
    f_limite(3, CL_droite, delta, 0);

    Vector CL_haut(NxGlobal, 0);
    f_limite(1, CL_haut, delta, 0);

    Vector CL_bas(NxGlobal, 0);
    f_limite(2, CL_bas, delta, 0);

    CL_MATRIX = {CL_bas, CL_haut, CL_gauche, CL_droite};

    Matrix Terme_Source = init_matrix(NyGlobal, NxGlobal, 0);
    f_source(0, Terme_Source, delta, NxGlobal, NyGlobal);
    
    Matrix Solution = init_matrix(NyGlobal, NxGlobal, 0);
    
    init_CL(CL_MATRIX, Solution, NxGlobal, NyGlobal);

    resolution_jacobi_sequentielle(Terme_Source, Solution, Lx, Ly, NxGlobal, NyGlobal, delta, N_ITER);

    write_matrix(Solution, "solution_sequentielle.txt");

    return 0;
}
