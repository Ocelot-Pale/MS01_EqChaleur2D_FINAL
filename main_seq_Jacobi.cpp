#include <iostream> 
#include <vector> 
#include <math.h> 
#include <string> 
#include <fstream> 
#include <ctime>

using namespace std ;
 using Matrix = vector<vector<double>>;
 using Vector = vector<double>;  
const double PI = 3.14159 ; 

Matrix init_matrix(int l , int c , double initial_value){
    Matrix result(l , Vector(c,initial_value));
    return  result;
}

void display(const Vector &  vec){
    cout << "[";
    for (double composante : vec){
        cout << composante << " ";
    }
    cout << "]\n";
}

void f_limite(int num_condition , Vector & vec , double delta ){
    //0 nulle 
    //1 -sin pi/2 t
    //2 +sin pi/2 t
    //3 +cos pi t 
    if (num_condition == 0){
        int k = 0;
        double t = 0 ; ; 
        for (double composante : vec){
           t = k*delta;
            vec[k] = 0;
            k++;
            }
            
    } 

    if (num_condition == 1){
        int k = 0;
        double t = 0 ; ; 
        for (double composante : vec){
           t = k*delta;
            vec[k] = sin(-PI *t/2);
            k++;
            }
            
    } 

    if (num_condition == 2){
        int k = 0;
        double t = 0 ; ; 
        for (double composante : vec){
           t = k*delta;
            vec[k] = sin(PI *t/2);
            k++;
            }
            
    } 

    if (num_condition == 3){
        int k = 0;
        double t = 0 ; ; 
        for (double composante : vec){
           t = k*delta;
            vec[k] = cos(PI *t);
            k++;
            }
            
    } 
}

void f_source(int num_condition , Matrix & mat , double delta , int Nx , int Ny){
        if (num_condition == 0){
            for (int ny = 0 ; ny < Ny ; ny++ ){
                for (int nx = 0  ; nx < Nx ; nx++){
                    mat[ny][nx] = (-5*PI*PI/4) * sin(PI*nx *delta) * cos(PI*ny*delta); 
                }
            }
        }
    return ; 
}

void init_CL( const Matrix &  CLs , Matrix & solution , int Nx , int Ny){
    solution[0] = CLs[0] ; 
    solution[Ny-1] = CLs[1] ;

    for (int ny = 0 ; ny < Ny ; ny++){
        solution[ny][0] = CLs[2][ny]; 
        solution[ny][Ny-1] = CLs[3][ny];
    } 

}


void resolution_jacobi_sequentielle( const Matrix &  T_SRC , Matrix &  solution ,
                                        int Lx , int Ly ,
                                        int Nx , int Ny ,
                                        double delta , int N_ITER ){
    Matrix sol_temp = solution;  
    for (int iter = 0; iter < N_ITER ; iter ++ ){
        //cout << iter << endl ; 
        for (int nx = 1 ; nx < Nx-1 ; nx ++){
            for(int ny = 1 ; ny<Ny-1 ; ny++){
                sol_temp[ny][nx] = ((sol_temp[ny][nx+1] +sol_temp[ny][nx-1] ) + (sol_temp[ny+1][nx] + sol_temp[ny-1][nx])  
                - delta * delta * T_SRC[ny][nx])/4 ; 
            }
        }
        swap(sol_temp , solution) ; 
    }
};

void write_matrix(const Matrix & m ,string path){
    ofstream file(path) ; 
    int n_line = m.size() ; 
    int n_row = m[0].size();
    for (int l = 0 ; l<n_line ; l++){
        for(int c = 0 ; c<n_row ; c++){
            file << m[l][c] << "; "; 
        }
        if (l!=n_line-1){file << endl ;} 
    }
}

//=============================//
int main(int argc , char** argv){
      clock_t start = clock();
    double Lx = 1. ; 
    double Ly = 1. ;
    double delta = 0.05;
    int N_ITER =50;
    int Nx = static_cast<int>(Lx/delta) ; 
    int     Ny = static_cast<int>(Ly/delta) ;


    Vector CL_droite(Ny , 3) ; 
    f_limite(3 , CL_droite , delta ); 

    Vector CL_gauche(Ny , 0) ;
    f_limite(0 , CL_gauche , delta ); 

    Vector CL_haut(Nx, 1) ;
    f_limite(1, CL_haut , delta ); 

    Vector CL_bas(Nx , 2) ;
    f_limite(2 , CL_bas , delta ); 

    Matrix CL_MATRIX = {CL_bas,CL_haut,CL_gauche,CL_droite} ; 

    Matrix Terme_Source = init_matrix(Ny,Nx,0);
    f_source(0,Terme_Source,delta,Nx,Ny) ; 

    Matrix Solution = init_matrix(Ny,Nx,1);
    write_matrix(Solution , "sollution_avant_tout.txt") ; 
    

    init_CL(CL_MATRIX , Solution , Nx , Ny) ; 
    write_matrix(Solution , "sollution_CL.txt") ; 
    
    resolution_jacobi_sequentielle(Terme_Source , Solution , Lx , Ly , Nx , Ny , delta , N_ITER ); 
    write_matrix(Solution , "solution_sortie.txt") ; 
    write_matrix(Terme_Source , "Source.txt");
   clock_t end = clock();    // Fin du chrono
    cout << "Execution time: " << double(end - start) / CLOCKS_PER_SEC << " seconds\n";  
    return 0; 
}

