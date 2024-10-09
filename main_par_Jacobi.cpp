#include <iostream> 
#include <vector> 
#include <math.h> 
#include <string> 
#include <fstream> 
#include <mpi.h> 
#include <stdexcept>
using namespace std ;
 using Matrix = vector<vector<double>>;
 using Vector = vector<double>;  
const double PI = 3.14159 ; 


Matrix init_matrix(int l , int c , double initial_value){
    Matrix result(l , Vector(c,initial_value));
    return  result;
}

void display(const Vector &  vec){
    int myRank ; 
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    cout << myRank<<"<------------------------------------[";
    for (double composante : vec){
        cout << composante << " ";
    }
    cout << "]\n";
}

void f_limite(int num_condition , Vector & vec , double delta ,double shift){
    //0 nulle 
    //1 -sin pi/2 t
    //2 +sin pi/2 t
    //3 +cos pi t 
    if (num_condition == 0){
        int k = 0;
        double t = shift;
        for (double composante : vec){
           t += delta;
            vec[k] = 0;
            k++;
            }
            
    } 

    if (num_condition == 1){
        int k = 0;
        double t = shift;
        for (double composante : vec){
           t += delta;
            vec[k] = sin(-PI *t/2);
            k++;
            }
            
    } 

    if (num_condition == 2){
        int k = 0;
        double t = shift;
        for (double composante : vec){
           t += delta;
            vec[k] = sin(PI *t/2);
            k++;
            }
            
    } 

    if (num_condition == 3){
        int k = 0;
        double t = shift; 
        for (double composante : vec){ 
           t += delta;
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
                    //cout << (-5*PI*PI/4) * sin(PI*nx *delta) * cos(PI*ny*delta) << endl; 
                }
            }
        }
    return ; 
}

void init_CL( const Matrix &  CLs , Matrix & solution , int Nx , int Ny , int rank , int nbTask ){
    //cout << "========== IMPORT DES CONDITION LIMITES PROCESSUS : " << rank << endl ; 
    solution[0] = CLs[0] ; 
    solution[Ny-1] = CLs[1] ;

    if (rank == 0 ){    
        for (int ny = 0 ; ny < Ny ; ny++){
        solution[ny][0] = CLs[2][ny]; 
        }}

    if (rank == nbTask- 1 ){
        for (int ny = 0 ; ny < Ny ; ny++){ 
            solution[ny][Nx-1] = CLs[3][ny];}   
    }
}


void write_matrix(const Matrix & m ,string path){
    
    ofstream file(path) ; 
    if (!file.is_open()){cerr << ("LE FICHIER " + path + " n'a pas pu être ouvert !");}
    int n_line = m.size() ; 
    int n_row = m[0].size();
    for (int l = 0 ; l<n_line ; l++){
        for(int c = 0 ; c<n_row ; c++){
            file << m[l][c] ;
            if (c != n_row-1){file<< "; ";} 
        }
        if (l!=n_line-1){file << endl ;} 
    }
    file.close() ;
    return ;
}
void resolution_jacobi_sequentielle( const Matrix &  T_SRC , Matrix &  solution ,
                                        int Lx , int Ly ,
                                        int Nx , int Ny ,
                                        double delta , int N_ITER ){
    Matrix sol_temp = solution;  
    for (int iter = 0; iter < N_ITER ; iter ++ ){
        for (int nx = 1 ; nx < Nx-1 ; nx ++){
            for(int ny = 1 ; ny<Ny-1 ; ny++){
                sol_temp[ny][nx] = ((solution[ny][nx+1] +sol_temp[ny][nx-1] ) + (solution[ny+1][nx] + sol_temp[ny-1][nx])  
                - delta * delta * T_SRC[ny][nx])/4 ; 
            }
        }
        swap(sol_temp , solution) ; 
    }
};

//RESOLUTION PARALLE PARTITION 1D
void resolution_jacobi_parallele( const Matrix &  T_SRC , Matrix &  solution ,
                                        int Lx , int Ly ,
                                        int Nx , int Ny ,
                                        double delta , int N_ITER,
                                        int rank, int maxRank){
    
    Matrix sol_temp = solution;
    MPI_Barrier(MPI_COMM_WORLD);
    //cout << "[PROCESUS " << rank << " ] début de résolution..." << endl ;   
    Vector CL_droite_to_recieve(Ny , 0); 
    Vector CL_droite_to_send(Ny , 1); 
    Vector CL_gauche_to_recieve(Ny , 0);
    Vector CL_gauche_to_send(Ny , 0);
    double t_inter_2=0;
    double t_inter_1=0;
    for (int iter = 0; iter < N_ITER ; iter ++ ){
        //MPI_Barrier(MPI_COMM_WORLD);
        if (rank == 0 && (iter%100==0)){
            t_inter_1 = MPI_Wtime() ; 
            //printf("Temps intermediaire total: %1.2fs\n", t_inter_1-t_inter_2);
            t_inter_2 = MPI_Wtime() ;
            //cout << "==== ITERATION "<< iter << " Résolution : " << static_cast<double>(iter)/N_ITER*100<< "\%" << endl <<endl;
             }



    if (rank != maxRank - 1) {
        MPI_Sendrecv(&CL_droite_to_send[0], Ny, MPI_DOUBLE, rank + 1, rank * 1000,
                     &CL_droite_to_recieve[0], Ny, MPI_DOUBLE, rank + 1, (rank + 1) * 100,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank != 0) {
        MPI_Sendrecv(&CL_gauche_to_send[0], Ny, MPI_DOUBLE, rank - 1, rank * 100,
                     &CL_gauche_to_recieve[0], Ny, MPI_DOUBLE, rank - 1, (rank - 1) * 1000,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }


//MPI_Barrier(MPI_COMM_WORLD);
// UPDATE DES CL
        
        for (int ny = 0 ; ny < Ny ; ny++){
            if (rank != 0){
                solution[ny][0] = CL_gauche_to_recieve[ny] ;}
            if (rank != maxRank-1){solution[ny][Nx-1] = CL_droite_to_recieve[ny];}
        } 

        for (int nx = 1 ; nx < Nx-1 ; nx ++){
            for(int ny = 1 ; ny<Ny-1 ; ny++){
                sol_temp[ny][nx] = ((solution[ny][nx+1] +solution[ny][nx-1] ) + (solution[ny+1][nx] + solution[ny-1][nx])  
                - delta * delta * T_SRC[ny][nx])/4 ; 
            }
        }
        swap(sol_temp , solution) ; 

        for (int ny = 0 ; ny < Ny ; ny++){
            if (rank != 0){
                CL_gauche_to_send[ny] = solution[ny][1]  ;
                }
            if (rank != maxRank - 1 ){
                CL_droite_to_send[ny] = solution[ny][Nx-2];
                } 
        } 
        //cout << "ITER:" << iter <<"[PROCESSUS " << rank << " Nx=" << Nx << "] TO SEND GAUCHE:" <<endl ; 
        //display(CL_gauche_to_send) ;
        //write_matrix(solution , "debug/intermediaire_iter" + to_string(iter+0.5) + "_tache" +  to_string(rank) );
    }
    //cout <<  "================= [PROCESUS " << rank << " ] FIN DE RESOLUTION" <<endl ; 
    return ; 
};


//=============================//
int main(int argc , char** argv){
    double t0 = MPI_Wtime(); 
    MPI_Init(&argc, &argv); 
    int nbTask ; 
    int myRank ; 
    MPI_Comm_size(MPI_COMM_WORLD, &nbTask);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    //cout<<"Processus " << myRank << "/" << nbTask-1 << " lancé" << endl ; 
    double Lx = 1. ; 
    double Ly = 1. ;
    double delta = 0.005;
    int N_ITER =10000;
    int NxGlobal = static_cast<int>(Lx/delta) ; 
    int NyGlobal = static_cast<int>(Ly/delta) ;

    int NxLocal = NxGlobal/nbTask ; 
    if (myRank==nbTask - 1){NxLocal +=  NxGlobal%nbTask;}
    Matrix CL_MATRIX ; 
    //cout << "Processus " << myRank << ": " << NxLocal << "/" << NxGlobal << " colonnes à charge" << endl ; 
    //Matrix Solution;
    double xShift = NxLocal * myRank * delta ; 

    if ( myRank == nbTask -1){
        NxLocal += 2;
        xShift = Lx - NxLocal * delta; 

        Vector CL_gauche(NyGlobal , 0) ; 
        f_limite(0, CL_gauche , delta,0 ); 

        Vector CL_droite(NyGlobal , 0) ; 
        f_limite(3, CL_droite , delta,0);


        Vector CL_haut(NxLocal, 0) ;
        f_limite(1, CL_haut , delta , xShift ); 

        Vector CL_bas(NxLocal , 0) ;
        f_limite(2 , CL_bas,delta , xShift ); 

        CL_MATRIX = {CL_bas,CL_haut,CL_gauche,CL_droite}; 
    }
     
    else {if (myRank==0){
        NxLocal += 2;
        Vector CL_gauche(NyGlobal , 3) ; 
        f_limite(0, CL_gauche , delta ,0); 

        Vector CL_droite(NyGlobal , 0) ; 
        f_limite(3 , CL_droite , delta , 0 );

        Vector CL_haut(NxLocal, 0) ;
        f_limite(1, CL_haut , delta , xShift ); 

        Vector CL_bas(NxLocal , 0) ;
        f_limite(2 , CL_bas , delta,  xShift ); 

        CL_MATRIX = {CL_bas,CL_haut,CL_gauche,CL_droite} ; 
        
    }
    else{
        NxLocal += 4; 

        Vector CL_gauche(NyGlobal , 0) ; 
        f_limite(0, CL_gauche , delta,0 ); 

        Vector CL_droite(NyGlobal , 0) ; 
        f_limite(3, CL_droite , delta , 0 );

        Vector CL_haut(NxLocal, 0) ;
        f_limite(1, CL_haut , delta , xShift ); 

        Vector CL_bas(NxLocal , 0) ;
        f_limite(2 , CL_bas ,delta, xShift ); 

        CL_MATRIX = {CL_bas,CL_haut,CL_gauche,CL_droite} ; 
    }}
    
    Matrix Terme_Source = init_matrix(NyGlobal,NxLocal,0);
    f_source(0,Terme_Source,delta,NxLocal,NyGlobal) ; 
    Matrix Solution = init_matrix(NyGlobal,NxLocal,0) ;
    string path_sol =  "debug/sollution_avant_tout_task" + to_string(myRank) +".txt" ;     
    //write_matrix(Solution , path_sol ) ; 
    
    init_CL(CL_MATRIX , Solution , NxLocal , NyGlobal, myRank,nbTask);
    //path_sol =  "debug/sollution_CL" + to_string(myRank) +".txt" ; 
    //write_matrix(Solution ,path_sol) ;  
    resolution_jacobi_parallele(Terme_Source , Solution , Lx , Ly , NxLocal , NyGlobal , delta , N_ITER, myRank , nbTask); 

    Matrix Solution_Finale ; 
    Solution_Finale = init_matrix(NyGlobal, NxGlobal,0) ; 
     
    write_matrix(Solution , "sorties/" + to_string(nbTask)+ "process_solution_sortie" + to_string(myRank) +".txt") ; 
/*    write_matrix(Terme_Source , "Source.txt");
    */
   double t_fin = MPI_Wtime() ; 
   if (myRank ==0){printf("Temps exécution total: %1.2f\n", t_fin-t0);}
    MPI_Finalize();

    return 0;
}

