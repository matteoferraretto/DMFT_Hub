#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "FermionicFockState.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>

/* Eigen names */
using Eigen::MatrixXcf;
using Eigen::MatrixXf;
using Eigen::SelfAdjointEigenSolver;
using Eigen::VectorXf;
using Eigen::VectorXcf;
using Eigen::VectorXi;
using Eigen::SparseMatrix;
/* std names */
using std::cout;
using std::endl;
using namespace std::literals;
/* chrono names */
using std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int BinarySearch(Eigen::RowVectorXi, int);

int main(){
    int dim = 1000;
    MatrixXcf H(dim, dim);
    H = MatrixXcf::Zero(dim, dim);
    for(int i=0; i<dim; i++){
        H(i, i) = 1.0f + 0.0if;
        if(i != dim-1) { H(i, i+1) = 1.0f + 2.0if; }
        if(i != 0) { H(i, i-1) = 1.0f - 2.0if; }
    }
    
    SparseMatrix<std::complex<float>> H_sparse(dim, dim);
    H_sparse.reserve(VectorXi::Constant(dim, 4)); // 4: estimated number of non-zero enties per column
    for (int i = 0; i < dim; i++){
        H_sparse.coeffRef(i, i) = 1.0f;
        if(i != dim-1) { H_sparse.coeffRef(i, i+1) = 1.0f + 2.0if; }
        if(i != 0) { H_sparse.coeffRef(i, i-1) = 1.0f - 2.0if; }
    }
    H_sparse.makeCompressed();
  
    //      DENSE DIAGONALIZATION
    // start clock
    auto start_time = steady_clock::now();

    SelfAdjointEigenSolver<MatrixXcf> eigensolver;
    //eigensolver.computeFromTridiagonal(v, w);
    eigensolver.compute(H);
    VectorXf evals = eigensolver.eigenvalues();
    MatrixXcf evecs = eigensolver.eigenvectors();
    //MatrixXf diag = evals.asDiagonal();

    // end clock
    auto end_time = steady_clock::now();
/*    
    cout << "Eigenvalues of H = \n" << evals << endl;
    cout << "Eigenvectors of H = \n" << evecs << endl;
    cout << "Consistency check: \n" << (evecs.adjoint())*diag*evecs << endl;
*/
    cout << "matrix dimension: " << dim << "\n";
    cout << "elapsed time: " << duration_cast<milliseconds>(end_time - start_time).count() << " milliseconds.\n";

    //      SPARSE DIAGONALIZATION
    // start clock
    start_time = steady_clock::now();

    SelfAdjointEigenSolver<SparseMatrix<std::complex<float>>> eigensolver_sparse;
    //eigensolver.computeFromTridiagonal(v, w);
    eigensolver_sparse.compute(H_sparse);
    evals = eigensolver_sparse.eigenvalues();
    evecs = eigensolver_sparse.eigenvectors();
    //MatrixXf diag = evals.asDiagonal();

    // end clock
    end_time = steady_clock::now();

    cout << "matrix dimension: " << dim << "\n";
    cout << "elapsed time: " << duration_cast<milliseconds>(end_time - start_time).count() << " milliseconds.\n";

/*
    Eigen::RowVectorXi vec(8);
    vec << 10,20,30,40,50,60,70,80;
    cout << vec.cols() << endl;
    for (int i=0; i<8; i++){
        cout << "value = " << vec(i) << "; index = " << BinarySearch(vec, vec(i)) << endl;
    }
*/
    std::vector<FermionicFockState> states;
    states.reserve(2);
    FermionicFockState psi(3,2,2);
    for(int i=0; i<2; i++){
        states.push_back(psi);
    }

    return 0;
}


int BinarySearch(Eigen::RowVectorXi vec, int element){
        // assuming that vec is sorted in ascending order
        int dim = vec.cols();
        int index = (dim-1)/2; // choose index in the middle of the list 
        int value = vec(index);
        while(value != element && (index >= 0 || index < dim)){
            if (value > element){
                dim = index; 
                index = index/2; 
            }
            else {
                index += dim/2; 
            }
            value = vec(index);
        }
        return index;
    }
