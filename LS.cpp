#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

void print(vector<vector<double> > A) {
    int i, j;
    for (i = 0 ; i < A.size() ; i++) {
        for (j = 0 ; j < A[1].size() ; j++) {
            cerr << setw(12) <<  A[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << endl;
}

void print(vector<double> A) {
    int i;
    for (i = 0 ; i < A.size() ; i++) {
        cerr << setw(12) <<  A[i] << endl;
    }
    cerr << endl;
}

void lsnormal (vector<vector<double> > A, vector<double> B, vector<double> *X) {
    
    vector<vector<double> > C, G;
    vector<double> D(A[1].size(),0);
    C.resize(A[1].size());
    G.resize(A[1].size());
    double sum;
    // Compute lower triangular portion of (A^T)A
    int i,j,k;
    for (i = 0 ; i < A[1].size() ; i++) {
        C[i].resize(A[1].size());
        for (j = 0 ; j <= i ; j++) {
            for (k = 0 ; k < A.size() ; k++) {
                C[i][j]=C[i][j]+A[k][i]*A[k][j];
                if (k == A.size() - 1) (C[j][i]=C[i][j]);
            }
        }
    }
    
    // Form (A^T)b
    for (i = 0 ; i < A[1].size() ; i++) {
        for (j = 0 ; j < A.size() ; j++) {
            D[i] = D[i]+A[j][i]*B[j];
        }
    }
    
    // Cholesky Factorization
    for(i = 0 ; i < C.size() ; i++){
        G[i].resize(C[1].size());
        for(j = 0 ; j <= i ; j++){
            sum=0;
            for(k = 0 ; k < j ; k++) sum += G[i][k]*G[j][k];
            G[i][j]= (i == j) ? (sqrt(C[i][j]-sum)) : ((C[i][j]-sum)/G[j][j]);
        }
    }
    
    // Solve System
    D[0] = D[0]/G[0][0];
    for (i = 1 ; i < D.size() ; i++) {
        sum = 0;
        for (j = 0 ; j < i ; j++) {
            sum += G[i][j]*D[j];
        }
        D[i] = (D[i] - sum)/G[i][i];
    }
    
    D[D.size()-1] = D[D.size()-1]/G[D.size()-1][D.size()-1];
    for (i = D.size() - 2 ; i > -1 ; i--) {
        sum = 0;
        for (j = i + 1 ; j < D.size() ; j++) {
            sum += G[j][i]*D[j];
        }
        D[i] = (D[i] - sum)/G[i][i];
    }
    *X = D;
}

void lsmgs (vector<vector<double> > A, vector<double> B, vector<double> *Y) {
    
    vector<vector<double> > R, Q;
    R.resize(A[1].size());
    Q.resize(A.size());
    vector<double> D(A[1].size(),0);
    int i,j,k;
    double sum;
    
    for (i = 0 ; i < A.size() ; i++) {
        R[i].resize(A[1].size());
        Q[i].resize(A[1].size());
    }
    
    // Modified Gram-Schmidt Algorithm
    for (i = 0 ; i < A[1].size() ; i++) {
        for (j = 0 ; j < A.size() ; j++) {
            R[i][i] = R[i][i] + A[j][i]*A[j][i];
        }
        R[i][i] = sqrt(R[i][i]);
        for (j = 0 ; j < A.size() ; j++) {
            Q[j][i] = A[j][i]/R[i][i];
        }
        
        for (j = i + 1 ; j < A[1].size() ; j++) {
            for (k = 0 ; k < A.size() ; k++) {
                R[i][j] = R[i][j] + Q[k][i] * A[k][j];
            }
            for (k = 0 ; k < A.size() ; k++) {
                A[k][j] = A[k][j] - Q[k][i] * R[i][j];
            }
        }
    }
  
    // Form (Q^T)b
    for (i = 0 ; i < A[1].size() ; i++) {
        for (j = 0 ; j < A.size() ; j++) {
            D[i] = D[i]+Q[j][i]*B[j];
        }
    }
    
    //Solve System
    D[D.size()-1] = D[D.size()-1]/R[D.size()-1][D.size()-1];
    for (i = D.size() - 2 ; i > -1 ; i--) {
        sum = 0;
        for (j = i + 1 ; j < D.size() ; j++) {
            sum += R[i][j]*D[j];
        }
        D[i] = (D[i] - sum)/R[i][i];
    }
    *Y = D;
}
    
int main (int argc, char** argv) {
    
    int m,n,temp_m,i,j;
    double e1 = 0, e2 = 0;
    vector<vector<double> > A;
    vector<double> B, X, Y, sol(6,.4082);
    
    // reads matrix and vector from text file
    ifstream matrix, vector;
    matrix.open("A.txt");
    vector.open("B.txt");
    if (!matrix || !vector) cerr << "\nError opening file." << endl;
        matrix >> m >> n;
        vector >> temp_m;
    if (m < n || (m != temp_m)) {
        cerr << "Dimension error. " << endl;
        exit(1);
    } else {
        A.resize(m);
        B.resize(m);
        for (i = 0 ; i < m ; i++) {
            A[i].resize(n);
            vector >> B[i];
            for (j = 0 ; j < n ; j++) {
                matrix >> (A[i])[j];
            }
        }
    }
    matrix.close();
    vector.close();
    cerr << "Matrix A: " << endl;
    print(A);
    cerr << "Vector B: " << endl;
    print(B);
    
    lsnormal(A,B,&X);
    cerr << "\nNormal Equations Solution" << endl;
    print(X);
    lsmgs(A,B,&Y);
    cerr << "\nModified Gram Schmidt Solution" << endl;
    print(Y);
    
    for (i = 0 ; i < X.size() ; i++) {
        e1 += pow(X[i] - sol[i],2);
        e2 += pow(Y[i] - sol[i],2);
    }
    e1 = sqrt(e1);
    e2 = sqrt(e2);
    cerr << "\ne1 = " << e1 << endl;
    cerr << "\ne2 = " << e2 << endl;
    cerr << "\nThe difference vector: " << endl;
    for (i = 0 ; i < X.size() ; i++) {
        cerr << X[i] - Y[i] << endl;
    }
    
    return 0;
}
