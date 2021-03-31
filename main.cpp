#include <iostream>
#include <cmath>
#include <vector>
using namespace std;


template <typename T>
class Matrix{
public:
    T **data;
    int n;

    void clear(){
        for (int i = 0; i < this->n; i++)
        {
            delete [] this->data[i];
        }
        delete [] this->data;
    }

    Matrix(int n, int m){
        this->n = n;
        this->m = m;
        data = new T *[n];
        for(int i = 0;i<n;i++){ data[i] = new T[m]; }
    }

    ~Matrix(){
        clear();
    }

    Matrix(const Matrix &other){
        this->n = other.n;
        this->m = other.m;
        this->data = new T *[this->n];
        for(int i = 0;i<this->n;i++){ this->data[i] = new T[this->m]; }
        for(int i = 0;i<this->n;i++){
            for(int j =0;j<this->m;j++){
                this->data[i][j] = other.data[i][j];
            }
        }
    }

    int getN(){return this->n;}
    int getM(){return this->m;}

    template<typename U>
    friend ostream &operator << (ostream &out, Matrix<U> &matrix);
    template<typename U>
    friend istream &operator >> (istream &in, Matrix<U> &matrix);



    Matrix & operator = (const Matrix &other){
        clear();
        this->n = other.n;
        this->m = other.m;
        this->data = new T *[this->n];
        for(int i = 0;i<this->n;i++){ this->data[i] = new T[this->m]; }

        for(int i = 0;i<this->n;i++){
            for(int j = 0;j<this->m;j++){
                this->data[i][j] = other.data[i][j];
            }
        }
        return *this;
    }

    Matrix & operator + (const Matrix &other){
        Matrix *sum = new Matrix(this->n, this->m);
        for(int i = 0;i<this->n;i++){
            for(int j = 0;j<this->m;j++){
                sum->data[i][j] = this->data[i][j] + other.data[i][j];
            }
        }

        return *sum;
    }

    Matrix & operator - (const Matrix &other){
        Matrix *dif = new Matrix(this->n, this->m);
        for(int i = 0;i<this->n;i++){
            for(int j = 0;j<this->m;j++){
                dif->data[i][j] = this->data[i][j] - other.data[i][j];
            }
        }

        return *dif;
    }

    Matrix & operator * (const Matrix &other){
        Matrix *product = new Matrix(this->n, other.m);
        for(int i = 0; i < this->n; i++)
            for(int j = 0; j < other.m; j++)
            {
                product->data[i][j] = 0;
                for(int k = 0; k < this->m; k++){
                    product->data[i][j] += this->data[i][k] * other.data[k][j];
                }
            }

        return *product;
    }


    void transpose(){
        T **copy;
        copy = new T *[this->n];
        for(int i = 0;i<this->n;i++){ copy[i] = new T[this->m]; }
        for(int i = 0;i<this->n;i++){
            for(int j = 0;j<this->m;j++){
                copy[i][j] = this->data[i][j];
            }
        }

        clear();
        int tmp = this->n;
        this->n = this->m;
        this->m = tmp;
        data = new T *[n];
        for(int i = 0;i<n;i++){ data[i] = new T[m]; }

        for(int i = 0;i< this->n;i++){
            for(int j = 0;j<this->m;j++){
                this->data[i][j] = copy[j][i];
            }
        }

    }
private:
    int m;

};

template <typename T>
class SquareMatrix : public Matrix<T> {
public:
    SquareMatrix(int n) : Matrix<T>(n,n) {}

    SquareMatrix(const Matrix<T> & other) : Matrix<T>(other.n,other.n) {
        for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < this->n; ++j) {
                this->data[i][j] = other.data[i][j];
            }
        }
    }
};

template <typename T>
class IdentityMatrix : public SquareMatrix<T> {
public:
    IdentityMatrix(int n) : SquareMatrix<T>(n) {
        for (int i = 0; i < this->n; ++i) {
            for (int j = 0; j < this->n; ++j) {
                if(i == j) this->data[i][j] = 1;
                else this->data[i][j] = 0;
            }
        }
    }


};

template <typename T>
class EliminationMatrix : public IdentityMatrix<T> {
public:
    EliminationMatrix(int i, int j, const SquareMatrix<T> & other) : IdentityMatrix<T>(other.n) {
        i--; j--;
        T x = 0;
        for (int k = 0; k < other.n; ++k) {
            if(k != j) x -= this->data[i][k] * other.data[k][j];
        }
        x = x / other.data[j][j];
        this->data[i][j] = x;
    }

};

template <typename T>
class PermutationMatrix : public IdentityMatrix<T> {
public:
    PermutationMatrix(int i, int j, const SquareMatrix<T> & other) : IdentityMatrix<T>(other.n) {
        i--; j--;
        for (int k = 0; k < other.n; ++k) {
            swap(this->data[i][k], this->data[j][k]);
        }
    }
};

template <typename T>
ostream &operator << (ostream &out, Matrix<T> &matrix){
    for(int i = 0;i<matrix.n;i++){
        for(int j = 0;j<matrix.m - 1;j++){
            if(matrix.data[i][j] == 0) out << 0.00 << " ";
            else out <<  matrix.data[i][j] << " ";
        }
        if(matrix.data[i][matrix.m - 1] == 0) out << 0.00 << endl;
        else out <<  matrix.data[i][matrix.m - 1] << endl;
    }

    return out;
}

template <typename T>
istream &operator >> (istream &in, Matrix<T> &matrix){
    for(int i = 0;i<matrix.n;i++){
        for(int j = 0;j<matrix.m;j++){
            in >> matrix.data[i][j];
        }
    }

    return in;
}

template <typename T>
SquareMatrix<T> & inverse(SquareMatrix<T> & A, int n){
    SquareMatrix<double> A_inverse = IdentityMatrix<double>(n);


    for (int i = 0; i < n-1; ++i) {
        for (int j = i+1; j < n; ++j) {
            if(A.data[j][i] != 0){
                EliminationMatrix<double> E(j+1, i+1, A);
                A = E * A;
                A_inverse = E * A_inverse;
            }
        }
    }

    for (int i = n-1; i > 0; i--) {
        for (int j = i-1; j > -1; j--) {
            if(A.data[j][i] != 0){
                EliminationMatrix<double> E(j+1, i+1, A);
                A = E * A;
                A_inverse = E * A_inverse;
            }
        }
    }

    for (int i = 0; i < n; ++i) {
        double div = A.data[i][i];
        for (int j = 0; j < n; ++j) {
            A_inverse.data[i][j] = A_inverse.data[i][j] / div;
        }
        A.data[i][i] = A.data[i][i] / div;
    }

    A = A_inverse;
    return A;
}

template <typename T>
Matrix<T> & leastSquareApproximation(int n, int m, vector<T> t_i, vector<T> b_i){
    Matrix<T> A(m, n+1);
    Matrix<T> b(m, 1);

    for (int i = 0; i < m; ++i) {
        A.data[i][0] = 1;
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 1; j < n+1; ++j) {
            A.data[i][j] = pow(t_i[i], j);
        }
        b.data[i][0] = b_i[i];
    }

    cout << "A:\n"<<A;

    Matrix<T> A_copy = A;
    A.transpose();

    SquareMatrix<T> ATA_1 = A * A_copy;
    SquareMatrix<T> ATA = ATA_1;
    cout << "A_T*A:\n" << ATA_1;
    inverse(ATA_1, n+1);
    cout << "(A_T*A)^-1:\n" << ATA_1;

    Matrix<T> ATb = A * b;
    cout << "A_T*b:\n" << ATb;

    Matrix<T> *X = &(ATA_1 * ATb);
    cout << "x~:\n" << *X;

    return *X;

}

int main(){
    std::cout.setf(std::ios::fixed);
    std::cout.precision(4);/////

    cout << "Size of dataset:\n";
    int m; cin >> m;
    vector<double> t_i, b_i;
    cout << "Enter data like this : Data1 Data2\n";
    for (int i = 0; i < m; ++i) {
        double input_t, input_b;
        cin >> input_t >> input_b;
        t_i.push_back(input_t);
        b_i.push_back(input_b);
    }
    cout << "polynom degree\n";
    int n; cin >> n;

    Matrix<double> X = leastSquareApproximation(n,m,t_i,b_i);
    cout << "y = ";
    for (int i = 0; i < X.getN(); ++i) {
        if(i>1)
            if(X.data[i][0] > 0) cout << " + "<<X.data[i][0]<<"*x^"<<i;
            else cout << " - "<<abs(X.data[i][0])<<"*x^"<<i;
        else
            if(i == 0) cout << X.data[i][0];
            else  if(X.data[i][0] > 0) cout << " + "<<X.data[i][0]<<"*x";
            else cout << " - "<<abs(X.data[i][0])<<"*x";
    }

}