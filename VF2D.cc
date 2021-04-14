#include <iostream>
#include <math.h>
#include <ctime>
#include <fstream>
#include <string>
#include<vector>
#include "Dense"
#include "Sparse"

using namespace std;
using namespace Eigen;


int main()
{

int Nx(4);
int Ny(3);
VectorXd X(Nx+1);
VectorXd Y(Ny+1);
double hx(1./Nx);
double hy(1./Ny);

double mu(1.); // a modifier

// Vecteur des points equidistants
X(0)=0;
Y(0)=0;
for (int i=1;i<=Nx;i++)
{
  X(i)=X(i-1)+hx;
}

for (int i=1;i<=Ny;i++)
{
  Y(i)=Y(i-1)+hy;
}


// Vecteurs des demi-points equidistants
VectorXd X_demi(Nx+2);
VectorXd Y_demi(Ny+2);
// On remplit les points milieux
for (int i=1;i<=Nx;i++)
{
  X_demi(i)= (X(i)+X(i-1))/2;
}
for (int i=1;i<=Ny;i++)
{
  Y_demi(i)= (Y(i)+Y(i-1))/2;
}
// points fantômes du bords
X_demi(0) = X(0) - (X_demi(1)-X(0)); // X(-1/2)
X_demi(Nx+1) = X(Nx) + (X(Nx) - X_demi(Nx)); // X(N+1/2)
Y_demi(0) = Y(0) - (Y_demi(1)-Y(0)); // Y(-1/2)
Y_demi(Ny+1) = Y(Ny) + (Y(Ny) - Y_demi(Ny)); // Y(N+1/2)

// Schema numérique
VectorXd T((Nx+1)*(Ny+1));
VectorXd T0((Nx+1)*(Ny+1)); // condition initiale

for (int i=0;i<(Nx+1)*(Ny+1);i++)
{
  T0(i) = 20. ;
}


double tmax(1.), dt(0.001) ;
int nb_iteration(tmax/dt) ;

// initialisation du vecteur à t=0
T=T0;

// cout << "Température initiale =" << endl;
// for (int i=0;i<(Nx+1)*(Ny+1);i++)
// {
//   cout <<
// }


// Construction de la matrice H du laplacien
SparseMatrix<double> H((Nx+1)*(Ny+1),(Nx+1)*(Ny+1));
vector<Triplet<double>> triplets;
for (int i=0; i<=Nx ; ++i)
{
  for (int j=0; j<=Ny ; j++)
  {
    triplets.push_back({j*Nx+i,j*Nx+i,(2./pow(hx,2)+2./pow(hy,2))});
    if (i > 0)
    triplets.push_back({j*Nx+i,j*Nx+i-1,-1./pow(hx,2)});
    if (i < Nx-1)
    triplets.push_back({j*Nx+i,j*Nx+i+1,-1./pow(hx,2)});
    if (j > 0)
    triplets.push_back({j*Nx+i,(j-1)*Nx+i,-1./pow(hy,2)});
    if (j < Ny-1)
    triplets.push_back({j*Nx+i,(j+1)*Nx+i,-1./pow(hy,2)});
  }
}
H.setFromTriplets(triplets.begin(), triplets.end());
cout << "Laplacien =" << H << endl;

for (int n = 0; n < nb_iteration; n++)
  { // Boucle en temps
    T += dt*(-mu*H*T);
  }

  cout << "Température finale =" << T << endl;


}
