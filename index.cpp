#include <bits/stdc++.h>
#include "antoine_solver.h"
#include "linear_system_solver.h"
#include "energy_solver.h"
using namespace std;
double R, D, Cond; //Cond == L4
double P_total = 1;
int no_of_comp = 2;
int no_of_stages = 4;
double T_Feed;
double F1, F2;
vector<vector<double>> H_param(no_of_comp,vector<double>(3));
vector<vector<double>> h_param(no_of_comp,vector<double>(3));
vector<double>A(no_of_comp+1),B(no_of_comp+1),C(no_of_comp+1);
  

int main()
{
  cout<<"Enter Reflux_ratio, Distillate, Condensate, Total_pressure, Temperature_of_Feed, Flow_rate_of_1_in_Feed, Flow_rate_of_2_in_Feed,"<<endl;
  cin>>R>>D>>Cond>>P_total>>T_Feed>>F1>>F2;

  cout<<"Enter the 'H' parameters for each component in order of A, B, C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    H_param[i][0] = A1;
    H_param[i][1] = B1;
    H_param[i][2] = C1;
  }

  cout<<"Enter the 'h' parameters for each component in order of A, B, C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    h_param[i][0] = A1;
    h_param[i][1] = B1;
    h_param[i][2] = C1;
  }

  cout<<"Enter the Antoine parameters for each component in order of A, B and C"<<endl;
  for(size_t i = 1; i<=no_of_comp; i++)
  {
    double A1, B1, C1;
    cin>>A1>>B1>>C1;
    A[i] = A1;
    B[i] = B1;
    C[i] = C1;
  }

  





  
  // L1 = D*R;
  // V2 = D+L1;
  vector<double> T(no_of_stages+1),V(no_of_stages+1),
  L(no_of_stages+1);
  vector<vector<double>>Sij,lij,xij;



  
  while (true)
  {
    //what to initialise here??



   //loop2 needs T,V,L
    while (true)
    {
      // calculation of K (y=K*x);
      vector<vector<double>> K = calculate_K(A,B,C,T, P_total,no_of_stages);

      // calcution of S (S[1][2]=k[1][2]*V[2]/L[2]);
      vector<vector<double>> S = calculate_S(K, V, L,D);

      // make matrix A,X and B for A*X=C
      vector<vector<double>>l=matrix_solver(S,lf);

      // we got X all l
      // now we calculate all xij;
      vector<vector<double>>x=calculate_x(l);

      //  Find the temperature using Newton's method
      vector<double> Tnew= temp_solver(T,A,B,C,x);
      //breaking condition for T
      if(accurate(Tnew,T)){Sij=S,lij=l; T=Tnew;xij=x;break;}
    }
    // Initializing and calculating vij.
    vector<vector<double> > vij = calculate_vij(Sij,lij);
    // Initializing and calculating Yij.
    vector<vector<double> > Yij = calculate_Yij(vij);

    // Assuming that we've taken Enthalpy parameters for every component in 2D vectors named H_param and h_param.

    // Initializing and calculating Hij.

    //make H_param in A,B,C form if possible??
    vector<vector<double> > Hij = calculate_Hij(H_param,T);

    // Initializing and calculating hij.
    vector<vector<double> > hij = calculate_hij(h_param,T);

    // Initializing and calculating Hi.
    vector<double> Hi = calculate_Hi(Hij,Yij);

    // Initializing and calculating hi.
    vector<double> hi = calculate_hi(hij,xij);

    // Calculating hfi
    vector<double> hfi = calculatehfi(h_param, Feed_T);

    // Assuming distillate as 'distillate'
    vector<double> Vnew = calculate_Vnew(Hi, hi, Li, Vi, distillate, F1, F2, hfi);


  }
  return 0;
}