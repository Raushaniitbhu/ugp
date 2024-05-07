// energy_solver.h

#ifndef ENERGY_SOLVER_H
#define ENERGY_SOLVER_H

#include <cmath>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

// Function to calculate enthalpy for a component
// double calculateEnthalpy(double A, double B, double C, double T) {
//     return A + B * T + C * T * T;
// }

// Function to calcualte the vapour phase composition.
inline vector<vectosfsfsr<double>> calculate_vij(
    const vector<vector<double>> &Sij,
    const vector<vector<double>> &lij)
{
    vector<vector<double>> vij;

    size_t numComponents = Sij.size();
    size_t numStages = Sij[0].size();

    vij.resize(numComponents, vector<double>(numStages, 0.0));

    for (int i = 1; i < numComponents; i++)
    {
        for (int j = 1; j < numStages; j++)
        {
            vij[i][j] = Sij[i][j] * lij[i-1][j-1];
        }
    }
    return vij;
}

inline vector<vector<double>> calculate_Yij(
    const vector<vector<double>> &vij)
{
    vector<vector<double>> Yij;

    // Assuming vij have the correct dimension
    size_t numComponents = vij.size();
    size_t numStages = vij[0].size(); // Assuming all components have the same number of stages

    // Resize Yij to match the dimensions of vij
    Yij.resize(numComponents, vector<double>(numStages, 0.0));

    // Calculate Yij
    for (size_t j = 1; j < numStages; ++j)
    {
        double sumofall = 0.0;
        for (size_t i = 1; i < numComponents; ++i)
        {
            sumofall += vij[i][j];
        }
        for (size_t i = 1; i < numComponents; ++i)
        {
            Yij[i][j] = vij[i][j] / sumofall;
        }
    }
    return Yij;
}

inline double calculate_enthalpy(vector<double> &param, double temp)
{
    double q1 = param[0];
    double q2 = param[1] * temp;
    double q3 = param[2] * temp * temp;
    return q1 + q2 + q3;
}

inline vector<vector<double>> calculate_Hij(vector<vector<double>> &H_param, vector<double> &T)
{
    vector<vector<double>> Hij;

    // Assuming H_param and T have the correct dimension
    size_t numComponents = H_param.size();
    size_t numStages = T.size(); // Assuming all components have the same number of stages

    // Resize Hij to match the dimensions of vij
    Hij.resize(numComponents, vector<double>(numStages, 0.0));

    for (size_t j = 1; j < numStages; ++j)
    {
        for (size_t i = 1; i < numComponents; ++i)
        {
            Hij[i][j] = calculate_enthalpy(H_param[i], T[j]);
        }
    }
    return Hij;
}

inline vector<vector<double>> calculate_hij(vector<vector<double>> &h_param, vector<double> &T)
{
    vector<vector<double>> hij;

    // Assuming h_param and T have the correct dimension
    size_t numComponents = h_param.size();
    size_t numStages = T.size(); // Assuming all components have the same number of stages

    // Resize Hij to match the dimensions of vij
    hij.resize(numComponents, vector<double>(numStages, 0.0));

    for (size_t j = 1; j < numStages; ++j)
    {
        for (size_t i = 1; i < numComponents; ++i)
        {
            hij[i][j] = calculate_enthalpy(h_param[i], T[j]);
        }
    }
    return hij;
}

inline vector<double> calculate_Hi(vector<vector<double>> &Hij, vector<vector<double>> &Yij)
{
    vector<double> Hi;

    size_t numComponents = Hij.size();
    size_t numStages = Hij[0].size();
    Hi.resize(numStages,0.0);

    for(size_t j = 1; j<numStages; j++)
    {
        double sumofallcomp = 0.0;
        for(size_t i = 1; i<numComponents; i++)
        {
            sumofallcomp += (Yij[i][j]*Hij[i][j]);
        }
        Hi[j] = sumofallcomp;
    }
    return Hi;
}

inline vector<double> calculate_hi(vector<vector<double>> &hij, vector<vector<double>> &xij)
{
    vector<double> hi;

    size_t numComponents = hij.size();
    size_t numStages = hij[0].size();
    hi.resize(numStages,0.0);

    for(size_t j = 1; j<numStages; ++j)
    {
        double sumofallcomp = 0.0;
        for(size_t i = 1; i<numComponents; ++i)
        {
            sumofallcomp += (xij[i][j]*hij[i][j]);
        }
        hi[j] = sumofallcomp;
    }
    return hi;
}

vector<double> calculatehfi(vector<vector<double>> h_param, double& Feed_T)
{
    vector<double> hfi;

    size_t numComponents = h_param.size();
    size_t numStages = h_param[0].size();

    hfi.resize(numComponents,-1);

    hfi[1] = calculate_enthalpy(h_param[1],Feed_T);
    hfi[2] = calculate_enthalpy(h_param[2],Feed_T);
    return hfi;
}

vector<double> calculate_Vnew (vector<double> &Hi, vector<double> &hi, vector<double> &Li, vector<double> &Vi, double &distillate, double F1, double F2, vector<double> hfi)
{
    double Feed = F1+F2;
    double hf = ((F1*hfi[1])+(F2*hfi[2]))/Feed;
    vector<double> Vnew;

    size_t numComponents = Hi.size();
    
    Vnew.resize(numComponents,-1);

    Vnew = Vi;

    Vnew[3] = ((Vi[2]*Hi[2])-(distillate*hi[2])-(Li[1]*hi[1]))/(Hi[3]-hi[2]);

    Vnew[4] = ((Vi[3]*Hi[3])+(Li[4]*hi[3])-(Feed*(hf))-(Vi[3]*hi[2])+(distillate*hi[2]))/(Hi[4]-hi[3]);

    return Vnew;
}



#endif // ENERGY_SOLVER_H
