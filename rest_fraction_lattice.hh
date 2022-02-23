#include "rest_fraction_lattice.h"
namespace plb {
namespace descriptors {

// AdvectionDiffusion D2Q5 //////////////////////////////////////////////
template<typename T>
const T rest_fraction_D2Q5Constants<T>::invD = (T)1 / (T) d;

template<typename T>
const int rest_fraction_D2Q5Constants<T>::vicinity = 1;

template<typename T>
const T rest_fraction_D2Q5Constants<T>::J0 = 0.0;

template<typename T>
const int rest_fraction_D2Q5Constants<T>::c
    [rest_fraction_D2Q5Constants<T>::q][rest_fraction_D2Q5Constants<T>::d] =
{
    { 0, 0},
    {-1, 0},
    {0, -1},
    {1,0}, 
    { 0,1}
};

template<typename T>
const int rest_fraction_D2Q5Constants<T>::cNormSqr[rest_fraction_D2Q5Constants<T>::q] =
{   0,
    1,
    1,
    1,
    1 
    };

template<typename T>
const T rest_fraction_D2Q5Constants<T>::t[rest_fraction_D2Q5Constants<T>::q] =
    {
        (T)0.0,
        (T)(1-J0)/4,
        (T)(1-J0)/4,
        (T)(1-J0)/4,
        (T)(1-J0)/4

    };

template<typename T>
const T rest_fraction_D2Q5Constants<T>::cs2 = (T)1 / (T)3;

template<typename T>
const T rest_fraction_D2Q5Constants<T>::invCs2 = (T)3;

template<typename T>
const char rest_fraction_Descriptor<T>::name[] = "rest_fraction_Descriptor";

   }  // namespace descriptors
}  // namespace plb
