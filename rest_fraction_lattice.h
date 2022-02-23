#include "../../src/core/globalDefs.h"
#include "../../src/latticeBoltzmann/externalFields.h"
#include "../../src/latticeBoltzmann/roundOffPolicy.h"
#include "../../src/latticeBoltzmann/nearestNeighborLattices3D.h"
#include <vector>

namespace plb {

namespace descriptors
{
//===========================================================================//
//===================   Rest Fraction Lattice Descriptors   =================//
//===========================================================================//
    
    /// D2Q5 lattice
    template <typename T> 
    struct rest_fraction_D2Q5Constants
    {
        enum { d = 2, q = 5 };          ///< number of dimensions/distr. functions
        static const T invD;            ///< 1 / (number of dimensions)
        static const int vicinity;      ///< size of neighborhood
        static const int c[q][d];       ///< lattice directions
        static const int cNormSqr[q];   ///< norm-square of the vector c
        static const T t[q];            ///< lattice weights
        static const T cs2;             ///< lattice constant cs2 (in BGK, this is the square-speed-of-sound)
        static const T invCs2;          ///< 1 / cs2
        static const T J0;              ///Rest Fraction
    };

    template <typename T> 
    struct rest_fraction_D2Q5DescriptorBase: public rest_fraction_D2Q5Constants<T>, public DefaultRoundOffPolicy<T>
    {
        typedef D2Q5DescriptorBase<T> BaseDescriptor;
        enum { numPop=rest_fraction_D2Q5Constants<T>::q };
    };
    
    /// AD D2Q5 lattice
    template <typename T> 
    struct rest_fraction_Descriptor : public rest_fraction_D2Q5DescriptorBase<T>, public descriptors::Velocity2dDescriptorBase
    {
        static const char name[];
    };
    
    
}  // namespace descriptor_expansion

}  // namespace plb


