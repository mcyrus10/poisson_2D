#include "poisson_2D_processor.h"


T sinh(T x)
{
    return((exp(x)-exp(-x))/2.0);
}

T west_boundary_helmholtz(plint iY,
                SimulationParams<T> const& simParams)
{
    plint ny = simParams.getNy()-1;
    return(sinh(mu*(1.0-(T)iY/(T)ny))/sinh(mu));  // Equation 3.6
}

T east_boundary_helmholtz(plint iY,
                SimulationParams<T> const& simParams)
{
    plint ny = simParams.getNy()-1;
    pcout << "east --> iY = " << iY << "; y =" << (T)iY/(T)ny << "; "  << sinh(mu*(1.0-(T)iY/(T)ny))/sinh(mu) << endl;
    return(-sinh(mu*(1.0-(T)iY/(T)ny))/sinh(mu));  // Equation 3.6
}

T south_boundary_helmholtz(plint iX,
                SimulationParams<T> const& simParams)
{
    plint nx = simParams.getNx()-1;
    return(cos(pi*(T)iX/(T)nx));                // Equaiton 3.6
}

template<typename T>
class helmholtz_east{
    public: helmholtz_east(SimulationParams<T> simParams_)
        : simParams(simParams_)
        {}
    T operator()(plint iX, plint iY) const {
        return east_boundary_helmholtz(iY, simParams);
    }
private:
    SimulationParams<T> simParams;
};

template<typename T>
class helmholtz_west{
    public: helmholtz_west(SimulationParams<T> simParams_)
        : simParams(simParams_)
        {}
    T operator()(plint iX, plint iY) const {
        return west_boundary_helmholtz(iY, simParams);
    }
private:
    SimulationParams<T> simParams;
};

template<typename T>
class helmholtz_south{
    public: helmholtz_south(SimulationParams<T> simParams_)
        : simParams(simParams_)
        {}
    T operator()(plint iX, plint iY) const {
        return south_boundary_helmholtz(iX, simParams);
    }
private:
    SimulationParams<T> simParams;
};


void phi_setup( MultiBlockLattice2D<T,ADESCRIPTOR>& phiLattice,
                SimulationParams<T> & simParams)
{
    plint nx = simParams.getNx();
    plint ny = simParams.getNy();
    Box2D south_boundary(0,nx-1,0,0);
    Box2D west_boundary(0,0,0,ny-1);
    Box2D north_boundary(0,nx-1,ny-1,ny-1);
    Box2D east_boundary(nx-1,nx-1,0,ny-1);

    plint processorLevel = 5;

    initializeAtEquilibrium(    phiLattice,
                                phiLattice.getBoundingBox(),
                                simParams.getPhi_init(),
                                Array<T,2>((T)0.0,(T)0.0));

    OnLatticeAdvectionDiffusionBoundaryCondition2D<T, ADESCRIPTOR>* 
        bc = createLocalAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();

    bc->addTemperatureBoundary1N(south_boundary,phiLattice);
    setBoundaryDensity(phiLattice,south_boundary,helmholtz_south<T>(simParams));

    bc->addTemperatureBoundary1N(east_boundary,phiLattice);
    setBoundaryDensity(phiLattice,east_boundary,helmholtz_east<T>(simParams));

    bc->addTemperatureBoundary1N(west_boundary,phiLattice);
    setBoundaryDensity(phiLattice,west_boundary,helmholtz_west<T>(simParams));


    // North boundary == 0
    bc->addTemperatureBoundary1N(north_boundary,phiLattice);
    setBoundaryDensity(phiLattice,north_boundary,(T)0.0);



    integrateProcessingFunctional(  new poisson_2D<T,ADESCRIPTOR>(simParams),
                                    phiLattice.getBoundingBox(),
                                    phiLattice,
                                    processorLevel);


    phiLattice.initialize();

    delete bc;
}


int main(int argc, char* argv[])
//{{{
    {
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    // Instantiate variables
    SimulationParams<T> simParams = assign_params("params.xml");

    IncomprFlowParam<T> parameters( 1.0,
                                    1.0,
                                    simParams.getResolution(),
                                    simParams.getLx(),
                                    simParams.getLy());

    T phiOmega = 1.0/simParams.getTau_phi();
    pcout << "phiOmega = " << phiOmega << endl;
    plint nx = simParams.getNx();
    plint ny = simParams.getNy();

    pcout << "mu = " << mu << "; lambda = " << lambda << "; pi = " << pi << endl;
    T var = 0.0;
    for (int i = 0; i < 100; i+=1){
        var = (T)i/(T)simParams.getNx();
        pcout << "sinh(" << var << ") = " << sinh(var) << endl;;
    }
    writeLogFile(parameters,simParams,"output testing square domain");

    // Instantiate Charge Carrier Lattice
    MultiBlockLattice2D<T,ADESCRIPTOR> phiLattice(
                                    nx,
                                    ny,
                                    new ADYNAMICS<T, ADESCRIPTOR>(phiOmega));

    phi_setup(  phiLattice,
                simParams);

    plint convergenceIter = simParams.getConvergenceIter();
    // =========================================================================
    // Time Stepping
    // =========================================================================
    for (plint iT = 0; iT<=simParams.getMaxIter(); iT++)
    {
        if (iT%convergenceIter==0) {
            pcout << "Writing vtk_instance at iT = " << iT << endl;
            writeVTK(   phiLattice,
                        parameters,
                        simParams,
                        iT,
                        "phi Lattice");
        }
        phiLattice.collideAndStream();
    }

    plb_ofstream succesiveProfiles("concentration_final.dat");
    succesiveProfiles << std::setprecision(7)
        << *computeDensity(phiLattice,phiLattice.getBoundingBox()) << endl <<endl;
    /*
    */
}
// }}}
