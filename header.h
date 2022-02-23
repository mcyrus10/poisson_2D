# include "output.h"


void coupleAdvectionDiffusion(  MultiBlockLattice2D<T,ADESCRIPTOR>& vII,
                                MultiBlockLattice2D<T,ADESCRIPTOR>& vIII,
                                MultiScalarField2D<int>& boundary_domain,
                                MultiScalarField2D<bool>& solid_domain,
                                SimulationParams<T> simParams)
{
    plint processorLevel = 4;
    plint ny = simParams.getNy();
    plint nx = simParams.getNx();
    T C_boundary_vII = simParams.getSoC()*simParams.getC_init();
    T C_boundary_vIII = (1.0-simParams.getSoC())*simParams.getC_init();
    //--------------------------------------------------------------------------
    // Constrict the domain so that the edges don't get counted
    //--------------------------------------------------------------------------
    //              lets say the domain is nx = 19 -> this means 
    //              the index goes from 0 to 18,so the nx-2 (17th) 
    //              element is the last element with an element to the right
    //                |      |
    //                |      |
    //                v      v
    Box2D domain(1,nx-2,1,ny-2);
    integrateProcessingFunctional(  new KangBoundaryExperiment<T,ADESCRIPTOR>(simParams), 
                                    //domain,
                                    vII.getBoundingBox(),
                                    vII,
                                    vIII,
                                    boundary_domain,
                                    processorLevel);
    pcout << "IntegrateProcessingFunctional Completed" << endl;
    // This OneCellIndexedFunctional Discriminates between solid and electrolyte
    // domain in order to have zero concentration in solid and NONZERO in
    // electolyte phase. Use initializeAtEquilibrium to get zero concentration
    // everywhere
    //applyIndexed(   vII,
    //                vII.getBoundingBox(),
    //                new ConcentrationInitializer<T,ADESCRIPTOR>(solid_domain,C_boundary_vII,simParams));
    //applyIndexed(   vIII,
    //                vIII.getBoundingBox(),
    //                new ConcentrationInitializer<T,ADESCRIPTOR>(solid_domain,C_boundary_vIII,simParams));

    defineDynamics( vII,
                    solid_domain,
                    new BounceBack<T,ADESCRIPTOR>, 
                    true);

    defineDynamics( vIII,
                    solid_domain,
                    new BounceBack<T,ADESCRIPTOR>, 
                    true);


    Array<T,2> u0((T)0,(T)0);                       // 2D velocity = 0 = (0,0)

    initializeAtEquilibrium(vII,
                            vII.getBoundingBox(),
                            (T)simParams.getC_init()*simParams.getSoC(),
                            u0);

    initializeAtEquilibrium(vIII,
                            vIII.getBoundingBox(),
                            (T)simParams.getC_init()*(1.0-simParams.getSoC()),
                            u0);

    vIII.initialize();
    vII.initialize();
    pcout << "Lattices Initalized" << endl;
}