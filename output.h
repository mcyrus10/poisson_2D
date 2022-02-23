# include "base.h" 

//{{{ void writeVTK
void writeVTK(  MultiBlockLattice2D<T,NSDESCRIPTOR>& lattice,
                IncomprFlowParam<T> const& parameters,
                plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk_velocity", iter, 6), dx);
    vtkOut.writeData<float>(*computeDensity(lattice), "density", 1.);
    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
}
//}}}
//{{{ void writeVTK
void writeVTK(  MultiBlockLattice2D<T,ADESCRIPTOR>& adLattice,
                IncomprFlowParam<T> & parameters,
                SimulationParams<T> & simParams,
                plint iter,
				string f_name)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
	plint nx = parameters.getNx();
	plint ny = parameters.getNy();
	plint resolution = parameters.getResolution();
	unique_ptr<MultiScalarField2D<T> > density = computeDensity(adLattice);
	unique_ptr<MultiTensorField2D<T,2> > gradient = computeGradient(*density);
	//unique_ptr<MultiScalarField2D<T> > flux = computeMassFlux(adLattice,boundary_domain,simParams);
    VtkImageOutput2D<T> vtkOut(createFileName(f_name.c_str(), iter, 6), dx);
    vtkOut.writeData<float>(*density, "density", 1.);
    vtkOut.writeData<2,float>(*gradient, "potential flux",resolution);
    vtkOut.writeData<float>(*computeNorm(*gradient), "norm flux",resolution);
    vtkOut.writeData<2,float>(*computeGradient(*computeNorm(*gradient)), "Laplacian",resolution);
    pcout << "Writing vtk_instance at iT = " << iter << endl;
    pcout << "Average Flux = " << computeAverage(*computeGradientNorm(*density), adLattice.getBoundingBox()) << endl;
}
//}}}
