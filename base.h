#include "palabos2D.h"
#include "palabos2D.hh"
#include "rest_fraction_lattice.hh"
#include <cmath>

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::D2Q9Descriptor
#define ADESCRIPTOR descriptors::rest_fraction_Descriptor
#define ADYNAMICS AdvectionDiffusionBGKdynamics

// Gas Constant and Faraday's Constant
T const Ru = 8.3144598;
T const F = 96485.33289;
T const lambda = 2.0; 						// below equation 3.6
T const pi = 4.0*atan(1.0); 				// below equation 3.6
T const mu = pow(lambda*lambda+pi*pi,0.5);	// below equation 3.6

//------------------------------------------------------------------------------
//  These are debugging functions for outputting they are separate from IO so
//  that IO can be defined downstream
//------------------------------------------------------------------------------
// {{{ printPopulations
void printPopulations(	plint iX,
						plint iY,
						BlockLattice2D<T,ADESCRIPTOR>& adLattice,
						string popName)
{
	Array<T,5>& cell = adLattice.get(iX,iY).getRawPopulations();
	pcout << "printing " << popName << " @ ix,iy = " << iX << "," << iY << endl;
	T density = 0.0;
	for (int i = 0; i < 5; i++)
	{
		pcout << "i = " << i << "; population = " << cell[i]+ADESCRIPTOR<T>::t[i] << endl;
		density+=cell[i]+ADESCRIPTOR<T>::t[i];
	}
	pcout << "concentration =  " << density << endl;
}
// }}}
// {{{  print_array 
template<int nrows, int ncols>
void print_array(   T const X[nrows][ncols],
                    string arr_name){
    pcout << arr_name << " = " << endl;
    for (int i = 0; i<nrows; i++){
        for (int j = 0; j<ncols; j++){
            pcout << "--> " << arr_name << "[" << i << "][" << j << "] = " << X[i][j] << endl;
        }
    }
}
// }}}


//------------------------------------------------------------------------------
// Essentially All the methods need SimulationParams so it is defined at the top
// of the stream as well as its I/O functions
//------------------------------------------------------------------------------
// {{{ start SimulationParams
template<typename T>
class SimulationParams
{
	private:
		plint lx, ly, resolution, maxIter, convergenceIter, nx, ny;
		T phi_init, K_0, Temperature, epsilon, tau_phi, inletFlux, conductivity;
	public:
		SimulationParams<T>(
					plint lx_,
					plint ly_,
					plint resolution_,
					plint maxIter_,
					plint convergenceIter_,
					T phi_init_,
					T K_0_,
					T Temperature_,
					T epsilon_,
					T tau_phi_,
					T inletFlux_
					):
			lx(lx_),
			ly(ly_),
			resolution(resolution_),
			maxIter(maxIter_),
			convergenceIter(convergenceIter_),
			phi_init(phi_init_),
			K_0(K_0_),
			Temperature(Temperature_),
			epsilon(epsilon_),
			tau_phi(tau_phi_),
			inletFlux(inletFlux_)
		{
			//Calculated values
			nx = lx*resolution+1;
			ny = ly*resolution+1;
			conductivity = (tau_phi-0.5)/3;
		}
		//Methods
		plint	getLx() const { return lx; }
		plint	getLy() const { return ly; }
		plint	getResolution() const { return resolution; }
		plint	getMaxIter() const { return maxIter; }
		plint	getConvergenceIter() const { return convergenceIter; }
		plint	getNx() const { return nx; }
		plint	getNy() const { return ny; }
		T		getPhi_init() const { return phi_init; }
		T		getK_0() const { return K_0; }
		T		getTemperature() const { return Temperature; }
		T		getEpsilon() const { return epsilon; }
		T		getTau_phi() const { return tau_phi; }
		T		getInletFlux() const { return inletFlux; }
		T		getConductivity() const { return conductivity; }
};
// }}} end SimulationParams
// {{{ start assign_params
SimulationParams<T> assign_params(string f_name)
{
	plint lx, ly, resolution, maxIter, convergenceIter;
	T phi_init, K_0, Temperature, epsilon, tau_phi, inletFlux;
	try{
		XMLreader xmlFile(f_name);
		pcout << "CONFIGURATION" << endl;
		pcout << "=============" << endl;
		xmlFile.print(0);
		xmlFile["inputs"]["lx"].read(lx);
		xmlFile["inputs"]["ly"].read(ly);
		xmlFile["inputs"]["resolution"].read(resolution);
		xmlFile["inputs"]["maxIter"].read(maxIter);
		xmlFile["inputs"]["convergenceIter"].read(convergenceIter);
		xmlFile["inputs"]["phi_init"].read(phi_init);
		xmlFile["inputs"]["K_0"].read(K_0);
		xmlFile["inputs"]["Temperature"].read(Temperature);
		xmlFile["inputs"]["epsilon"].read(epsilon);
		xmlFile["inputs"]["tau_phi"].read(tau_phi);
		xmlFile["inputs"]["inletFlux"].read(inletFlux);
		pcout << "=============" << endl << endl;
		} catch (PlbIOException& exception) { 
		    pcout << exception.what() << endl;
		}
		return SimulationParams<T> (
					lx,
					ly,
					resolution,
					maxIter,
					convergenceIter,
					phi_init,
					K_0,
					Temperature,
					epsilon,
					tau_phi,
					inletFlux);
};
//}}} end assign_params
// {{{ start writeLogFile
template<typename T>
void writeLogFile(  IncomprFlowParam<T> const& parameters,
		SimulationParams<T> const& simParams,
		std::string const& title)
{
	std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
	plb_ofstream ofile(fullName.c_str());
	ofile << title << "\n\n";
	ofile << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
	ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
	ofile << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
	ofile << "Relaxation frequency:      omega=" << parameters.getOmega() << "\n";
	ofile << "LatticeViscosity:          nu=" << parameters.getLatticeNu() << "\n";
	ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
	ofile << "nx:                        nx=" << parameters.getNx() << "\n";
	ofile << "ny:                        ny=" << parameters.getNy() << "\n";
	ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
	ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
	ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
	ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
	ofile << "==============================" << "\n";
	ofile << "lx			lx = " << simParams.getLx() << "\n";
	ofile << "ly			ly = " << simParams.getLy() << "\n";
	ofile << "resolution			resolution = " << simParams.getResolution() << "\n";
	ofile << "maxIter			maxIter = " << simParams.getMaxIter() << "\n";
	ofile << "convergenceIter			convergenceIter = " << simParams.getConvergenceIter() << "\n";
	ofile << "nx			nx = " << simParams.getNx() << "\n";
	ofile << "ny			ny = " << simParams.getNy() << "\n";
	ofile << "phi_init			phi_init = " << simParams.getPhi_init() << "\n";
	ofile << "K_0			K_0 = " << simParams.getK_0() << "\n";
	ofile << "Temperature			Temperature = " << simParams.getTemperature() << "\n";
	ofile << "epsilon			epsilon = " << simParams.getEpsilon() << "\n";
	ofile << "tau_phi			tau_phi = " << simParams.getTau_phi() << "\n";
	ofile << "inletFlux			inletFlux = " << simParams.getInletFlux() << "\n";
	ofile << "conductivity			conductivity = " << simParams.getConductivity() << "\n";
}
// }}} end writeLogFile



