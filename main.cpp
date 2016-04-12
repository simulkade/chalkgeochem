#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <vector>
#include "PhreeqcRM.h"
#include "IPhreeqc.hpp"
#include "IPhreeqcPhast.h"

using namespace std;

int main()
{
// libphreeqcrm.so
// libphreeqcrm.a
// I have one cell, I want to change the sulphate concentration and obtain the surface charges
// next step is to obtain the wettability
    int nxyz=1;
    int nthreads=3;
    int status;
    PhreeqcRM phreeqc_rm(nxyz, nthreads);
    phreeqc_rm.SetComponentH2O(true);
    phreeqc_rm.SetFilePrefix("surf_cmplx");
    phreeqc_rm.OpenFiles();

    // Set concentration units
    status = phreeqc_rm.SetUnitsSolution(2);      // 1, mg/L; 2, mol/L; 3, kg/kgs
    status = phreeqc_rm.SetUnitsPPassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    status = phreeqc_rm.SetUnitsExchange(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    status = phreeqc_rm.SetUnitsSurface(2);       // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    status = phreeqc_rm.SetUnitsGasPhase(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    status = phreeqc_rm.SetUnitsSSassemblage(1);  // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    status = phreeqc_rm.SetUnitsKinetics(1);      // 0, mol/L cell; 1, mol/L water; 2 mol/kg rock
    // Set conversion from seconds to user units (days)
    double time_conversion = 1.0 / 86400;
    status = phreeqc_rm.SetTimeConversion(time_conversion);
    // Set representative volume
    vector<double> rv;
    rv.resize(nxyz, 1.0);
    status = phreeqc_rm.SetRepresentativeVolume(rv);
    // Set initial porosity
    vector<double> por;
    por.resize(nxyz, 0.2);
    status = phreeqc_rm.SetPorosity(por);
    // Set initial saturation
    vector<double> sat;
    sat.resize(nxyz, 1.0);
    status = phreeqc_rm.SetSaturation(sat);

//here we map the transport grids to the chemistry solver
//in this case, I'm solving chemistry in all the grids.
//however, we can have fewer chemistry cells to save calculation time
    vector<int> grid2chem;
    grid2chem.resize(nxyz, -1);
    for (int i=0; i<nxyz; i++)
    {
        grid2chem[i]=i;
    }
    status= phreeqc_rm.CreateMapping(grid2chem);
    phreeqc_rm.LoadDatabase("phreeqc_oil.dat");
    // workers, initial phreeqc, utility
    phreeqc_rm.RunFile(true, true, true, "chalkoil.pqi");
    cout << status << endl;

    vector<int> ic1;
    ic1.resize(nxyz*7, -1);
    for (int i = 0; i < nxyz; i++)
    {
        ic1[i] = 1;              // Solution 1
        ic1[nxyz + i] = -1;      // Equilibrium phases none
        ic1[2*nxyz + i] = -1;     // Exchange none
        ic1[3*nxyz + i] = 1;    // Surface 1
        ic1[4*nxyz + i] = -1;    // Gas phase none
        ic1[5*nxyz + i] = -1;    // Solid solutions none
        ic1[6*nxyz + i] = -1;    // Kinetics none
    }
    status = phreeqc_rm.InitialPhreeqc2Module(ic1);
    cout << status;

    int ncomps = phreeqc_rm.FindComponents(); // accumulates a list of elements
    cout << endl << "Find componets output: "<< ncomps << endl;

// This is not necessary; just to show that the molar weight can be shown for elements
// GetComponents return the elements that are listed by the FindComponents call
    const vector<string> &components = phreeqc_rm.GetComponents();
    const vector<double> & gfw = phreeqc_rm.GetGfw();
    for (int i = 0; i < ncomps; i++)
    {
        ostringstream strm;
        strm.width(10);
        strm << components[i] << "    " << gfw[i] << "\n";
        phreeqc_rm.OutputMessage(strm.str());
        phreeqc_rm.ScreenMessage(strm.str()); // prints the message on the screen
    }
// prints a list of elements with their molecular weight

    cout << endl << "number of components: "<< components.size() << endl;
    // *** NOW THE IMPORTANT PART ***
    // working with concentrations
    // we run the model once to equilibrate the system at time zero
    vector<double> c;
    phreeqc_rm.SetTime(0.0);
    phreeqc_rm.SetTimeStep(0.0);
    phreeqc_rm.RunCells();
    phreeqc_rm.GetConcentrations(c);

    cout << endl << "number of components: "<< components.size() << endl;

    for (int i = 0; i < ncomps; i++)
    {
        ostringstream strm;
        strm.width(10);
        strm << components[i] << "    " << c[i] << "\n";
        //phreeqc_rm.OutputMessage(strm.str());
        phreeqc_rm.ScreenMessage(strm.str()); // prints the message on the screen
    }

// working with species
    //status = phreeqc_rm.InitialPhreeqc2Module(ic1);
    status = phreeqc_rm.SetSelectedOutputOn(true);
    status = phreeqc_rm.SetSpeciesSaveOn(true);
    ncomps = phreeqc_rm.FindComponents();
    int nspecies = phreeqc_rm.GetSpeciesCount();
    const vector<string> &species = phreeqc_rm.GetSpeciesNames();
    vector<double> cs;
    phreeqc_rm.RunCells();
    phreeqc_rm.GetSpeciesConcentrations(cs);
    cout << endl << "number of aqueous species: "<< nspecies << endl;
    for (int i = 0; i < nspecies; i++)
    {
        ostringstream strm;
        strm.width(10);
        strm << species[i] << "    " << cs[i] << "\n";
        //phreeqc_rm.OutputMessage(strm.str());
        phreeqc_rm.ScreenMessage(strm.str()); // prints the message on the screen
    }

// one issue is that the surface species are not among the species. we have to obtain
//their values with getselectedoutput function. This is the way it works.
    status = phreeqc_rm.SetCurrentSelectedOutputUserNumber(1); //SELECTED_OUTPUT user number
    vector<double> so;
    int col = phreeqc_rm.GetSelectedOutputColumnCount();
    status = phreeqc_rm.GetSelectedOutput(so);
    vector<string> headings;
    headings.resize(col);
    std::cerr << "     Selected output: " << "\n";
    for (int j = 0; j < col; j++)
    {
      status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
      std::cerr << "          " << j << " " << headings[j] << ": " << so[j] << "\n";
    }

// Now I'm going to plot the surface charge changes versus sulfate concentration
    status = phreeqc_rm.SetSpeciesSaveOn(false);
    double t=0;
    int n_point=10;
    vector<double> s_charge;
    s_charge.resize(n_point, -1);
    for (int i=0; i<n_point; i++)
    {
        t+=0.1;
        phreeqc_rm.SetTime(t);
        phreeqc_rm.SetTimeStep(0.0);
        c[9]+=0.000005; // sulfur concentration
        phreeqc_rm.SetConcentrations(c);
        phreeqc_rm.RunCells();
        status = phreeqc_rm.GetSelectedOutput(so);
        s_charge[i]=so[0]/(so[0]+so[1]+so[2]); //sum of surface charge
        cout << s_charge[i] << endl;
//        for (int j = 0; j < col; j++)
//        {
//          status = phreeqc_rm.GetSelectedOutputHeading(j, headings[j]);
//          std::cerr << "          " << j << " " << headings[j] << ": " << so[j] << "\n";
//        }

    }

    for (int i = 0; i < ncomps; i++)
    {
        ostringstream strm;
        strm.width(10);
        strm << components[i] << "    " << c[i] << "\n";
        //phreeqc_rm.OutputMessage(strm.str());
        phreeqc_rm.ScreenMessage(strm.str()); // prints the message on the screen
    }

}
