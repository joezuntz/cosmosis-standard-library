/* main.cxx
*  =========
*  This file is part of EuclidEmulator2 
*  Copyright (c) 2020 Mischa Knabenhans
*
*  EuclidEmulator2 is free software: you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  EuclidEmulator2 is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <typeinfo>
#include <vector>
#include "emulator.h"
#include "parse.h"

int main(int argc, char *argv[]) {
    int verbosity_level = 0, n_cosmologies;
	const int nPCA = 15; //15 principal components
	std::vector<double> nlc, kmodes;
	csmpars CSM; // this is just a container for easier parsing of cosmological parameters
	string outfilename;


	/* GET INPUT PARAMETERS */
	CSM = ee2_parser(argc, argv);
	n_cosmologies = CSM.Omega_b.size();
	std::cout << "\nThere is/are " << n_cosmologies << " cosmologies specified in the input." << std::endl;

	/* INITIALIZE EE2 SESSION */
    EuclidEmulator ee2 = EuclidEmulator();
    std::cout << "\nEuclidEmulator2 >> Session started... " << std::endl;

	/* ITERATE THROUGH COSMOLOGIES */
	for (int cntr = 0; cntr < n_cosmologies; cntr++){
		/* Instatiate the Cosmology class */
		/*printf("Omega_b: ");
		printf("Omega_b[%d] = %f\n", cntr, CSM.Omega_b[cntr]);
		printf("Omega_m: ");
		printf("Omega_m[%d] = %f\n", cntr, CSM.Omega_m[cntr]);
		printf("Sum_m_nu: ");
		printf("Sum_m_nu[%d] = %f\n", cntr, CSM.Sum_m_nu[cntr]);
		printf("n_s: ");
		printf("n_s[%d] = %f\n", cntr, CSM.n_s[cntr]);
		printf("h: ");
		printf("h[%d] = %f\n", cntr, CSM.h[cntr]);
		printf("w_0: ");
		printf("w_0[%d] = %f\n", cntr, CSM.w_0[cntr]);
		printf("w_a: ");
		printf("w_a[%d] = %f\n", cntr, CSM.w_a[cntr]);
		printf("A_s: ");
		printf("A_s[%d] = %.5e\n", cntr, CSM.A_s[cntr]);
		printf("Entering cosmology class..\n");
		*/
		Cosmology cosmo = Cosmology(CSM.Omega_b[cntr], CSM.Omega_m[cntr], \
								    CSM.Sum_m_nu[cntr], CSM.n_s[cntr], CSM.h[cntr], \
									CSM.w_0[cntr], CSM.w_a[cntr], CSM.A_s[cntr]);

		//printf("Cosmology struct created successfully\n");
		/* compute NLCs for each cosmology */
		/*for(int i = 0; i < CSM.n_redshift[cntr]; i++){
			printf("z[%d][%d] = %f\n", cntr, i, CSM.zvec.at(cntr).at(i));
		}*/
		ee2.compute_nlc(cosmo, CSM.zvec[cntr], CSM.n_redshift[cntr]);
		//printf("NLC computed successfully.\n");
		/* Write result to file. Format: k [h/Mpc] B(k,z0) B(k,z1) ... B(k,zn) */
		string filename = CSM.outdir+CSM.outfilename+to_string(cntr)+".dat";
		std::cout << "Filename-->"<< filename << std::endl;
		ee2.write_nlc2file(filename, CSM.zvec[cntr], CSM.n_redshift[cntr]);
	}

	/* CLOSE EE2 SESSION (no explicit destructor required)*/
	std::cout << "EuclidEmulator2 >> Session closed... " << std::endl;
	return 0;
}
