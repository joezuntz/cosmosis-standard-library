/* parse.h
*  =======
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
#include <vector>
#include <stdlib.h>
#include <string>
#include <cxxopts.hpp>

struct csmpars {
	int verbosity_level;
	std::vector<double> Omega_b;
    std::vector<double> Omega_m;
    std::vector<double> Sum_m_nu;
    std::vector<double> n_s;
    std::vector<double> h;
    std::vector<double> w_0;
    std::vector<double> w_a;
    std::vector<double> A_s;
    std::vector<int> n_redshift;
    std::vector< std::vector<double> > zvec;
	std::string outfilename;
	std::string outdir;
};

csmpars ee2_parser(int n_args, char * vec_args[]);
void read_cosmo_from_cmdline(cxxopts::ParseResult result, csmpars &CSM);
void read_classfile(std::string class_file_name, csmpars &CSM);
void read_cambfile(std::string class_file_name, csmpars &CSM);
void read_parfile(std::string par_file_name, csmpars &CSM);
