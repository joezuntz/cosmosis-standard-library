#include "cosmosis/datablock/datablock.hh"
#include "cosmosis/datablock/section_names.h"
#include <vector>
#include <algorithm>
#include "gsl/gsl_spline.h"
#include "src/jla.h"
#include <cstring>

extern "C" {


int get_option(cosmosis::DataBlock * options, const std::string &name, std::string &parameter){
	
	auto status = options->get_val(OPTION_SECTION, name, parameter);
	if (status!=DBS_SUCCESS) {
		parameter = "";
		std::cerr<< "Could not find or understand parameter in JLA section: " << name << std::endl; 
		return 1;
	}
	return 0;
}

void * setup(cosmosis::DataBlock * options){
	int verbosity;
	int default_verbosity = 2;
	int error = 0;
	auto status = options->get_val(OPTION_SECTION, 
		"verbosity", default_verbosity, verbosity);
	if (status!=DBS_SUCCESS) {
		std::cerr<< "JLA Parameter error: verbosity not set right" << std::endl; 
		error|=status;
	}
	Configuration config;
	std::string parameter;

		
	options->get_val(OPTION_SECTION, "scriptmcut", config.scriptmcut);

	error|=get_option(options, "data_dir", parameter);
	std::string data_dir = parameter;
	if (data_dir!="") data_dir += "/";

	error|=get_option(options, "data_file", parameter);
	parameter = data_dir + parameter;
	config.data_file = strdup(parameter.c_str());

	error|=get_option(options, "mag_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C00 = strdup(parameter.c_str());

	error|=get_option(options, "stretch_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C11 = strdup(parameter.c_str());

	error|=get_option(options, "colour_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C22 = strdup(parameter.c_str());

	error|=get_option(options, "mag_stretch_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C01 = strdup(parameter.c_str());

	error|=get_option(options, "mag_colour_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C02 = strdup(parameter.c_str());

	error|=get_option(options, "stretch_colour_covmat_file", parameter);
	parameter = data_dir + parameter;
	config.C12 = strdup(parameter.c_str());

	if (error) exit(1);

	auto calculator = new JLALikelihood(verbosity);
	calculator->configure(config);


	return (void*) calculator;
}

gsl_spline * spline_from_vectors(std::vector<double> &x, std::vector<double> &y)
{
	int n = x.size();
	gsl_spline * s = gsl_spline_alloc(gsl_interp_akima, x.size());
	double* px = &x[0];
	double* py = &y[0];
	gsl_spline_init(s, px, py, n);
	return s;
}

#define CHECK_ERR(p) if (status!=DBS_SUCCESS) {std::cerr<< "JLA error loading required thing: " << p << std::endl; return status;}

int execute(cosmosis::DataBlock * block, void * config){
	auto calculator = (JLALikelihood*) config;

	// Space for Z and MU
	std::vector<double> Z;
	std::vector<double> MU;

	//Load Z
	auto status = block->get_val(DISTANCES_SECTION, "z", Z);
	CHECK_ERR("z");

	//Load MU
	status =  block->get_val(DISTANCES_SECTION, "mu", MU);
	CHECK_ERR("mu");

	// Load the four nuisance params
	double nuisance[4];
	const char * section = SUPERNOVA_PARAMS_SECTION;
	status =  block->get_val(section, "alpha", nuisance[0]);
	CHECK_ERR("alpha");

	status =  block->get_val(section, "beta", nuisance[1]);
	CHECK_ERR("beta");

	status =  block->get_val(section, "M", nuisance[2]);
	CHECK_ERR("M");

	status =  block->get_val(section, "deltaM", nuisance[3]);
	CHECK_ERR("deltaM");

	// Reverse if needed
	if (Z[1]<Z[0]) {std::reverse(Z.begin(), Z.end()); std::reverse(MU.begin(), MU.end());}

	// Get rid of the first element at z=0 since there mu=-inf
	Z.erase(Z.begin());
	MU.erase(MU.begin());


	// Make an interpolating spline
	gsl_spline * mu_of_z = spline_from_vectors(Z, MU);

	// Get the number of SN and their redshifts
	int n_sn = calculator->size();
	double * z_sn = calculator->getZ();
	double mu_sn[n_sn];

	// Evaluate our theory mu(z) for each supernova
	for (int i=0; i<n_sn; i++){
		mu_sn[i] = gsl_spline_eval(mu_of_z, z_sn[i], NULL);
	}

	gsl_spline_free(mu_of_z);

	//Get the chi^2 and convert to log-like
	double chi2 = calculator->computeLikelihood(mu_sn, nuisance);
	double like = -chi2/2.0;


	// Save the likelihood
	status = block->put_val(LIKELIHOODS_SECTION, "JLA_LIKE", like);
	if (status!=DBS_SUCCESS) return status;
	status = block->put_metadata(LIKELIHOODS_SECTION, "JLA_LIKE", "number_supernova", std::to_string(n_sn));

	return status;

}


} // end extern C