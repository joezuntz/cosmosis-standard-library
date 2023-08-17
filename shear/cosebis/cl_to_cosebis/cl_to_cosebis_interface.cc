///for c++ use the .hh and for c the .h version
//This deals with the inputs and outputs
#include "datablock/datablock.hh"
//This is just a header file which defines the different section names
#include "datablock/section_names.h"
#include <typeinfo>

/*CosmoSIS interface file for going from shear C(l) to E/B - COSEBIs
*/
#include "COSEBIs.h"

//could add an option to read B-modes in and use it to constrain c-term

extern "C" {

	const string shear_cl = SHEAR_CL_SECTION;
	const string cterm_section= "shear_c_bias";

	typedef struct COSEBIs_config 
	{
		string input_section_name; //input section name the default is: shear_cl
		string output_section_name;// where cosebis results are saved, default: cosebis
		int n_max; // largest n-mode to be caluclated for COSEBIs: 1,2,...,n_max
		number theta_min; //theta_min
		number theta_max; //theta_max

		COSEBIs *cosebis; //cosebis object
		
		//Constant c-term modelling. Will need input_sin4phi and input_cos4phi either through a file or in the block 
		string cterm_section;
		bool add_c_term; // Set to true if constant c-term modelling is required
		matrix En_cos4phi;// needed for constant cterm modelling 
		matrix En_sin4phi;// needed for constant cterm modelling 
		bool inputSinCos_given; //if input files are given in the ini file then will look at the matrices above
		//default is cosebis_cos4phi
		string input_cos4phi_section_name; //if input files not given will look in the block in this section for cos4phi
		//default is cosebis_sin4phi
		string input_sin4phi_section_name; //if input files not given will look in the block in this section for sin4phi

		//2D c-term modelling. Will need input_2Dcterm either through a file: input_2Dcterm_filename or in the block
		bool add_2D_cterm;
		matrix En_2D;
		bool input2D_given; // if the input file is given in the ini file then will look at the matrix above
		//deafult is cosebis_2D
		string input_2D_section_name; //if input files not given will look in the block in this section for En_2D

	} COSEBIs_config;


	///define a structure with everything that is needed to be read in setup and sent to execute
	//Some of these are read from the ini file. For example input_section_name and n_max
	//Some are initialized in the setup, such as cosebis. COSEBIs is a class that produces En/Bn and 
	//their covariance, etc. 
	int get_option(cosmosis::DataBlock * options, const string &name, string &parameter)
	{
	
		auto status = options->get_val(OPTION_SECTION, name, parameter);
		if (status!=DBS_SUCCESS) 
		{
			parameter = "";
			cerr<< "Could not find or understand parameter in cosebis section: " 
				<< name << std::endl; 
			return 1;
		}
		return 0;
	}
	
	void * setup(cosmosis::DataBlock * options, cosmosis::DataBlock * block)
	{
		//options reads the ini file
		//define config here and then read from options the relevant input quantities
		COSEBIs_config * config = new COSEBIs_config;

  		string sectionName=OPTION_SECTION;

		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		clog<<endl<<endl;
		clog<<"*********in COSEBIs interface setup*********"<<endl;
		clog<<endl;

		config->cterm_section=cterm_section;

		//get input section name, default= shear_cl
		status=options->get_val<string>(sectionName, string("input_section_name"), config->input_section_name);
		if (status) 
		{
			clog<<"You did not define the input_section_name,";
			clog<<"setting to the default: shear_cl"<<endl;
			config->input_section_name=string("shear_cl");
		}
		else
			clog<<"The input_section_name:"<<config->input_section_name<<endl;

		//get output_section_name, default= "cosebis"
		status=options->get_val<string>(sectionName, string("output_section_name"), config->output_section_name);
		if (status) 
		{
			clog<<"You did not define the out_section_name, ";
			clog<<"setting to the default: cosebis"<<endl;
			config->output_section_name=string("cosebis");
		}
		else
			clog<<"The output_section_name:"<<config->output_section_name<<endl;

		//get theta_min value, default=0.5
		status=options->get_val<number>(sectionName, string("theta_min"), config->theta_min);
		if (status) 
		{
			clog<<"You did not define the COSEBIs range theta_min,"<<endl;
			clog<<"setting to the default value of"<<0.5<<endl;
			config->theta_min=0.5;
		}
		else
			clog<<"COSEBIs range theta_min="<<config->theta_min<<endl;

		//get theta_max value, default=300
		status=options->get_val<number>(sectionName, string("theta_max"), config->theta_max);
		if (status) 
		{
			clog<<"You did not define the COSEBIs range theta_max,"<<endl;
			clog<<"setting to the default value of"<<300.<<endl;
			config->theta_max=300.;
		}
		else
			clog<<"COSEBIs range theta_max="<<config->theta_max<<endl;

		//get n_max value, default=10
		status=options->get_val<int>(sectionName, string("n_max"), config->n_max);
		if (status) 
		{
			clog<<"You did not define the number of COSEBIs modes that you require, n_max,"<<endl;
			clog<<"setting to the default value of"<<5<<endl;
			config->n_max=5;
		}
		else
			clog<<"Calculating COSEBI nodes up to n_max="<<config->n_max<<endl;

		//To add the c-term effect we either need input files or will read from the block if add_c_term is true
		string input_cterm_cos4phi_filename,input_cterm_sin4phi_filename;
		int add_c_term=0;
		config->inputSinCos_given=false;
		config->add_c_term=false;
		status=options->get_val<int>(sectionName, string("add_c_term"), add_c_term);
		if(add_c_term<1)
		{
			clog<<"This calculation does not include a spatially varying c_term correction"<<endl;
			config->add_c_term=false;
		}
		else
		{
			status=options->get_val<string>(sectionName, string("cterm_section_name"), config->cterm_section);
			if(status)
			{
				clog<<"You've requested a c_term correction without providing the cterm_section_name.  Setting it to the default value:";
				clog<<config->cterm_section<<endl;
			}
			else
			{
				clog<<"The cterm_section_name: "<<config->cterm_section<<endl;
			}
			clog<<"This calculation includes a spatially varying c_term correction using the "<< config->cterm_section<<" section for values"<<endl;
			config->add_c_term=true;
			config->input_cos4phi_section_name="cosebis_cos4phi";

			//read in the cterm section name for cos4phi, default is cosebis_cos4phi
			status=options->get_val<string>(sectionName, string("input_cos4phi_section_name"), config->input_cos4phi_section_name);
			if(status)
			{
				clog<<"The c-term input_cos4phi_section_name was not given.  Setting it to the default value:";
				clog<<config->input_cos4phi_section_name<<endl;
			}
			else
			{
				clog<<"The input_cos4phi_section_name: "<<config->input_cos4phi_section_name<<endl;
			}

			//read in the cterm section name for sin4phi, default is cosebis_sin4phi
			config->input_sin4phi_section_name="cosebis_sin4phi";
			status=options->get_val<string>(sectionName, string("input_sin4phi_section_name"), config->input_sin4phi_section_name);
			if(status)
			{
				clog<<"The c-term input_sin4phi_section_name was not given.  Setting it to the default value:";
				clog<<config->input_sin4phi_section_name<<endl;
				
			}
			else
			{
				clog<<"The input_sin4phi_section_name: "<<config->input_sin4phi_section_name<<endl;
			}

			//get input cos4phi cterm file name, if not given look into the block
			status=options->get_val<string>(sectionName, string("InputCos4phi"), input_cterm_cos4phi_filename);
			if (status) 
			{
				clog<<"Could not find the name of the input_cterm_cos4phi_filename for COSEBIs file";
				clog<<", going to look in the block at execute level"<<endl;
				config->inputSinCos_given=false;
			}
			else
			{
				clog<<"The name of the input_cterm_cos4phi_filename for COSEBIs file=";
				clog<<input_cterm_cos4phi_filename<<endl;
				config->inputSinCos_given=true;
				///read En_cos4phi from file
				config->En_cos4phi.readFromASCII_marika((input_cterm_cos4phi_filename).c_str());
				
				//get input cos4phi cterm file name
				status=options->get_val<string>(sectionName, string("InputSin4phi"), input_cterm_sin4phi_filename);
				if (status)
				{ 
					clog<<"Could not find the name of the input_cterm_sin4phi_filename for COSEBIs file"<<endl;
					clog<<", going to look in the block at execute level"<<endl;
					config->inputSinCos_given=false;
				}
				else
				{
					clog<<"The name of the input_cterm_sin4phi_filename for COSEBIs file=";
					clog<<input_cterm_sin4phi_filename<<endl;
					config->inputSinCos_given=true;
					///read En_sin4phi from file
					config->En_sin4phi.readFromASCII_marika((input_cterm_sin4phi_filename).c_str());
				}
			}
		}

		//To add 2D cterm modelling set add_2D_cterm to true
		int add_2D_cterm=0;
		config->add_2D_cterm=false;
		config->input2D_given=false;
		status=options->get_val<int>(sectionName, string("add_2D_cterm"), add_2D_cterm);
		if(add_2D_cterm>0)
		{
			config->add_2D_cterm=true;
			config->input_2D_section_name="cosebis_2D";
			status=options->get_val<string>(sectionName, string("cterm_section_name"), config->cterm_section);
			if(status)
			{
				clog<<"cterm_section_name not given setting it to the default value:";
				clog<<config->cterm_section<<endl;
			}
			else
			{
				clog<<"got the cterm_section_name: "<<config->cterm_section<<endl;
			}
			clog<<"using c_term correction, will look in "<< config->cterm_section<<" section for values"<<endl;

			//input section name for 2D cterm, default is cosebis_2D
			status=options->get_val<string>(sectionName, string("input_2D_section_name"), config->input_2D_section_name);
			if(status)
			{
				clog<<"input section name for cosebis_2D not given setting it to the default value:";
				clog<<config->input_2D_section_name<<endl;
			}
			else
			{
				clog<<"got the input section name for cosebis_2D: "<<config->input_2D_section_name<<endl;
			}
			

			string input_2Dcterm_filename;
			status=options->get_val<string>(sectionName, string("input_2Dcterm_filename"), input_2Dcterm_filename);
			if (status)
			{ 
				clog<<"Could not find the name of the input_2Dcterm_filename for COSEBIs file";
				clog<<"going to look in the block at execute level"<<endl;
				config->input2D_given=false;
			}
			else
			{
				clog<<"got the name of the input_2Dcterm_filename for COSEBIs file="<<input_2Dcterm_filename<<endl;
				config->input2D_given=true;
				///read input_2Dcterm for COSEBIs from file
				config->En_2D.readFromASCII_marika((input_2Dcterm_filename).c_str());
			}
		}

		//get Wn, Tn and output Tn folder names
		string WnFolderName,TnFolderName,OutputTnFolderName,WnFileName;
		WnFolderName = COSEBIS_DIR  "WnLog/";
		TnFolderName = COSEBIS_DIR  "TLogsRootsAndNorms/";
		OutputTnFolderName= COSEBIS_DIR  "TpnLog/";
		WnFileName   = "WnLog";

		// int precision = 20;
		// int Nlbins    = 1000000;

		int precision = 10;
		int Nlbins    = 100000;

		clog << "WnFolderName = " << WnFolderName << endl;

		status=options->get_val<string>(sectionName, string("Wn_Output_FolderName"), WnFolderName);
		if(status)
		{
			clog<<"Could not find WnLog folder name in Wn_Output_FolderName, ";
			clog<<"setting to default: "<<WnFolderName<<endl;
		}
		else
			clog<<"WnLog folder name is:"<<WnFolderName<<endl;

		status=options->get_val<string>(sectionName, string("Wn_file_name"), WnFileName);
		if(status)
		{
			clog<<"Could not find Wn file name in Wn_file_name, ";
			clog<<"setting to default: "<<WnFileName<<endl;
		}
		else
			clog<<"Wn file name is:"<<WnFileName<<endl;

		status=options->get_val<int>(sectionName, string("table_precision"), precision);
		if(status)
		{
			clog<<"Could not find the number of digits in table_precision, ";
			clog<<"setting to default: "<<precision<<endl;
		}
		else
			clog<<"table_precision is:"<<precision<<endl;

		status=options->get_val<int>(sectionName, string("number_of_Wn_l_bins"), Nlbins);
		if(status)
		{
			clog<<"Could not find number_of_Wn_l_bins, ";
			clog<<"setting to default: "<<Nlbins<<endl;
		}
		else
			clog<<"number_of_Wn_l_bins is:"<<Nlbins<<endl;


		status=options->get_val<string>(sectionName, string("Roots_n_Norms_FolderName"), TnFolderName);
		if(status)
		{
			clog<<"Could not find Root and Norms folder name in Roots_n_Norms_FolderName,"; 
			clog<<" setting to default: "<<TnFolderName<<endl;
		}
		else
			clog<<"Root and Norms folder name is:"<<TnFolderName<<endl;


		status=options->get_val<string>(sectionName, string("Tn_Output_FolderName"), OutputTnFolderName);
		if(status)
		{
			clog<<"Could not find T_pm folder name in Tn_Output_FolderName, ";
			clog<<"setting to default: "<<OutputTnFolderName<<endl;
	  	}
		else
			clog<<"T_pm folder name is:"<<OutputTnFolderName<<endl;

		//   initialize COSEBIs
		COSEBIs *cosebis = new COSEBIs();
		cosebis->initialize(config->n_max,config->theta_min,config->theta_max,1 //npair set to one for now, will be set seperately in execute to the correct value
				,WnFolderName,TnFolderName,OutputTnFolderName,WnFileName,precision);

		cosebis->setWns(config->n_max,Nlbins);
		config->cosebis=cosebis;
		return (void *) config;
		// config is sent to execute 
	}

	DATABLOCK_STATUS execute(cosmosis::DataBlock *block, void *config_in) 
	{
		enable_gsl_error_handling();
		// Config is whatever you returned from setup above
		// Block is the collection of parameters and calculations for
		// this set of cosmological parameters
		COSEBIs_config *config= (COSEBIs_config*) config_in;
		DATABLOCK_STATUS status = (DATABLOCK_STATUS)0;
		const DATABLOCK_STATUS failure = (DATABLOCK_STATUS)1;

		//lets save some of the atributes to block
		status = block->put_val<double>(config->output_section_name, string("theta_max"), config->theta_max);
		status = block->put_val<double>(config->output_section_name, string("theta_min"), config->theta_min);

		//get the number of redshift bins from cosmosis
		int num_z_bin_A;
		int num_z_bin_B;
		status = block->get_val(config->input_section_name, string("nbin_a"), num_z_bin_A);
		if(status)
		{
			status = block->get_val(config->input_section_name, string("nbin"), num_z_bin_A);
			if(status)
			{
				clog<<"looked for nbin_a first and couldn't find it. Then looked for nbin and it wasn't there"<<endl;
				return status;
			}
			num_z_bin_B = num_z_bin_A;
		}
		else
		{
			status = block->get_val(config->input_section_name, string("nbin_b"), num_z_bin_B);
			if(status)
			{
				clog<<"looked for nbin_b it is not set"<<endl;
				return status;
			}
		}

		vector<number> ell,logell;
		//get ell vector
		status = block->get_val(config->input_section_name, string("ell"), ell);
		if (status) 
		{
			clog<<"Could not load ell in C_ell to COSEBIs"<<endl;
			return status;
		}

		int nell=ell.size();
		//make logell to send to COSEBIs
		for(int i=0; i<nell; i++)
		{
			logell.push_back(log(ell[i]));
		}

		int n_max=config->n_max;
		vector<int> n_vals(n_max);
		for(int n=0;n<n_max;n++)
		 	n_vals[n]=n+1;

		status = block->put_val<vector<int> >(config->output_section_name, string("n"), n_vals);


		//here is where we read in the Cls and calculate COSEBIs
		int start_ind=0;
		int start_ind_2D=0;
		for (int i_bin=1; i_bin<=num_z_bin_A; i_bin++) 
		{
			for (int j_bin=1; j_bin<=num_z_bin_B; j_bin++) 
			{
				// read in C(l)
				vector<number> C_ell;
				string name_in=string("bin_")+toString(j_bin)+string("_")+toString(i_bin);
				//check if C(l) exists for this z-bin combination
				bool has_val = block->has_val(config->input_section_name, name_in);
				if (has_val) 
				{
					//if C(l) exists then read in
					status = block->get_val<vector<number> >(config->input_section_name, name_in, C_ell);
					
					matrix En_mat;

					//this also extrapolates the input Cl to LHIGH
					config->cosebis->setPower_single_withExtrapolation(logell,C_ell);

					//Calculate COSEBIs for the given C(l)
					En_mat=config->cosebis->calEn();
					

					//now add c-term contributions
					if(config->add_c_term)
					{
						number c1=0.;
						number c2=0.;

						status=block->get_val(config->cterm_section, "c1", c1);
						if(status)
						{
							clog<<"could not read c1 from "<<config->cterm_section<<endl;
							return status;
						}

						status=block->get_val(config->cterm_section, "c2", c2);
						if(status)
						{
							clog<<"could not read c2 from "<<config->cterm_section<<endl;
							return status;
						}

					 	matrix En_sin4phi(n_max),En_cos4phi(n_max);
					 	if(config->inputSinCos_given)
						{
						 	for(int ind=0; ind<n_max; ind++)
						 	{
						 		En_cos4phi.load(ind,config->En_cos4phi.get(ind+start_ind));
						 		En_sin4phi.load(ind,config->En_sin4phi.get(ind+start_ind));
						 	}
						 	start_ind+=n_max;
						}
						else
						{
							bool has_val1 = block->has_val(config->input_cos4phi_section_name, name_in);
							bool has_val2 = block->has_val(config->input_sin4phi_section_name, name_in);
							vector<number> cos4phi_vec, sin4phi_vec;
							if (has_val1 & has_val2) 
							{
								status = block->get_val<vector<number> >(config->input_cos4phi_section_name, name_in, cos4phi_vec);
								status = block->get_val<vector<number> >(config->input_sin4phi_section_name, name_in, sin4phi_vec);
								//now check if the size is correct, needs to have at least the same number of n-modes
								if((cos4phi_vec.size()>=n_max) &(sin4phi_vec.size()>=n_max))
								{
									for(int i=0; i<n_max; i++)
									{
										En_cos4phi.load(i,cos4phi_vec[i]);
										En_sin4phi.load(i,sin4phi_vec[i]);
									}
								}
							}
							else
							{
								clog<<"****WARNING no value for cos4phi or sin4phi given in the block, exiting now ..."<<endl;
								exit(1);
							}
						}
					 	//this is the additional contribution due to constant cterm see appendix D of https://arxiv.org/abs/2007.15633
						En_mat+=(c1*c1-c2*c2)*En_cos4phi+2.*c1*c2*En_sin4phi;
					}

					if(config->add_2D_cterm)
					{
						
						number Ac=0.;
						status=block->get_val(config->cterm_section, "Ac", Ac);
						if(status)
						{
							clog<<"could not read Ac from "<<config->cterm_section<<endl;
							return status;
						}

						matrix En_2D(n_max);
						if(config->input2D_given)
						{
							for(int ind=0; ind<n_max; ind++)
						 	{
						 		En_2D.load(ind,config->En_2D.get(ind+start_ind_2D));
						 	}
				        	start_ind_2D+=n_max;
						}
						else
						{
							bool has_val = block->has_val(config->input_2D_section_name, name_in);
							vector<number> En_2D_vec;
							if (has_val) 
							{
								status = block->get_val<vector<number> >(config->input_2D_section_name, name_in, En_2D_vec);
								//now check if the size is correct
								if(En_2D_vec.size()>=n_max)
								{
									for(int i=0; i<n_max; i++)
									{
										En_2D.load(i,En_2D_vec[i]);
									}
								}
							}
							else
							{
								clog<<"*****WARNING no value for xi_2D is given in the block exiting now..."<<endl;
								exit(1);
							}
				        }
				        //this is the contribution due to 2D cterm see section 5.1.2 of https://arxiv.org/abs/1810.02353
						En_mat+=Ac*En_2D;
					}

					// turning the matrix into a vector to be saved in block
					vector<number> En_vec(En_mat.rows);
					for(int m=0 ;m<En_mat.rows ;m++)
					{
						En_vec[m]=En_mat.get(m);
					}
					status = block->put_val<vector<number> >(config->output_section_name, name_in, En_vec);

				}
			}
		}
		status = block->put_val<double>(config->output_section_name, string("nbin_a"), num_z_bin_A);
		status = block->put_val<double>(config->output_section_name, string("nbin_b"), num_z_bin_B);
		status = block->put_val(config->output_section_name, "theta_min", config->theta_min);
    	status = block->put_val(config->output_section_name, "theta_max", config->theta_max);
		status = block->put_val(config->output_section_name, "input_section_name", config->input_section_name);

		
		disable_gsl_error_handling();
		return status;
	}
}// end of extern C



    
