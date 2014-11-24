from cosmosis import option_section, names
import wigglez_pk
import os

DEFAULT_KMAX=0.3
SCRIPT_DIRNAME=os.path.split(__file__)[0]
DEFAULT_DATA_DIRNAME=os.path.join(SCRIPT_DIRNAME, "data")



def setup(options):
	raise ValueError("The WiggleZ code is not yet complete due to the NL power handling")
	data_dir = options.get_string(option_section, "data_dir", default=DEFAULT_DATA_DIRNAME)
	root_dir = options.get_string(option_section, "root_dir", default=SCRIPT_DIRNAME)
	common = os.path.join(data_dir, "wigglez_nov11_common.dataset")
	datasets = []
	for redshift_bin in 'abcd':
		data_file = "%s/wigglez_nov11%s.dataset" % (data_dir,redshift_bin)
		dataset = wigglez_pk.WigglezDataset(common, data_file, root_dir)
		datasets.append(dataset)
	return datasets

def execute(block, config):
	datasets = config
	z1 = block[names.distances, "z"]
	d_a = block[names.distances, "d_a"]
	h = block[names.distances, "h"]
	d_v = ((1+z1)**2 * d_a**2 * z1 / h)**(1./3)

	if z1[1]<z1[0]:
		z1 = z1[::-1]
		d_v = d_v[::-1]

	k = block["matter_power_gal", "k_h"]
	z = block["matter_power_gal", "z"]
	p = block["matter_power_gal", "p_k"]

	H0 = block[names.cosmological_parameters, "h0"]*100

	L = 0
	for dataset in datasets:
		like = dataset.likelihood(H0, z1, d_v, k, z, p)
		block[names.likelihoods, dataset.name+"_like"] = like
		L += like

	block[names.likelihoods, "wigglez_pk_like"] = L
	return 0


def cleanup(config):
	pass