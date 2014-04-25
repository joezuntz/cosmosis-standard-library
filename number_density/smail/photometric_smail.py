#Code by Donnacha Kirk
#Needs updating!  Some things should move to setup function

import numpy as np
import cosmosis_py

def gaussian(z,mu,sigma):
	return np.exp(-0.5*(z-mu)**2/sigma**2) / np.sqrt(2*np.pi) / sigma

def smail_distribution(z, alpha, beta, z0):
	return (z**alpha) * np.exp(-(z/z0)**beta)

def photometric_error(z, Nz, sigma_z, bias):
	nz = len(z)
	output = np.zeros((nz,nz))
	for i in xrange(nz):
		p = gaussian(z,z[i]-bias,sigma_z*(1+z[i]))
		#Could add f_cat or other terms here.
		output[:,i] = p * Nz[i]
	return output

def find_bins(z, nz_true, nbin):
	nz_true = nz_true/nz_true.sum()*nbin
	cum = np.cumsum(nz_true)
	bin_edges = [0.0]
	for i in xrange(1,nbin):
		edge = np.interp(1.0*i,cum,z)
		bin_edges.append(edge)
	bin_edges.append(z.max())	
	return np.array(bin_edges)

def compute_bin_nz(z_prob_matrix, z, edges, ngal):
	NI = []
	nbin = len(edges)-1
	dz = z[1]-z[0]
	for low,high in zip(edges[:-1], edges[1:]):
		w = np.where((z>low) & (z<high))[0]
		ni = z_prob_matrix[w,:].sum(axis=0)
		#make this integrate to ngal/nbin, since we
		#have divided the ngal total up evenly
		ni/=(ni.sum() * dz) 
		ni*=ngal/(nbin*ni.sum()*dz)
		assert(len(ni)==len(z))
		NI.append(ni)
	return NI
	

def compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias):
	#Set up Smail distribution of z vector
	nz_true = smail_distribution(z, alpha, beta, z0)
	z_prob_matrix = photometric_error(z, nz_true, sigma_z, bias)
	edges = find_bins(z,nz_true,nbin)
	bin_nz = compute_bin_nz(z_prob_matrix, z, edges, ngal)
	return edges,bin_nz
		
def execute(fits_handle):
	#read fits file data
	try:
		data = pydesglue.DesDataPackage.from_fits_handle(fits_handle)

		#get parameters from fits data
		section = pydesglue.section_names.number_density_params
		dz      = data.get_param(section,"DZ")
		zmax    = data.get_param(section,"ZMAX")
		alpha   = data.get_param(section,"ALPHA")
		beta    = data.get_param(section,"BETA")
		z0      = data.get_param(section,"Z0")
		sigma_z = data.get_param(section,"SIGZ")
		nbin    = data.get_param(section,"NBIN")
		ngal    = data.get_param(section,"NGAL",default=20.0)
		bias    = data.get_param(section,"BIAS",default=0.0)
		
		#Compute the redshift vector
		z = np.arange(0,zmax+dz/2,dz)
		
		#Run the main code for getting n(z) in bins
		edges,bins = compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias)

		#Save them back to memory
		nz_section = pydesglue.section_names.wl_number_density
		#Save the basic parameters
		data.set_param(nz_section,"NBIN",nbin)
		data.set_param(nz_section,"NZ",len(z))
		data.set_data(nz_section, "Z", z)
		#Loop through the bins
		for i,bin in enumerate(bins):
			#The bin numbering starts at 1
			b=i+1
			name = "BIN_%d" % b
			#Save the bin edges as parameters
			data.set_param(nz_section,"EDGE_%d"%b,edges[i])
			#And save the bin n(z) as a column
			data.set_data(nz_section, name, bin)
		#Also save the upper limit to the top bin
		data.set_param(nz_section,"EDGE_%d"%(nbin+1),edges[-1])
		#Finally save everything back to the fits file.
		data.write_to_fits_handle(fits_handle)
	except KeyboardInterrupt:
		raise KeyboardInterrupt
	except Exception as error:
		print "Failed to generate spectra in kirk_nz"
		print error
		return 1
	return 0
		
	

# function [n_z_tot ngals_bin pzpzs]=get_nofz_amarar06(s,z_vals,verb)
# % n_z=z_vals.^alpha .* exp(-(z_vals/z_0).^beta);
# % This convolves the n(z) distributions from get_nofz_equal_tophat
# % with a Gaussian corresonding to s.deltaz
# % NB. convolving n(z) is not quite the right thing to do, since it
# % does not conserve the total number of galaxies at a given z
# % In fact, should do some backwards thing. (Probably index kernel
# % opposite?)
# % Then it adds an outlier fraction.
# %
# % based on get_nofz_equal_th_photoz 
# % SLB 27 Feb 2007
# 
# % set defaults
# if (~exist('verb')) verb=0; end % verbosity level
# if (~isfield(s,'verb')) s.verb=0; end
# verb=max([s.verb verb]);
# if (~isfield(s,'ng')) s.ng=1; end % arbitrary normalisation if not reqd
# 
# % make overall Smail et al n(z)
# % n_z=z_vals.^alpha .* exp(-(z_vals/z_0).^beta);
# n_z=z_vals.^s.alpha .* exp(-(z_vals/s.z_0).^s.beta);
# 
# % make P(z_phot|z_spec) on a 2d grid
# nz=length(z_vals);
# dz=z_vals(2)-z_vals(1); % should really check they are evenly spaced
# pzpzs=zeros(nz,nz);
# vbx=3;
# for iz=1:nz
#     z_t=z_vals(iz);
#     if (verb>=vbx) clf; end
#     
#     % P_stat
#     z_m=z_t;
#     sigma_z=s.deltaz * (1+z_m); % red-shift errors are deltaz*(1+z)
#     p_stat=(1/(sqrt(2*pi)*sigma_z))*exp(-(z_vals-z_m).^2 / (2*sigma_z^2));    
#     %z_m_test=sum(z_vals.*p_stat)*dz, % assumes equally spaced z values
#     %sqrt(sum(z_vals.^2.*p_stat)*dz- z_m_test^2), % I assume A2 should be this!?
#     if (verb>=vbx) plot(z_vals,p_stat); end
#     
#     % P_cat_minus
#     z_cat=z_t-s.Deltaz;
#     sigma_z=s.deltaz * (1+z_cat); % red-shift errors are deltaz*(1+z)
#     p_cat_minus=(1/(sqrt(2*pi)*sigma_z))*exp(-(z_vals-z_cat).^2 / (2*sigma_z^2));    
#     %z_m_test=sum(z_vals.*p_cat_minus)*dz
#     %sqrt(sum(z_vals.^2.*p_cat_minus)*dz-z_m_test^2)
#     if (verb>=vbx)
#         hold on
#         plot(z_vals,p_cat_minus,'r-')
#     end
#     
#     % P_cat_plus
#     z_cat=z_t+s.Deltaz;
#     sigma_z=s.deltaz * (1+z_cat); % red-shift errors are deltaz*(1+z)
#     p_cat_plus=(1/(sqrt(2*pi)*sigma_z))*exp(-(z_vals-z_cat).^2 / (2*sigma_z^2));    
#     %z_m_test=sum(z_vals.*p_cat_plus)*dz
#     %sqrt(sum(z_vals.^2.*p_cat_plus)*dz-z_m_test^2)
#     if (verb>=vbx) plot(z_vals,p_cat_plus,'g-'); end
#     
#     % add it all up
#     p_tot=(1-s.fcat)*p_stat + s.fcat/2 * (p_cat_minus+p_cat_plus); % ?? must surely be a factor of 2 here?
#     % normalise, in case e.g. p_cat_minus is outside range ????
#     p_tot=p_tot/(sum(p_tot)*dz);
#     if (verb>=vbx)
#         plot(z_vals,p_tot,'k:')
#         axis([0 3 0 6])
#         pause(0.0001)
#     end
#     
#     % get P(z_t,z_phot)=P(z_phot|z_t)*P(z_t)
#     %    pzpzs(:,iz)=p_tot;
#     pzpzs(:,iz)=p_tot*n_z(iz);
# 
# end
# 
# if (verb>=2)
#     imagesc(pzpzs); colorbar; set(gca,'ydir','normal')
# end
# 
# 
# %% Now bin it
# 
# % find bin divisions
# % normalise in a way that makes it easier to find z bin divisions
# n_z_tmp = s.nbin* n_z/sum(n_z);
# cn_z=cumsum(n_z_tmp);
# % plot(z_vals,cn_z)
# eps=1e-5;
# z_cuts=interp1(cn_z,z_vals,1:(s.nbin-1));
# z_cuts_l = [0 z_cuts];
# z_cuts_u = [z_cuts max(z_vals)];
# 
# % cut up the z distribution
# for ibin=1:s.nbin
#     iz=find((z_vals>z_cuts_l(ibin))&(z_vals<z_cuts_u(ibin)));
#     n_z_tot(:,ibin)=sum(pzpzs(iz,:),1);
#     ngals_bin(ibin)=sum(n_z_tot(:,ibin));
#     % normalise wrt chi later
# end
# ngals_bin=ngals_bin*s.ng/sum(ngals_bin);
# 
# % check everything
# if (verb>=1)
#     %clf
#     for ibin=1:s.nbin
#         if (mod(ibin,2)==1) col='r-'; else col='b--'; end
#         plot(z_vals,n_z_tot(:,ibin),col)
#         hold on
#     end
# end
# 
# % nomalise
# for ibin=1:s.nbin
#     norm_nz=sum(n_z_tot(:,ibin));
#     n_z_tot(:,ibin)=n_z_tot(:,ibin)/norm_nz;
# end
# 
# %debugging
# ibin = ibin;