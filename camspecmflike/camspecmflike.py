from cobaya.likelihoods.base_classes import InstallableLikelihood
import numpy as np
from scipy import linalg
import sys
from .fgspectra import cross as fgc
from .fgspectra import power as fgp
from .fgspectra import frequency as fgf

class CamspecMFLike(InstallableLikelihood):
	def initialize(self):
		self.expected_params = [
			# TT parameters
			'amp_100', 'amp_143', 'amp_217', 'amp_143x217', 'n_100', 'n_143', 'n_217', 'n_143x217',

			# calibration parameters
			'calTE', 'calEE', 'A_planck',
		]
		
		self.enable_tt = True
		self.enable_te = True
		self.enable_ee = True
		
		self.prepare_data()

	
	def prepare_data(self):
		
		self.load_ranges()
		
		self.nbintt = [ b - a + 1 for a, b in zip(self.nmin[0:-2], self.nmax[0:-2]) ]
		self.nbinte = [ self.nmax[-2] - self.nmin[-2] + 1 ]
		self.nbinee = [ self.nmax[-1] - self.nmin[-1] + 1 ]

		self.ellmin = np.array([30,   500,  500,  30,   30])
		self.ellmax = np.array([2000, 2500, 2500, 2000, 2000])
		
		self.crosstt = [ (0,0), (1,1), (0,1) ]
		self.crosste = [ (0,0) ]
		self.crossee = [ (0,0) ]
		
		self.freqs = [ 143, 217 ]
		
		self.sys_vec = None
		self.inv_cov = None
		
		self.b_ell = None
		self.b_dat = None
		self.win_func = None
		

		self.log.debug('Loading spectrum data.')
		self.b_dat = self.load_spectra()

		self.log.debug('Loading inv cov.')
		self.covmat = self.load_covmat()

		self.inv_cov = np.linalg.inv(self.covmat)
		
		
		self.log.debug('Done preparing all data!')


	def load_ranges(self):
		
		# NEW FUNCTION FOR CAMSPEC_MFLIKE!
		# --------------------------------
		# This function loads the data ranges that will be used

		self.nmin = np.loadtxt(self.data_folder + self.dataranges, dtype = int, usecols=range(1,3))[1:, 0]
		self.nmax = np.loadtxt(self.data_folder + self.dataranges, dtype = int, usecols=range(1,3))[1:, 1]
		

	def load_spectra(self):
		
		# NEW FUNCTION FOR CAMSPEC_MFLIKE!
		# --------------------------------
		# This function loads the power spectra data and puts
		# it in the same format that was used in PlikMFLike

		spec = np.loadtxt(self.data_folder + self.specfile, dtype = float, usecols=range(1,6))

		b_dat = []

		for i in range(0, self.nspec):
    			b_dat += list(spec[self.nmin[i]:self.nmax[i]+1, i])
		
		return b_dat
	
	
	def load_covmat(self):
		
		# NEW FUNCTION FOR CAMSPEC_MFLIKE!
		# --------------------------------
		# This function loads the COVMAT data and puts
		# it in the same format that was used in PlikMFLike
		
		with open((self.data_folder + self.covfile), "rb") as cov_f:
			self.fullcov = np.fromfile(cov_f, dtype=[np.float32, np.float64]["64.bin" in self.covfile])
		
		self.fullcov.shape = (int(np.sqrt(len(self.fullcov))), int(np.sqrt(len(self.fullcov))))

		self.fullcov = self.fullcov[(len(self.fullcov) - np.sum((self.ellmax+1)-self.ellmin)):, (len(self.fullcov) - np.sum((self.ellmax+1)-self.ellmin)):]

		if (self.nmin == self.ellmin).all() and (self.nmax == self.ellmax).all():
			return self.fullcov

		else:
			new_cov = self.trim_covmat()
			return new_cov


	def trim_covmat(self):
		
		# NEW FUNCTION FOR CAMSPEC_MFLIKE!
		# --------------------------------
		# This function trims the COVMAT according to
		# the ell ranges specified for the run

		offsets = []
		current_offset = 0

		for i in range(self.nspec):
		    N = self.ellmax[i] - self.ellmin[i] + 1
		    offsets.append(current_offset)
		    current_offset += N

		offsets = np.array(offsets)
		selected_indices = []

		for i in range(self.nspec):
		    ell_start = self.ellmin[i]
		    ell_stop = self.ellmax[i]
		    block_offset = offsets[i]

		    ells = np.arange(ell_start, ell_stop + 1)

		    mask = (ells >= self.nmin[i]) & (ells <= self.nmax[i])
		    selected_ell_indices = np.where(mask)[0]

		    global_indices = selected_ell_indices + block_offset
		    selected_indices.append(global_indices)

		selected_indices = np.concatenate(selected_indices)

		new_covmat = self.fullcov[np.ix_(selected_indices, selected_indices)]

		return new_covmat
	
	
	def get_model(self, cl, **params_values):
		self.log.debug('Start calculating model.')
		l0 = int(2 - cl['ell'][0])
		ls = cl['ell'][l0:self.shape+l0]
		dl_tt = cl['tt'][l0:self.shape+l0]
		dl_te = cl['te'][l0:self.shape+l0]
		dl_ee = cl['ee'][l0:self.shape+l0]
		
		fg_tt = self.get_foregrounds(params_values, ls)
		fg_te = np.zeros((self.nspecte, self.shape))
		fg_ee = np.zeros((self.nspecee, self.shape))
		
		self.log.debug('Summing theory = CMB + foreground.')
		
		x_theory = np.zeros((self.nspec, self.shape))
		x_theory[0                         : self.nspectt                          ,:self.shape] = np.tile(dl_tt, (self.nspectt, 1)) + fg_tt
		x_theory[self.nspectt              : self.nspectt+self.nspecte             ,:self.shape] = np.tile(dl_te, (self.nspecte, 1)) + fg_te
		x_theory[self.nspectt+self.nspecte : self.nspectt+self.nspecte+self.nspecee,:self.shape] = np.tile(dl_ee, (self.nspecee, 1)) + fg_ee
		
		self.log.debug('Completed theory vector. Calibrating')
		
		x_theory[self.nspectt              : self.nspectt+self.nspecte             ,:self.shape] /= params_values['calTE']
		x_theory[self.nspectt+self.nspecte : self.nspectt+self.nspecte+self.nspecee,:self.shape] /= params_values['calEE']
		x_theory /= (params_values['A_planck'] ** 2.0)
		
		x_model = []
		
		for i in range(0, self.nspec):
    			x_model += list(x_theory[i, self.nmin[i]-2:self.nmax[i]-1])
			
		self.log.debug('Done calculating model.')
		
		return x_model
	
	def logp(self, **params_values):
		cl = self.provider.get_Cl(ell_factor = True)
		return self.loglike(cl, **params_values)
	
	def loglike(self, cl, **params_values):
		x_model = self.get_model(cl, **params_values)
		diff_vec = np.array(self.b_dat) - np.array(x_model)
		tmp = self.inv_cov @ diff_vec
		return -0.5 * np.dot(tmp, diff_vec)
	
	@property
	def use_tt(self):
		return self.enable_tt
	
	@use_tt.setter
	def use_tt(self, val):
		self.enable_tt = val
	
	@property
	def use_te(self):
		return self.enable_te
	
	@use_te.setter
	def use_te(self, val):
		self.enable_te = val
	
	@property
	def use_ee(self):
		return self.enable_ee
	
	@use_ee.setter
	def use_ee(self, val):
		self.enable_ee = val
	
	@property
	def frequencies(self):
		return self.freqs
	
	@property
	def nspectt(self):
		return len(self.nbintt)
	
	@property
	def nspecte(self):
		return len(self.nbinte)
	
	@property
	def nspecee(self):
		return len(self.nbinee)
	
	@property
	def nbin(self):
		# total number of bins
		return sum(self.nbintt) + sum(self.nbinte) + sum(self.nbinee)
	
	@property
	def nspec(self):
		# total number of spectra
		return self.nspectt + self.nspecte + self.nspecee
	
	@property
	def shape(self):
		return self.lmax_win-1
	


	def get_foregrounds(self, fg_params, ell, powerlaw_pivot=1500):
		
		self.fg_TT_amp = np.zeros(3)
		self.fg_TT_amp[0] = fg_params['amp_143']
		self.fg_TT_amp[1] = fg_params['amp_217']
		self.fg_TT_amp[2] = fg_params['amp_143x217']
		
		self.fg_TT_n = np.zeros(3)
		self.fg_TT_n[0] = fg_params['n_143']
		self.fg_TT_n[1] = fg_params['n_217']
		self.fg_TT_n[2] = fg_params['n_143x217']
		
		self.foregrounds = np.array([self.fg_TT_amp[ii] * (ell / powerlaw_pivot) ** self.fg_TT_n[ii] for ii in range(len(self.fg_TT_amp))])  # Powerlaw Foregrounds
		
		return self.foregrounds
