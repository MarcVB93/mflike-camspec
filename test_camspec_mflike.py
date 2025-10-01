import os
import tempfile
import unittest
import sys
import numpy as np
import matplotlib.pyplot as plt

packages_path = os.environ.get("COBAYA_PACKAGES_PATH") or os.path.join(tempfile.gettempdir(), "camspecmflike")

cosmo_params = {
	# PR3_12.6 TTTEEE best fit parameters.
	# [arXiv:2205.10869]

    "ombh2" : 0.02225,
    "omch2" : 0.1197,
    "cosmomc_theta" : 0.0104102,
    "tau" : 0.0533,
    "As" : 2.09262e-9,
    "ns" : 0.9667
}

nuisance_params = {
	# Planck 2018 best fit parameters.
	# [arXiv:2205.10869]

	# TT parameters
	"amp_100" : 0.0,
	"amp_143" : 17.2,
	"amp_217" : 12.8,
	"amp_143x217" : 8.0,
	"n_100" : 1.0,
	"n_143" : 1.12,
	"n_217" : 1.40,
	"n_143x217" : 1.53,
	
	
	# calibration parameters
	"calTE" : 1.0005,
	"calEE" : 1.0013,
	"A_planck" : 1.0005
}

# Calculated using the original CamSpec NPIPE likelihood with the parameters above.
chi2_camspec_npipe = 10573.39061

class CamspecMFLikeTest(unittest.TestCase):
	#def setUp(self):
	#	from cobaya.install import install
		
	#	install({
	#		"likelihood" : {
	#			"camspecmflike.CamspecMFLike" : None
	#		}
	#	}, path = packages_path, skip_global = True)
	
	def test_camspecmflike(self):
		import camb
		
		camb_params = cosmo_params.copy()
		camb_params.update({ "lmax" : 2500, "lens_potential_accuracy" : 1 })
		
		pars = camb.set_params(**camb_params)
		results = camb.get_results(pars)
		powers = results.get_cmb_power_spectra(pars, CMB_unit = "muK")
		cl_dict = { k : powers["total"][2:2501, v] for k, v in { "tt" : 0, "te" : 3, "ee" : 1 }.items() }
		cl_dict["ell"] = np.arange(2, 2501)
		
		from camspecmflike import CamspecMFLike
		
		camspec = CamspecMFLike({
			"packages_path" : packages_path,
			"data_folder" : "data/",
			"covfile" : "like_NPIPE_12.6_unified_cov.bin",
			"specfile" : "like_NPIPE_12.6_unified_spectra.txt",
			"dataranges" : "like_NPIPE_12.6_unified_data_ranges.txt"
		})
		
		chi2 = -2.0 * camspec.loglike(cl_dict, **nuisance_params)

		print("-----------------------------------------------------------------")
		print("Calculated chi square value          : {:10.5f}".format(chi2))
		print("Chi square value from CamSpec NPIPE  : {:10.5f}".format(chi2_camspec_npipe))
		print("-----------------------------------------------------------------")
		

if __name__ == "__main__":
	unittest.main()
