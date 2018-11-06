import logging
from operator import itemgetter
import pandas as pd
from scipy.stats import hypergeom
import numpy as np

class ScoringFunction():
	def __init__(self):
		pass

	def calc_hypergeom(self,
			aligned_peaks,
			theoretical_spectrum,
			experimental_spectrum,
			num_bins,
			log10_transform=True):

		num_aligned_peaks = len(aligned_peaks)
		num_exp_peaks = len(experimental_spectrum)
		num_theor_peaks = len(theoretical_spectrum)
		hyper = hypergeom(num_bins,num_theor_peaks,num_exp_peaks)

		if log10_transform: return(np.log10(hyper.pmf(np.arange(0,num_aligned_peaks+1))[-1])*-1)
		else: return(hyper.pmf(np.arange(0,num_aligned_peaks+1))[-1])

	def calc_intensity_explained(self,
			aligned_peaks,
			theoretical_spectrum,
			experimental_spectrum):

		explained_intensity = 0.0
		already_assigned = []
		for ap in aligned_peaks:
			if ap[1][0] in already_assigned:
				continue
			explained_intensity += ap[1][1]
			already_assigned.append(ap[1][0])
		total_intensity_exp = sum(experimental_spectrum)

		return((explained_intensity/total_intensity_exp)*100)

	def ms1_error_perc(self,ppm_error):
		pass

	def ms2_error_perc(self,ppm_error):
		pass