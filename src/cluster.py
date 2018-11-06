"""
 Copyright (c) 2017 Robbin Bouwmeester
 Permission is hereby granted, free of charge, to any person
 obtaining a copy of this software and associated documentation
 files (the "Software"), to deal in the Software without
 restriction, including without limitation the rights to use,
 copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the
 Software is furnished to do so, subject to the following
 conditions:

 The above copyright notice and this permission notice shall be
 included in all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE."""

__author__ = "Robbin Bouwmeester"
__copyright__ = "Copyright 2017"
__credits__ = ["Robbin Bouwmeester"]
__license__ = "MIT"
__version__ = "0.1"
__maintainer__ = "Robbin Bouwmeester"
__email__ = "Robbin.bouwmeester@ugent.be"
__status__ = "nightly funzies"

import pandas as pd
from itertools import groupby
import logging
from operator import itemgetter
import re
from matplotlib import pyplot as plt

class Clusters():
	def __init__(self,ppm=5,max_dist=0.15,rt_dist=3.0):
		self.ppm = ppm
		self.max_dist = max_dist
		self.clusters = []
		self.tot_clust = 0
		self.tot_added = 0
		self.rt_dist = rt_dist
		
	def add_to_clusters(self,prec_mass,mz_list,intensity_list,rt=0.0):
		self.tot_added += 1
		#print("----")
		if self.tot_clust == 0:
			self.tot_clust += 1
			name = "Cluster_%s" % (self.tot_clust)
			self.clusters.append(Cluster(prec_mass,mz_list,intensity_list,name,rt=rt))
			return()
			
		mass_error_threshold = (prec_mass*self.ppm)/1000000
		
		
		need_for_new_clust = True
		for c in self.clusters:
			#print(mass_error_threshold)
			#print(abs(c.prec_mass-prec_mass))
			if abs(c.prec_mass-prec_mass) < mass_error_threshold and abs(c.rt-rt) < self.rt_dist:
				if not need_for_new_clust: continue
				score,aligned_list = self.align_spectra(zip(c.mz_list,c.intensity_list),zip(mz_list,intensity_list))
				
				#raw_input()
				if score < self.max_dist:
					#print("----")
					#print(self.tot_added)
					#print(len(self.clusters))
					#print(c.num_in_clust)
					#print("score: %s" % (score))
					need_for_new_clust = False
					c.add_spectrum(aligned_list,prec_mass,rt=rt)
		
		if need_for_new_clust:
			self.tot_clust += 1
			name = "Cluster_%s" % (self.tot_clust)
			self.clusters.append(Cluster(prec_mass,mz_list,intensity_list,name,rt=rt))
		return()
					
	
	def normalize_spectrum(self,spectrum):
		tot_intens = sum([intens for peak,intens in spectrum])
		return(zip([peak for peak,intens in spectrum],[intens/tot_intens for peak,intens in spectrum]))
	
	def align_spectra(self,spectrum_one,spectrum_two):
		spectrum_one = self.normalize_spectrum(spectrum_one)
		spectrum_two = self.normalize_spectrum(spectrum_two)
		
		max_dist = sum([intensity_one for mz_one,intensity_one in spectrum_one]) + sum([intensity_two for mz_two,intensity_two in spectrum_two])
		tot_dist = 0.0
		aligned_spectrum_two = 0.0
		
		aligned_list = []
		
		#print(len(spectrum_two))

		for peak_one,intensity_one in spectrum_one:
			peak_has_aligned = False
			index_aligned = -1
			for peak_two,intensity_two in spectrum_two:
				if peak_has_aligned: continue
				index_aligned += 1
				if abs(peak_one-peak_two) < 0.005:
					aligned_list.append([[peak_one,intensity_one],[peak_two,intensity_two]])
					dist_intensity = abs(intensity_one-intensity_two)
					aligned_spectrum_two += intensity_two
					peak_has_aligned = True
			if not peak_has_aligned:
				aligned_list.append([[peak_one,intensity_one],[0.0,0.0]])
				tot_dist += intensity_one
				
			#peak_aligned,error = self.peak_to_spec(spectrum_two,peak_one)
			else:
				tot_dist += dist_intensity
				del spectrum_two[index_aligned]
		
		for peak_two,intensity_two in spectrum_two:
			aligned_list.append([[0.0,0.0],[peak_two,intensity_two]])
			tot_dist += intensity_two
		
		#print(sum([intensity_one for mz_one,intensity_one in spectrum_one]))
		#print(sum([intensity_two for mz_two,intensity_two in spectrum_two]))
		
		#print([mz_one for mz_one,intensity_one in spectrum_one])
		#print([mz_two for mz_two,intensity_two in spectrum_two])
		#print(sum([intensity_one for mz_one,intensity_one in spectrum_one]))
		#print(sum([intensity_two for mz_two,intensity_two in spectrum_two]))

		return((tot_dist/max_dist),aligned_list)
		
	def get_clusters(self,n_min=2):
		return([c for c in self.clusters if c.num_in_clust > n_min])
			
	def peak_to_spec(self,spectrum,peak,threshold=False,mz_loc=1,orig_spec=[],pos_in_orig_spec=0,return_all_peaks=False):
		def _look_around_matched_peak(spectrum,peak,pos,rev=False):
			error = 0.0
			ret_list = []
			while error < threshold and pos < len(spectrum)-1 and pos > 1:
				if rev: pos = pos-1
				else: pos = pos+1
				error = abs(spectrum[pos][mz_loc]-peak)
				if error < threshold:
					ret_list.append(spectrum[pos])
			if len(ret_list) > 0:
				print("no return list")
			return(ret_list)

		if not threshold:
			threshold = (peak*self.ppm)/1000000		

		curr_sel = int(len(spectrum)/2)
		if return_all_peaks:
			if orig_spec == []:
				pos_in_orig_spec = curr_sel
				orig_spec = spectrum

		if abs(spectrum[curr_sel][mz_loc] - peak) < threshold:
			_look_around_matched_peak(orig_spec,peak,pos_in_orig_spec)
			return(spectrum[curr_sel],abs(spectrum[curr_sel][mz_loc] - peak))
		elif curr_sel == 0 or curr_sel == int(len(spectrum)):
			return(False,False)
		elif (spectrum[curr_sel][mz_loc] - peak) > 0.0:
			if return_all_peaks:
				if curr_sel % 2 != 0: pos_in_orig_spec = pos_in_orig_spec-int(len(spectrum[0:curr_sel])/2)-1
				else: pos_in_orig_spec = pos_in_orig_spec-int(len(spectrum[0:curr_sel])/2)
			return(self.peak_to_spec(spectrum[0:curr_sel],peak,mz_loc=mz_loc,pos_in_orig_spec=pos_in_orig_spec,orig_spec=orig_spec))
		else:
			if return_all_peaks:
				if curr_sel % 2 != 0: pos_in_orig_spec = pos_in_orig_spec+int(len(spectrum[curr_sel:])/2)
				else: pos_in_orig_spec = pos_in_orig_spec+int(len(spectrum[curr_sel:])/2)
			return(self.peak_to_spec(spectrum[curr_sel:],peak,mz_loc=mz_loc,pos_in_orig_spec=pos_in_orig_spec,orig_spec=orig_spec))
					
					
class Cluster():
	def __init__(self,prec_mass,mz_list,intensity_list,name,rt=0.0,num_in_clust=1):
		self.prec_mass = prec_mass
		self.mz_list = mz_list
		self.intensity_list = self.normalize_spectrum(intensity_list)
		self.name = name
		self.num_in_clust = num_in_clust
		self.rt = rt
		
	def normalize_spectrum(self,mz_list):
		return([peak/sum(mz_list) for peak in mz_list])
		
	def remove_low_intensity(self,min_intensity=0.005):
		new_mz = []
		new_intens = []
		for mz,intens in zip(self.mz_list,self.intensity_list):
			if intens > min_intensity:
				new_mz.append(mz)
				new_intens.append(intens)
		self.mz_list = new_mz
		self.intensity_list = self.normalize_spectrum(new_intens)
		
		
	def add_spectrum(self,aligned_list,prec_mass_aligned,rt=0.0):
		new_mz_list = []
		new_intensity_list = []
		
		#print("======")
		#print(self.mz_list)
		#print(self.intensity_list)
		#print(aligned_list)
		
		
		for spec_clust,spec_new in aligned_list:
			peak_clust,intensity_clust = spec_clust
			peak_new,intensity_new = spec_new
			if peak_new == 0.0:
				new_mz_list.append(peak_clust)
				#print("peak_new: %s" % (peak_clust))
				new_intensity_list.append(intensity_clust*((self.num_in_clust-1)/float(self.num_in_clust)))
			elif peak_clust == 0.0:
				new_mz_list.append(peak_new)
				new_intensity_list.append(intensity_clust+((intensity_new-intensity_clust)*(1/float(self.num_in_clust))))
				#print("intens_new2: %s" % (intensity_clust))
				#print("intens_new2: %s" % (intensity_new))
				#print("peak_new2: %s" % (peak_new))
			else:
				#print("peak_new3: %s" % (peak_clust+((peak_new-peak_clust)*(1/float(self.num_in_clust)))))
				new_mz_list.append(peak_clust+((peak_new-peak_clust)*(1/float(self.num_in_clust))))
				new_intensity_list.append(intensity_clust+((intensity_new-intensity_clust)*(1/float(self.num_in_clust))))
		self.mz_list = new_mz_list
		self.intensity_list = self.normalize_spectrum(new_intensity_list)
		self.prec_mass = self.prec_mass+((prec_mass_aligned-self.prec_mass)*(1/float(self.num_in_clust)))
		self.rt = self.rt+((rt-self.rt)*(1/float(self.num_in_clust)))
		
		#print(self.mz_list)
		#print(self.intensity_list)

		#if self.num_in_clust % 100 == 0:
		#	plt.scatter(self.mz_list,self.intensity_list,s=0.1)
		#	for mz,inten in zip(self.mz_list,self.intensity_list):
		#		if inten > 0.01:
		#			print(mz,inten)
		#		plt.axvline(x=mz,ymax=inten)
		#	plt.show()
		#	plt.close()
		
		self.remove_low_intensity()
		self.num_in_clust += 1
		
		
		
	
		
		
	
