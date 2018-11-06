from pyteomics import mzml,mgf
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from search_engine import MS1_extracter
import numpy as np
import copy

class Experiment():
	def __init__(self,filename,delta_precursor=1,delta_rt=1,mgf=False):

		self.filename = filename

		self.ms1 = []
		self.ms2 = []

		self.scan_to_spectrum = {}
		self.ms2_to_ms1 = {}

		self.delta_precursor = delta_precursor
		self.delta_rt = delta_rt

		self.prec_ms2 = {}
		self.rt_ms2 = {}

		if not mgf: self.read_mzml(self.filename)
		if mgf: self.read_mgf(self.filename)
		
	def __str__(self):
		ret_str = ""
		ret_str += "Number of ms1: %s\n" % (len(self.ms1))
		ret_str += "Number of ms2: %s\n" % (len(self.ms2))
		return(ret_str)
	
	def read_mgf(self,filename):
		infile = mgf.read(filename)

		for scan in infile:
			spec_obj = Spectrum(scan["m/z array"],
								scan["intensity array"],
								2,
								1,
								float(scan["params"]["rtinseconds"]),
								float(scan["params"]["rtinseconds"]),
								1,
								False,
								True)
			
			scan_id = scan["params"]["title"].split("_")[-2]

			self.scan_to_spectrum[scan_id] = spec_obj

			if 2 == 1:
				self.ms1.append(scan_id)
			if 2 == 2:
				self.ms2.append(scan_id)
				
				prec_scan_id = scan_id

				spec_obj.set_prec_mass(scan["params"]["pepmass"][0])
				prec_mass_round = int(round(spec_obj.prec_mass / self.delta_precursor))
				rt_time_round = int(round(spec_obj.scan_start_time  / self.delta_rt))

				if self.prec_ms2.has_key(prec_mass_round):
					self.prec_ms2[prec_mass_round].append(spec_obj)
				else:
					self.prec_ms2[prec_mass_round] = [spec_obj]

				if self.rt_ms2.has_key(rt_time_round):
					self.rt_ms2[rt_time_round].append(spec_obj)
				else:
					self.rt_ms2[rt_time_round] = [spec_obj]

				self.ms2_to_ms1[scan_id] = prec_scan_id

	def read_mzml(self,filename,use_prec_max_intens=False,isolation_window_bounds=0.49,filter_top=5,windowed_mode=True,window_size=50):
		infile = mzml.read(filename)
		for scan in infile:
			spec_obj = Spectrum(scan["m/z array"],
								scan["intensity array"],
								scan["ms level"],
								scan["total ion current"],
								scan["scanList"]["scan"][0]["scan start time"],
								scan["base peak intensity"],
								"positive scan" in scan.keys(),
								"negative scan" in scan.keys())

			scan_id = scan["id"].split("scan=")[1].split(" ")[0]
			
			if scan["ms level"] == 1:
				#print(list(spec_obj.mz_array))
				#print(len(spec_obj.mz_array))
				#print(scan_id)
				#raw_input("ms1")
				self.scan_to_spectrum[scan_id] = spec_obj
				self.ms1.append(scan_id)
		
		infile = mzml.read(filename)
		for scan in infile:
			spec_obj = Spectrum(scan["m/z array"],
								scan["intensity array"],
								scan["ms level"],
								scan["total ion current"],
								scan["scanList"]["scan"][0]["scan start time"],
								scan["base peak intensity"],
								"positive scan" in scan.keys(),
								"negative scan" in scan.keys())

			scan_id = scan["id"].split("scan=")[1].split(" ")[0]
			if scan["ms level"] == 2:
				if "function=" in scan["id"]: 
					prec_scan_id = scan_id
					scan_id = scan_id+"_"+scan["id"].split("function=")[1].split(" ")[0]
				else:
					self.scan_to_spectrum[scan_id] = spec_obj
					prec_scan_id = scan["precursorList"]["precursor"][0]["spectrumRef"].split("scan=")[1]
				self.ms2_to_ms1[scan_id] = prec_scan_id
				
				#print(list(spec_obj.mz_array))
				#print(len(spec_obj.mz_array))
				#print(scan_id)
				#raw_input("ms2")
				self.scan_to_spectrum[scan_id] = spec_obj
				
				try: self.scan_to_spectrum[self.ms2_to_ms1[scan_id]]
				except KeyError: continue 
				
				self.ms2.append(scan_id)

				spec_obj.filter_top_peaks(min_intensity=True)
				
				if use_prec_max_intens:
					spec_obj.set_prec_mass(scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"])
					prec_mz = self.scan_to_spectrum[scan_id].prec_mass
					lower_bound = prec_mz-isolation_window_bounds
					upper_bound = prec_mz+isolation_window_bounds
					retr_indexes = np.where(np.logical_and(self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array > lower_bound, self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array < upper_bound))[0]
					max_intensity_index = int(np.argmax(self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].intensity_array[retr_indexes]))
					isolated_precursor = self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array[retr_indexes][max_intensity_index]
					self.scan_to_spectrum[scan_id].set_prec_mass(float(isolated_precursor))
				else:
					spec_obj.set_prec_mass(scan["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]["selected ion m/z"])
					
				prec_mass_round = int(round(spec_obj.prec_mass / self.delta_precursor))
				rt_time_round = int(round(spec_obj.scan_start_time  / self.delta_rt))

				if prec_mass_round in self.prec_ms2.keys():
					self.prec_ms2[prec_mass_round].append(spec_obj)
				else:
					self.prec_ms2[prec_mass_round] = [spec_obj]

				if rt_time_round in self.rt_ms2.keys():
					self.rt_ms2[rt_time_round].append(spec_obj)
				else:
					self.rt_ms2[rt_time_round] = [spec_obj]

				
				
		#if use_prec_max_intens:
		#	for scan in self.ms2:
		#		prec_mz = self.scan_to_spectrum[scan].prec_mass
		#		lower_bound = prec_mz-isolation_window_bounds
		#		upper_bound = prec_mz+isolation_window_bounds
		#		retr_indexes = np.where(np.logical_and(self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array > lower_bound, self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array < upper_bound))[0]
		#		max_intensity_index = int(np.argmax(self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].intensity_array[retr_indexes]))
		#		isolated_precursor = self.scan_to_spectrum[self.ms2_to_ms1[scan_id]].mz_array[retr_indexes][max_intensity_index]
		#		self.scan_to_spectrum[scan].set_prec_mass(isolated_precursor)

	def get_XIC(self,mass,ppm=20,search_engine=False,positive_mode=True,negative_mode=True):
		if not search_engine:
			search_engine = MS1_extracter(ppm=ppm)

		ret_list = []
		for ms1_scan in self.ms1:
			ms1_spectrum = self.scan_to_spectrum[ms1_scan]
			if not ms1_spectrum.positive_scan == positive_mode: continue
			if not ms1_spectrum.negative_scan == negative_mode: continue
			ret_list.append([ms1_spectrum.scan_start_time,search_engine.retrieve_intensity(zip(ms1_spectrum.mz_array,ms1_spectrum.intensity_array),mass)[1]])
		
		return(ret_list)

class Spectrum():
	def __init__(self,
				 mz_array,
				 intensity_array,
				 ms_level,
				 total_ion_current,
				 scan_start_time,
				 base_peak_intensity,
				 positive_scan,
				 negative_scan,
				 prec_mass=0.0,
				 ion_injection_time=0.0,
				 top_peaks=100):
		
		self.mz_array = mz_array
		self.intensity_array = intensity_array
		self.ms_level = ms_level
		self.total_ion_current = total_ion_current
		self.scan_start_time = scan_start_time
		self.ion_injection_time = ion_injection_time
		self.prec_mass = prec_mass
		self.base_peak_intensity = base_peak_intensity
		self.positive_scan = positive_scan
		self.negative_scan = negative_scan
		
		# Test if list is sorted ... and sort if not the case. Needed later for binary search and stuff.
		if not all(a <= b for a, b in zip(self.mz_array[:-1],self.mz_array[1:])):
			self.intensity_array = [ia for ma,ia in sorted(zip(self.mz_array,self.intensity_array),key=itemgetter(0))]
			self.mz_array = sorted(self.mz_array)
		
		if not positive_scan and not negative_scan:
			positive_scan = True
			negative_scan = True

	def set_prec_mass(self,prec_mass):
		self.prec_mass = prec_mass
	
	# TODO add Minimal intensity!
	def filter_top_peaks(self,min_intensity=False,min_perc=False,windowed_mode=False,intensity_threshold=10.0,top=10,window_size=100,add_dummy_peak=False):
		"""
		Filter in multiple ways on the intensity of peaks.

	    Parameters
	    ----------
	    mz_list : list
	        The m/z values of a spectrum in a list; equal length to the intensity list
	    intensity_list : list
	        The intensity values of a spectrum in a list; equal length to the m/z list
		min_perc : bool
	        Flag to use a minimal percentage intensity to filter peaks
		windowed_mode : bool
	        Flag to use windowed mode to return the top intensity peaks
		top : int
	        The top intensity peaks to filter (in windowed mode it will return the top peaks within the window)
		window_size : int
	        The size of the window in windowed mode
		add_dummy_peak : bool
			Flag to add a dummy peak at 0.0 m/z
		
	    Returns
	    -------
		list
			the filtered m/z values from the spectrum
	    list
			the filtered intensity values from the spectrum  
	    """
		gr_intensity_list = []
		gr_mz_list = []
		
		#In the case of minimal percentage... calculate perc intensity and filter
		if min_perc:
			for i,mz in zip(self.intensity_array,self.mz_array):
				if i > min_perc:
					gr_intensity_list.append(i)
					gr_mz_list.append(mz)
		
		#In the case of windowed mode... iterate over the possible windows and intensity values; take the top per window
		if windowed_mode:
			start_index = 0
			for w in range(window_size,int(max(self.mz_array)),window_size):
				temp_mz = []
				temp_intens = []
				temp_start_index = 0
				
				#Iterate over all m/z values and see if they fall within the window
				for mz,intens in zip(self.mz_array[start_index:],self.intensity_array[start_index:]):
					if mz > w and mz <= w+window_size:
						temp_start_index += 1
						temp_mz.append(mz)
						temp_intens.append(intens)
					if mz > w+window_size:
						break
				#Next window ignore all these lower values
				start_index = temp_start_index
				
				#Use all if there are less peaks than the top number of peaks it should select
				if len(temp_mz) <= top:
					gr_mz_list.extend(temp_mz)
					gr_intensity_list.extend(temp_intens)
					continue
				
				#Get the indexes of the top peaks
				idxs = np.sort(np.argpartition(np.array(temp_intens), -top)[-top:])
				gr_mz_list.extend([temp_mz[idx] for idx in idxs])
				gr_intensity_list.extend([temp_intens[idx] for idx in idxs])
		
		if min_intensity:
			for i,mz in zip(self.intensity_array,self.mz_array):
				if i > intensity_threshold:
					gr_intensity_list.append(i)
					gr_mz_list.append(mz)

		#If not windowed, min perc or min intensity use a simple top peaks
		if not windowed_mode and not min_perc and not min_intensity:
			if len(self.intensity_array) > top:
				#Get the indexes of the top peaks
				idxs = np.sort(np.argpartition(np.array(intensity_list), -top)[-top:])
				gr_mz_list = [self.mz_array[idx] for idx in idxs]
				gr_intensity_list = [self.intensity_array[idx] for idx in idxs]
			else:
				#If there are less peaks than top peaks; return all
				gr_mz_list = self.mz_array
				gr_intensity_list = self.intensity_array
		
		#If needed add a dummy peak; this is important later since I want to take into account immonium ions and small fragments
		if add_dummy_peak:
			gr_mz_list.insert(0,0.0)
			gr_intensity_list.insert(0,1.0)
		
		self.mz_array = gr_mz_list
		self.intensity_array = gr_intensity_list


"""
	def filter_top_peaks(self,top=1000):
		if len(self.intensity_array) > top:
			#Get the indexes of the top peaks
			idxs = np.sort(np.argpartition(np.array(self.intensity_array), -top)[-top:])
			gr_mz_list = [self.mz_array[idx] for idx in idxs]
			gr_intensity_list = [self.intensity_array[idx] for idx in idxs]
		else:
			#If there are less peaks than top peaks; return all
			gr_mz_list = self.mz_array
			gr_intensity_list = self.intensity_array
		
		self.mz_array = gr_mz_list
		self.intensity_array = gr_intensity_list
"""

if __name__ == "__main__":
	exp = Experiment("23Mar17_HaCaT_ox15h_1.mzML")
	for scan_num in exp.ms2:
		print(exp.scan_to_spectrum[scan_num].prec_mass)