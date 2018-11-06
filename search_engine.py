import logging
from operator import itemgetter
import pandas as pd
from joblib import Parallel, delayed

class PrecursorFilter():
	def __init__(self,db,ppm=None,dalton=None):
		self.db = db
		self.ppm = ppm
		self.dalton = dalton
		if not ppm and not dalton: self.ppm = 20

	def retrieve_entry_pre_c_mass(self,pre_c_mass,lower_limit_mz=0.0,upper_limit_mz=10000.0):
		if self.ppm: mass_error_threshold = (pre_c_mass*self.ppm)/1000000
		else: mass_error_threshold = self.dalton

		if pre_c_mass <= (lower_limit_mz-mass_error_threshold) or pre_c_mass >= (upper_limit_mz+mass_error_threshold):
			return([])

		ret_entries = []

		loc_dict = int(pre_c_mass) - int(pre_c_mass) % self.db.dalt_diff_lookup_bin
		loc_dict_lower = (int(pre_c_mass-mass_error_threshold)) - (int(pre_c_mass-mass_error_threshold)) % self.db.dalt_diff_lookup_bin
		loc_dict_upper = (int(pre_c_mass+mass_error_threshold)) - (int(pre_c_mass+mass_error_threshold)) % self.db.dalt_diff_lookup_bin

		locs_to_search = list(set([loc_dict,loc_dict_lower,loc_dict_upper]))
		for loc in locs_to_search:
			try:
				for name,entr in self.db.ms1_dict_lookup[loc].items():
					mass_error = abs(entr.mw-pre_c_mass)
					if mass_error < mass_error_threshold:
						ret_entries.append([name,mass_error,entr])
			except KeyError:
				logging.warning("Could not find an entry in the database for prec mass: %s" % (pre_c_mass))
				continue
		return(ret_entries)

class MS1_extracter():
	def __init__(self,ppm=10):
		self.ppm = ppm

	def retrieve_intensity(self,spectrum,peak,threshold=False,mz_loc=0):
		if not threshold:
			threshold = (peak*self.ppm)/1000000
		spectrum = list(spectrum)
		curr_sel = int(len(spectrum)/2)

		if abs(spectrum[curr_sel][mz_loc] - peak) < threshold:
			return(spectrum[curr_sel])
		elif curr_sel == 0 or curr_sel == int(len(spectrum)):
			return([0.0,0.0])
		elif (spectrum[curr_sel][mz_loc] - peak) > 0.0:
			return(self.retrieve_intensity(spectrum[0:curr_sel],peak))
		else:
			return(self.retrieve_intensity(spectrum[curr_sel:],peak))

class IsotopeAlignment():
	def __init__(self,ppm=5):
		self.ppm = ppm

	def find_all_isotopics_dists(self,ms1_scan,peak_dist=1.0033548378):
		ms1_scan_temp = [[mz,intens] for mz,intens in ms1_scan if intens > 0.0]
		groups = []
		old_mass = 0.0
		

		while len(ms1_scan_temp) != 0:
			old_mass = 0.0
			t_group = []
			rem_indexes = []
			for index_list,spectrum in enumerate(ms1_scan_temp):
				mz,intens = spectrum
				if old_mass == 0.0:
					t_group.append([mz,intens])
					old_mass = mz
					threshold = (old_mass*self.ppm)/1000000
					rem_indexes.append(index_list)
					new_mass_lower = (old_mass+peak_dist)-threshold
					new_mass_upper = (old_mass+peak_dist)+threshold
					continue

				if mz > new_mass_lower and mz < new_mass_upper:
					t_group.append([mz,intens])
					
					old_mass = mz
					rem_indexes.append(index_list)
					new_mass_lower = (old_mass+peak_dist)-threshold
					new_mass_upper = (old_mass+peak_dist)+threshold
					continue

				elif mz > new_mass_upper:
					groups.append(t_group)
					for del_index in rem_indexes[::-1]:
						del ms1_scan_temp[del_index]
					print("Group size: %s" % (len(t_group)))
					print(len(ms1_scan_temp))
					break



		print(ms1_scan)
		raw_input("stp[")

	def align(self,isotopic_distribution,ms1_scan,normalize_exp=True,max_m = 3):
		counter = 0
		m_0_intens = 0.0

		aligned_m = {}
		exp_intens = []

		for m,isotope_peak in enumerate(isotopic_distribution):
			if m == 0:
				try: m_0_intens = self.peak_to_spec(ms1_scan,isotope_peak[0])[1]
				# Moving from ppm causes a problem... So this can be considered just outside of the ppm range from the experimental peaks' perspective
				except TypeError: return(1.0)
				continue

			if m > max_m:
				max_error = sum([max([isotopic_distribution[pos][1],aligned_m[pos]]) for pos in aligned_m.keys()])
				return(sum(aligned_m.values())/max_error)

			align_res = self.peak_to_spec(ms1_scan,isotope_peak[0])
			if not align_res:
				aligned_m[m] = isotope_peak[1]
				continue
			if normalize_exp:
				aligned_m[m] = abs((align_res[1]/m_0_intens)-isotope_peak[1])
			else:
				aligned_m[m] = abs(align_res[1]-isotope_peak[1])

		max_error = sum([max([isotopic_distribution[pos][1],aligned_m[pos]]) for pos in aligned_m.keys()])
		return(sum(aligned_m.values())/max_error)

	def peak_to_spec(self,spectrum,peak,threshold=False,mz_loc=0):
		if not threshold:
			threshold = (peak*self.ppm)/1000000
		curr_sel = int(len(spectrum)/2)

		if abs(spectrum[curr_sel][mz_loc] - peak) < threshold:
			return(spectrum[curr_sel])
		elif curr_sel == 0 or curr_sel == int(len(spectrum)):
			return(False)
		elif (spectrum[curr_sel][mz_loc] - peak) > 0.0:
			return(self.peak_to_spec(spectrum[0:curr_sel],peak))
		else:
			return(self.peak_to_spec(spectrum[curr_sel:],peak))


def peak_to_spec(spectrum,peak,threshold=False,mz_loc=1,orig_spec=[],pos_in_orig_spec=0,return_all_peaks=False,ppm=20):
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
			#print("no return list")
			print(str(ret_list))
		return(ret_list)

	if not threshold:
		threshold = (peak*ppm)/1000000		

	curr_sel = int(len(spectrum)/2)
	if return_all_peaks:
		if orig_spec == []:
			pos_in_orig_spec = curr_sel
			orig_spec = spectrum

	if abs(spectrum[curr_sel][mz_loc] - peak) < threshold:
		return(spectrum[curr_sel],abs(spectrum[curr_sel][mz_loc] - peak))
	elif curr_sel == 0 or curr_sel == int(len(spectrum)):
		return(False,False)
	elif (spectrum[curr_sel][mz_loc] - peak) > 0.0:
		if return_all_peaks:
			if curr_sel % 2 != 0: pos_in_orig_spec = pos_in_orig_spec-int(len(spectrum[0:curr_sel])/2)-1
			else: pos_in_orig_spec = pos_in_orig_spec-int(len(spectrum[0:curr_sel])/2)
		return(peak_to_spec(spectrum[0:curr_sel],peak,mz_loc=mz_loc,pos_in_orig_spec=pos_in_orig_spec,orig_spec=orig_spec))
	else:
		if return_all_peaks:
			if curr_sel % 2 != 0: pos_in_orig_spec = pos_in_orig_spec+int(len(spectrum[curr_sel:])/2)
			else: pos_in_orig_spec = pos_in_orig_spec+int(len(spectrum[curr_sel:])/2)
		return(peak_to_spec(spectrum[curr_sel:],peak,mz_loc=mz_loc,pos_in_orig_spec=pos_in_orig_spec,orig_spec=orig_spec))

def para_align_candidate(candidate,spectrum):
	temp_alignment = []
	for db_peak in candidate[2].ms2:
		peak = db_peak[0]
		res_alignment,peak_error = peak_to_spec(spectrum,peak,mz_loc=0)
		if res_alignment:
			temp_alignment.append([db_peak,res_alignment,peak_error])
	if len(temp_alignment) > 0:
		result_alignment = [candidate,temp_alignment,candidate[2].ms2]
	else:
		result_alignment = []
	return(result_alignment)
		
class MS2Filter():
	def __init__(self,fa_dict,class_spec_file="class_labels.csv",ppm=5):
		self.fa_dict = fa_dict
		self.sorted_fa = sorted(self.fa_dict.items(), key=itemgetter(1))

		self.min_fa = self.sorted_fa[0][1]
		self.max_fa = self.sorted_fa[-1][1]

		self.nl_dict,self.frag_dict = self.read_class_spec_file(class_spec_file)

		self.sorted_frag = sorted(self.frag_dict.items(), key=itemgetter(1))

		self.ppm = ppm

	def read_class_spec_file(self,class_spec_file):
		nl_dict = {}
		frag_dict = {}

		df = pd.read_csv(class_spec_file)
		for line in df.itertuples():
			if line[2] == "NL": nl_dict[line[7]+"|"+line[6]] = line[3]
			elif line[2] == "FRAG": frag_dict[line[7]+"|"+line[6]] = line[3]
			else: logging.warning("Unknown type for fragment: %s" % (line[2]))
		return(nl_dict,frag_dict)

	def find_ms2_fa(self,spectrum_mz):
		for peak in spectrum_mz:
			mass_error_threshold = (peak*self.ppm)/1000000
			if peak-mass_error_threshold < self.min_fa: continue
			if peak+mass_error_threshold > self.max_fa: continue
			res = self.peak_to_spec(self.sorted_fa,peak,threshold=mass_error_threshold)
			if res != False: print(res)

	def find_ms2_class(self,spectrum_mz,prec_mass):
		#print(prec_mass)
		sorted_nl = sorted([(name,prec_mass-nl) for name,nl in self.nl_dict.items()],key=itemgetter(1))

		#print(sorted_nl)

		for peak in spectrum_mz:
			mass_error_threshold = (peak*self.ppm)/1000000
			res = self.peak_to_spec(sorted_nl,peak,threshold=mass_error_threshold)
			#if res != False: 
			#	print(res)
				#raw_input()

		#print(self.sorted_frag)
		for peak in spectrum_mz:
			mass_error_threshold = (peak*self.ppm)/1000000
			res = self.peak_to_spec(self.sorted_frag,peak,threshold=mass_error_threshold)
			#if res != False: 
			#	print(res)
				#raw_input()
			
	def peak_to_spec(self,spectrum,peak,threshold=False,mz_loc=1,orig_spec=[],pos_in_orig_spec=0,return_all_peaks=True):
		def _look_around_matched_peak(spectrum,peak,pos,rev=False):
			error = 0.0
			ret_list = [[0.0,0.0]]
			while error < threshold and pos < len(spectrum)-1 and pos > 1:
				if rev: pos = pos-1
				else: pos = pos+1
				error = abs(spectrum[pos][mz_loc]-peak)
				if error < threshold:
					ret_list.append(spectrum[pos])
			if len(ret_list) > 1:
				index,val = max(enumerate(map(itemgetter(-1), ret_list)),key=itemgetter(1))
				ret_list = [ret_list[index]]
			return(ret_list)

		if not threshold:
			threshold = (peak*self.ppm)/1000000		

		curr_sel = int(len(spectrum)/2)
		if return_all_peaks:
			if orig_spec == []:
				pos_in_orig_spec = curr_sel
				orig_spec = spectrum

		if abs(spectrum[curr_sel][mz_loc] - peak) < threshold:
			if not return_all_peaks:
				candidate_peaks = [list(spectrum[curr_sel])]
				candidate_peaks.extend(_look_around_matched_peak(orig_spec,peak,pos_in_orig_spec))
				candidate_peaks.extend(_look_around_matched_peak(orig_spec,peak,pos_in_orig_spec,rev=True))
				index,val = max(enumerate(map(itemgetter(-1), candidate_peaks)),key=itemgetter(1))
				return(candidate_peaks[index],abs(spectrum[curr_sel][mz_loc] - peak))
			else:
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

	def align_db(self,spectrum,db_list,parallel=False,num_cores=6):	
		result_alignment = []
		peak_memory = {}
		
		if parallel == True:
			try:
				self.para_obj
			except:
				self.para_obj = Parallel(n_jobs=num_cores) #,backend="threading"

			res = self.para_obj(delayed(para_align_candidate)(candidate,spectrum) for candidate in db_list)
			result_alignment = [r for r in res if len(r) > 0]
			#print(result_alignment)
		else:
			for candidate in db_list:
				temp_alignment = []
				for db_peak in candidate[2].ms2:
					peak = db_peak[0]
					# Simply trying to access the key and otherwise align still seems to be the fastest option
					try:
						res_alignment,peak_error = peak_memory[peak]
					except KeyError: #KeyError
						res_alignment,peak_error = self.peak_to_spec(spectrum,peak,mz_loc=0)
						peak_memory[peak] = [res_alignment,peak_error]
					if res_alignment:
						temp_alignment.append([db_peak,res_alignment,peak_error])

				if len(temp_alignment) > 0:
					result_alignment.append([candidate,temp_alignment,candidate[2].ms2])

		return(result_alignment)
