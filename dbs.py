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

def formula_dict(chem):
	s = re.findall('([A-Z][a-z]?)([0-9]*)', chem)
	try: return(dict([[element, int(count)] for element,count in s]))
	except:
		temp_dict = {}
		for element,count in s:
			if count == "":
				count = 1
			temp_dict[element] = int(count)
		return(temp_dict)


def subtract_dicts(main_dict,diff_dict,subtract=True):
	for kdd in diff_dict.keys():
		if subtract: main_dict[kdd] -= diff_dict[kdd]
		else: 
			try: main_dict[kdd] += diff_dict[kdd]
			except KeyError: main_dict[kdd] = diff_dict[kdd]
	return(main_dict)
  
def ion_chem_form(chem_form,ion):
	chem_form_dict = formula_dict(chem_form)
	ion = ion.replace("[","").replace("]","").replace("M","")
	subtract_form = {}
	add_form = {}
	re_res = re.findall("[\+|\-][a-zA-Z0-9]*",ion)

	for form in re_res:
		if form[0] == "-": chem_form_dict = subtract_dicts(chem_form_dict,formula_dict(form[1:]))
		if form[0] == "+": chem_form_dict = subtract_dicts(chem_form_dict,formula_dict(form[1:]),subtract=False)

	return("".join([element+str(count) for element,count in chem_form_dict.items()]))
	
def determine_chunks_prec(f_name,chunks=10):
	with open(f_name) as infile:
		pre_c_masses = []
		pre_c_mass = 0.0

		for line in infile:
			line = line.strip()
			if len(line) == 0:
				if pre_c_mass == 0.0: continue

				pre_c_masses.append(pre_c_mass)
				pre_c_mass = 0.0

			elif ":" in line:
				if line.startswith("PRECURSORMZ"):
					pre_c_mass = float(line.split(": ")[1])
			else:
				if line=="\x1a": #EOF
					continue
	pre_c_masses.sort()
	ranges_mz = [pre_c_masses[i] for i in range(0,len(pre_c_masses),int(len(pre_c_masses)/chunks))]
	ranges_mz.append(pre_c_masses[-1])
	ranges_mz_twod = [[ranges_mz[i],ranges_mz[i+1]]for i in range(len(ranges_mz)-1)]
	return(ranges_mz_twod)


class DBSelf_entry():
	def __init__(self,
				 name="",
				 ion="",
				 mw=0.0,
				 chem_form="",
				 num_ms2_peaks=0,
				 f_acyl_lengths=[],
				 unsats=[],
				 ms2=[]):

		self.name = name
		self.ion = ion
		self.mw = mw
		self.chem_form = chem_form
		self.num_ms2_peaks = num_ms2_peaks
		self.ms2 = ms2
		self.f_acyl_lengths = f_acyl_lengths
		self.unsats = unsats


	def __str__(self):
		ret_string = []
		ret_string.append("================")
		ret_string.append("")
		ret_string.append("Lipid: %s" % (self.name))
		ret_string.append("MW: %s" % (self.mw))
		ret_string.append("Formula: %s" % (self.chem_form))
		ret_string.append ("")
		for f in self.ms2:
			ret_string.append("%s\t%s\t%s" % (f[0],f[1],f[2]))
		ret_string.append("")
		ret_string.append("================")

		return("\n".join(ret_string))

class DBSelf():
	def __init__(self,
				 f_names=["db.msp"],
				 min_acyl_length=0,
				 exclude_lyso=False,
				 include_ions=["[M-H]-"], #,"[M+]","[M+H]+","[M+NH4]+","[M-H]-","[M-2H](2-)","[M-Ac-H]-","[M+Na2-H]+","[M+]","[M+NH4]+","[M+Na]+","[M-2H](2-)","[M-Ac-H]-"                   "[M+]","[M+H]+","[M+NH4]+","[M-H]-","[M-2H](2-)","[M-Ac-H]-","[M+Na2-H]+","[M+]","[M+NH4]+","[M+Na]+","[M-2H](2-)","[M-Ac-H]-"
				 include_class=["PE","GPSer","GPCho","PC","GPA","PE","GPIns","GPEtn","GPGro"],  #,"SM","TG","CL",         #,"SM","TG","CL","GPSer","GPCho","PC","GPA","PE","GPIns","GPEtn","GPGro
				 aggregate_acyls=False,
				 use_simplified_names=True,
				 dalt_diff_lookup_bin=1,
				 pos_ions=True,
				 neg_ions=True,
				 lower_limit_mz=0.0,
				 upper_limit_mz=10000.0):

		self.f_names = f_names
		self.min_acyl_length = min_acyl_length
		self.exclude_lyso = exclude_lyso
		self.include_ions = include_ions
		self.include_class = include_class
		self.use_simplified_names = use_simplified_names
		self.dalt_diff_lookup_bin = dalt_diff_lookup_bin
		self.aggregate_acyls = aggregate_acyls
		self.pos_ions = pos_ions
		self.neg_ions = neg_ions
		self.lower_limit_mz = lower_limit_mz
		self.upper_limit_mz = upper_limit_mz
		
		self.dbs_dict = {}
		self.ms1_dict = {}
		self.ms1_dict_lookup = {}

		self.tot_entr_read = 0

		if len(self.f_names) > 0:
			for f_name in f_names:
				try: self.read_dbs(f_name,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz)
				except: continue
				#self.read_dbs(f_name)

	def __str__(self):
		ret_string = []

		ret_string.append("Filenames: %s" % (self.f_names))
		ret_string.append("Min. mz: %s" % (self.f_names))
		ret_string.append("Min acyl length: %s" % (self.min_acyl_length))
		ret_string.append("Exclude lyso: %s" % (self.exclude_lyso))
		ret_string.append("Include ions: %s" % (self.include_ions))
		ret_string.append("Include lipid classes: %s" % (self.include_class))
		ret_string.append("Use simplified names: %s" % (self.use_simplified_names))
		ret_string.append("Lookup diff: %s Da" % (self.dalt_diff_lookup_bin))
		ret_string.append("Total entries read: %s" % (self.tot_entr_read))
		
		return("\n".join(ret_string))

	def read_dbs(self,f_name,lower_limit_mz=0.0,upper_limit_mz=10000.0,max_counter=-1):
		with open(f_name) as infile:
			fragments = []
			pre_c_mass = 0.0
			name = ""
			ion_type = ""
			chem_form_ion = ""
			counter = 0
			ignore_entry = False

			for line in infile:
				line = line.rstrip()
				if ignore_entry and line != "":
					continue
				elif len(line) == 0:
					if len(name) == 0: continue
					if (ion_type[-1] == "-" and not self.neg_ions) or (ion_type[-1] == "+" and not self.pos_ions): 
						fragments = []
						pre_c_mass = 0.0
						name = ""
						ion_type = ""
						ignore_entry = False
						continue
					if not (pre_c_mass >= lower_limit_mz and pre_c_mass <= upper_limit_mz): 
						fragments = []
						pre_c_mass = 0.0
						name = ""
						ion_type = ""
						ignore_entry = False
						continue
					
					class_name = name.split("(")[0]
					#print(name) #.split("(")[1].split(")")[0]

					# TODO Make it possible to retrieve more than two FA
					#f_acyl_lengths = map(int,itemgetter(*[0,2])(re.split("\:|\_|\/",name.split("(")[1].split(")")[0])))
					#unsats = map(int,itemgetter(*[1,3])(re.split("\:|\_|\/",name.split("(")[1].split(")")[0])))

					#f_acyl_lengths_error = [a for a in f_acyl_lengths if a < self.min_acyl_length and a != 0]
					
					#if (len(class_name) == 0) or \
					#	(ion_type not in self.include_ions) or \
					#	(len([c for c in self.include_class if c in name]) == 0) or \
					#	(self.exclude_lyso and "/0:0" in name) or \
					#	(len(f_acyl_lengths_error) > 0):
#
#						fragments = []
#						pre_c_mass = 0.0
#						name = ""
#						ion_type = ""
#						continue
			
					#simplified_name = _simplify_name(class_name,f_acyl_lengths,unsats)

					new_entry = DBSelf_entry(name=name,
												 ion=ion_type,
												 mw=pre_c_mass,
												 chem_form=chem_form_ion,
												 num_ms2_peaks=num_peaks,
												 ms2=fragments,
												 f_acyl_lengths=20, #f_acyl_lengths
												 unsats=6) #unsats
					counter += 1
					self.dbs_dict["%s|%s" % (name,ion_type)] = new_entry

					loc_dict = int(pre_c_mass) - int(pre_c_mass) % self.dalt_diff_lookup_bin

					if loc_dict in self.ms1_dict_lookup.keys():
						self.ms1_dict_lookup[loc_dict]["%s|%s" % (name,ion_type)] = new_entry
					else:
						self.ms1_dict_lookup[loc_dict] = {}
						self.ms1_dict_lookup[loc_dict]["%s|%s" % (name,ion_type)] = new_entry

					self.tot_entr_read += 1
					if counter > max_counter and max_counter != -1:
						break
					fragments = []
					pre_c_mass = 0.0
					ignore_entry = False
					name = ""
					ion_type = ""

				elif ":" in line:
					if line.startswith("PRECURSORMZ"):
						pre_c_mass = float(line.split(": ")[1])
						if pre_c_mass <= lower_limit_mz or pre_c_mass >= upper_limit_mz:
							ignore_entry = True
					if line.startswith("Ion: "):
						ion_type = line.split(": ")[1]
					if line.startswith("Name: "):
						name = line.split(": ")[1]
					if line.startswith("Formula: "):
						chem_form_native = line.split(": ")[1]
						chem_form_ion = ion_chem_form(chem_form_native,ion_type)
					if line.startswith("Num Peaks:"):
						num_peaks = int(line.split(": ")[1])
				else:
					if line=="\x1a": #EOF
						continue
					splitline = line.split(" ")
					fragments.append([float(splitline[0]),float(splitline[1]),splitline[2].replace("\"","")])

# TODO this can be removed?
class PrecursorFilter():
	def __init__(self,db,ppm=10):
		self.db = db
		self.ppm = ppm

	def retrieve_entry_pre_c_mass(self,pre_c_mass):
		mass_error_threshold = (pre_c_mass*self.ppm)/1000000

		ret_entries = []

		loc_dict = int(pre_c_mass) - int(pre_c_mass) % self.db.dalt_diff_lookup_bin
		loc_dict_lower = (int(pre_c_mass-mass_error_threshold)) - (int(pre_c_mass-mass_error_threshold)) % self.db.dalt_diff_lookup_bin
		loc_dict_upper = (int(pre_c_mass+mass_error_threshold)) - (int(pre_c_mass+mass_error_threshold)) % self.db.dalt_diff_lookup_bin

		# TODO set does not have to be list
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

if __name__ == "__main__":
	logging.basicConfig(filename="prec_filter.log",
						level=logging.DEBUG,
						filemode="w",
						format="%(levelname)s:%(created)f:%(asctime)s:%(message)s")

	logging.info("Reading the LPB database ...")
	dbs = DBSelf()
	logging.info("Done reading the LPB database ...")
	logging.info(dbs)
	print(dbs)
	input()
	
	#step_three_df = pd.read_csv("stepone_new.csv")
	#precf = Precursor_filter(lpb)
	
	#prec_filt_result = []
	#for index,row in step_three_df.iterrows():
	#	if (index % 10000==0):
	#		logging.info("Analyzing row number and m/z: %s - %s" % (index,row["mz"]))
	#	prec_hits = precf.retrieve_entry_pre_c_mass(row["mz"])
	#	for hit in prec_hits:
	#		prec_filt_result.append([row["mz"],hit[2].mw,hit[1],hit[0].split("|")[0],hit[2].chem_form,hit[0].split("|")[1]])
	#
	#prec_filt_result = pd.DataFrame(prec_filt_result)
	#prec_filt_result.columns = ["Input Mass","Matched Mass","Delta","Abbreviation","Formula","Ion"]
	#prec_filt_result.to_excel("batch_results.xlsx",index=False)
