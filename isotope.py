from pyteomics.mass.mass import isotopic_composition_abundance
from pyteomics.mass.mass import isotopologues
from pyteomics import mass
from numpy.polynomial.polynomial import Polynomial
from numpy.polynomial.polynomial import polypow
from itertools import groupby
import ast

 # TODO same molecular formula, but different mw

class isotope_filters():
	def __init__(self,dalt_diff_lookup_bin=1,prec_calc_isotopic_distribution="isotope_lookup.txt"):
		self.dalt_diff_lookup_bin = dalt_diff_lookup_bin

		self.isotope_dict = {}
		self.isotope_fast_lookup = {}

		if len(prec_calc_isotopic_distribution) > 0:
			self.read_prec_isotopic_distribution(open(prec_calc_isotopic_distribution))

	def read_prec_isotopic_distribution(self,infile):
		for line in infile:
			split_line = line.strip().split(";")
			chem_form = split_line[0]
			mass_intens = split_line[1]
			try: mass_intens = ast.literal_eval(mass_intens)
			except ValueError: continue

			self.isotope_dict[chem_form] = mass_intens

			loc_dict = int(mass_intens[0][0]) - int(mass_intens[0][0]) % self.dalt_diff_lookup_bin

			if loc_dict in self.isotope_fast_lookup.keys():
				self.isotope_fast_lookup[loc_dict].append(chem_form)
			else:
				self.isotope_fast_lookup[loc_dict] = [chem_form]

	def add_isotope(self,chem_form):
		mass_intens = self._get_isotopes_new(chem_form)
		self.isotope_dict[chem_form] = mass_intens

		loc_dict = int(mass_intens[0][0]) - int(mass_intens[0][0]) % self.dalt_diff_lookup_bin

		if loc_dict in self.isotope_fast_lookup.keys():
			self.isotope_fast_lookup[loc_dict].append(chem_form)
		else:
			self.isotope_fast_lookup[loc_dict] = [chem_form]

	def get_isotope_distribution(self,chem_form):
		try: return(self.isotope_dict[chem_form])
		except KeyError: 
			self.add_isotope(chem_form)
			return(self.isotope_dict[chem_form])

	def _get_isotopes_new(self,chem_form,lim_atoms=10,min_prob=0.00000000001):
		mono_iso = mass.calculate_mass(chem_form)
		masses = [mono_iso+(i*1.0033548378) for i in range(lim_atoms)]


		chem_form_ion = ""
		for i,c in enumerate(chem_form):
			if i+1 >= len(chem_form):
				if c.isdigit(): chem_form_ion += c
				else: 
					chem_form_ion += c
					chem_form_ion += "1"
			elif c.isdigit(): chem_form_ion += c
			elif c.isupper() and chem_form[i+1].isdigit(): chem_form_ion += c
			elif c.isupper() and chem_form[i+1].isupper(): 
				chem_form_ion += c
				chem_form_ion += "1"
			#elif chem_form[i+1].isdigit():
			#	#print(chem_form,chem_form_ion)
			#	chem_form_ion += c
			elif chem_form[i].isupper() and not chem_form[i+1].isupper():
				print(chem_form,chem_form_ion)
				chem_form_ion += chem_form[i]
				chem_form_ion += chem_form[i+1]
				chem_form_ion += "1"

		list_chem= [''.join(g) for _, g in groupby(chem_form_ion, str.isalpha)]
		chem_form = dict(zip(list_chem[::2],map(int,list_chem[1::2])))

		print(chem_form)

		#print(chem_form["C"])
		#try: 
		if chem_form["C"] < 100:
			probs = Polynomial([prob for m,prob in mass.nist_mass["C"].values() if prob > min_prob][1:])**chem_form["C"]
		else:
			tot_c = chem_form["C"]
			probs = Polynomial([prob for m,prob in mass.nist_mass["C"].values() if prob > min_prob][1:])**99
			tot_c -= 99
			while tot_c > 99:
				probs *= probs**99
				tot_c -= 99
			probs *= probs**tot_c
		#print(list(probs))
		#print(mass.nist_mass["C"].values())
		#except ValueError: probs = Polynomial([prob for m,prob in mass.nist_mass["C"].values() if prob > 0.0001][1:])**99
		#probs = probs.truncate(num_atom)
		#input("first_stop")
		for atom,num_atom in chem_form.items():
			if atom == "C": continue
			if chem_form[atom] < 100:
				probs *= Polynomial([prob for m,prob in mass.nist_mass[atom].values() if prob > min_prob][1:])**chem_form[atom]
			else:
				tot_c = chem_form[atom]
				probs *= Polynomial([prob for m,prob in mass.nist_mass[atom].values() if prob > min_prob][1:])**99
				tot_c -= 99
				while tot_c > 99:
					probs *= probs**99
					tot_c -= 99
				probs *= probs**tot_c

		max_prob = max(list(probs)[:lim_atoms])
		probs = [i/max_prob for i in list(probs)[:lim_atoms]]
		return(list(zip(masses,probs)))

	def _get_isotopes(self,chem_form):
		masses = []
		probs = []
		for i,abundance in isotopologues(chem_form,overall_threshold=0.00001,isotope_threshold=0.0001,report_abundance=True): #,isotope_threshold=0.0005
			ab = isotopic_composition_abundance(i)
			masses.append(i.mass())
			probs.append(isotopic_composition_abundance(i))
			print(abundance)
			#print(masses,abundance)
			#print(probs)
		max_prob = max(probs)
		probs = [rel_intensity/max_prob for rel_intensity in probs]
		#print(probs)
		return(list(zip(masses,probs)))

class Isotopic_dists():
	def __init__(self,diff_peaks=1.0033548378,min_intens=0,mass_threshold=0.01):
		self.diff_peaks = diff_peaks
		self.min_intens = min_intens
		self.mass_threshold = mass_threshold

	def find_isotopic_distributions(self,mz_array,intensity_array):
		mz_array_filt = [mz for mz,intens in zip(mz_array,intensity_array) if intens > self.min_intens]
		for index,peak in enumerate(mz_array_filt):
			temp_diff_peak = self.diff_peaks
			temp_peak = peak
			counter = 0
			print("====++++=====")
			for peak_iso in mz_array_filt[index+1:]:
				if abs((temp_peak+temp_diff_peak)-peak_iso) < self.mass_threshold:
					print(index,peak,temp_peak,peak_iso,(temp_peak+temp_diff_peak)-peak_iso,intensity_array[index+counter+1],max(intensity_array))
					#temp_diff_peak += self.diff_peaks
					temp_peak = peak_iso
					counter += 1
				if ((temp_peak+temp_diff_peak)-peak_iso) < (self.mass_threshold*-1):
					break
			#if counter > 5:
			#	raw_input("=====")

import pickle
if __name__ == "__main__":
	itf = isotope_filters()
	infile = open("formula2.csv")
	#itf.add_isotope("H156C87O17P2")


	for line in infile:
		line = line.strip()
		print(itf.add_isotope(line))

	outfile = open("isotope_lookup.txt","w")
	for chem_formula,distribution in itf.isotope_dict.items():
		outfile.write("%s;%s\n" % (chem_formula,distribution))
	outfile.close()

	