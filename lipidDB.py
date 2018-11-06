from copy import deepcopy
from collections import Counter
import itertools
import exrex
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors
from pyteomics.mass.mass import calculate_mass
import operator
import re

mw_dict = {}
smiles_to_mw = {}
smiles_to_mol_form = {}

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
			except KeyError: diff_dict[kdd]
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

def update_progress(progress):
    print('\r[{0}] {1}%'.format('#'*(progress/10), progress))

class ConfigDB():
	def __init__(self,config_file):
		print(config_file)
		self.struct_config = {}
		self.structures = {}
		self.modifications = {}
		self.group = {}
		self.fragmentation_rules = []
		self.ions_ms1 = []
		self.name = []

		self.parse_config_file(config_file)

	def parse_config_file(self,config_file):
		with open(config_file) as infile:
			section = ""
			for line in infile:
				line = line.strip()
				if len(line) == 0: continue
				if line.startswith("#"): continue
				if line.startswith("!"):
					if line == "!structural_configuration":	section = "structural_configuration"
					if line == "!structure_def": section = "structure_def"
					if line == "!modifications": section = "modifications"
					if line == "!group_def": section = "group_def"
					if line == "!fragmentation_ms2": section = "fragmentation_ms2"
					if line == "!ions_ms1": section = "ions_ms1"
					if line == "!name": section = "name"
					continue
				
				if section == "structural_configuration":
					print(line)
					name_struct = line.split(" = ")[0].lstrip().rstrip()
					struct_descr = line.split(" = ")[1].lstrip().rstrip()
					self.struct_config[name_struct] = struct_descr
				elif section == "structure_def":
					name_def = line.split(" = ")[0].lstrip().rstrip()
					struct_descr = line.split(" = ")[1].lstrip().rstrip()
					if struct_descr.startswith("exrex:"):
						struct_descr = exrex.generate(struct_descr.split("exrex:")[1])
						
						#Next piece of code is to remove any duplicate molecular formulas
						unique_struct_set = set()
						unique_struct = []
						
						for x in struct_descr:
							if str(Counter(x)) not in unique_struct_set:
								unique_struct_set.add(str(Counter(x)))
								unique_struct.append(x)
						struct_descr = unique_struct

					self.structures[name_def] = struct_descr

				elif section == "modifications":
					name_struct = line.split(" = ")[0].lstrip().rstrip()
					struct_descr = line.split(" = ")[1].lstrip().rstrip()
					self.modifications[name_struct] = struct_descr

				elif section == "group_def":
					name_struct = line.split(" = ")[0].lstrip().rstrip()
					struct_descr = line.split(" = ")[1].lstrip().rstrip()
					self.group[name_struct] = struct_descr

				elif section == "fragmentation_ms2":
					limit_frag = []
					if " = " in line:
						try: limit_frag = line.split(" = ")[2].split("|")
						except IndexError: limit_frag = []
					frag_rule = line.split(" = ")[1]
					ion_rule = line.split(" = ")[0]
					self.fragmentation_rules.append([ion_rule,frag_rule,limit_frag])

				elif section == "ions_ms1":
					limit_ion = []
					if " = " in line:
						limit_ion = line.split(" = ")[1].split("|")
					ion_rule = line.split(" = ")[0]
					self.ions_ms1.append([ion_rule,limit_ion])
					
				elif section == "name":
					limit_name = []
					if " = " in line:
						limit_name = line.split(" = ")[1].split("|")
					self.name.append([line.split(" = ")[0],limit_name])


	def create_db(self,outfile_name="db_noether.msp"):
		lipid_db = []
		counter = 0
		fragmentar = ThaFragmenter(self.fragmentation_rules)
		outfile_msp = open(outfile_name,"w")
		outfile_smiles = open(outfile_name+"_tosmiles.csv","w")

		for name,struct in self.struct_config.items():
			temp_structs = [list([])]
			temp_structs_names = [list([])]
			
			
			

			for char in struct:
				if char == "[":
					str_element = ""
					start_read = True
				elif char == "]":
					
					if type(self.structures[str_element]) == list:
						#Take the product of two lists and flatten down one level per item in that list to get the list of smiles 
						temp_structs = [list(itertools.chain(*i)) for i in list(itertools.product(*[temp_structs,[[s] for s in self.structures[str_element]]]))]
						temp_structs_names = [list(itertools.chain(*i)) for i in list(itertools.product(*[temp_structs_names,[[str_element]]*len(self.structures[str_element])]))]
					else:
						[j.append(self.structures[str_element]) for j in temp_structs]
						[j.append(str_element) for j in temp_structs_names]
					str_element = ""
					start_read = False
				elif start_read:
					str_element += char
				elif char != " ":
					[j.append(char) for j in temp_structs]
					[j.append("-") for j in temp_structs_names]
			
			
			
			for n,s in zip(temp_structs_names,temp_structs):
				#temp_structs_names_group = []
				group_names = []
				
				for j in n:
					try:
						group_names.append(self.group[j])
					except KeyError:
						group_names.append("-")
				temp_structs_names_group = group_names
				
				possible_ions = []
				for ion,rules in self.ions_ms1:
					#rules = [r.split(":")[1] for r in rules]
					#print(" ***** ")
					#print(rules)
					
					viol_rules = [r for r in rules if r not in temp_structs_names_group and r not in n]
					#print(" viols: ")
					#print(viol_rules)
					# TODO something goes wrong here...
					#print(temp_structs_names_group)
					#print(n)
					#print(viol_rules)
					#print(rules)
					#raw_input()
					if len(viol_rules) != len(rules):
						possible_ions.append(ion)
				#print("possible_ions:")
				#print(possible_ions)

				temp_lipid = Lipid("".join(s),s,n,group_names,calc_name=self.name,possible_ions=possible_ions,modifications=self.modifications)
				counter += 1

				#print(self.fragmentation_rules)
				
				msp_entries = fragmentar.apply_fragmentation(temp_lipid)

				for msp in msp_entries:
					outfile_msp.write(msp)

				outfile_smiles.write("%s,%s\n" % (temp_lipid.name[0],temp_lipid.smiles))
				
				if counter % 1000 == 0:
					print(" Analysed total number of lipids: %s" % (counter))
				#print(temp_lipid.ions_mw_neg)
					#print(temp_lipid.mol_wt)

				# TODO check if this var should be reset

				#print(t1-t2)
				#print(t2-t3)
				#print(t3-t4)

				#print("----")

		outfile_msp.close()
		outfile_smiles.close()


class ThaFragmenter():
	def __init__(self,fragmentation_rules,intensity=100):
		self.fragmentation_rules = fragmentation_rules
		self.intensity = intensity



	def apply_fragmentation(self,lipid,return_msp_string=True):
		all_ions = list(lipid.ions_mw_neg.keys())
		all_ions.extend(list(lipid.ions_mw_pos.keys()))
		msp_entries = []

		for ion in all_ions:
			peaks = []
			peak_names = []
			intensity = []
			ion_mw = 0.0

			stripped_ion = ion.split("]")[0]+"]"
			for rule_ion,t_rule,limit_name in self.fragmentation_rules:
				# Todo make the double mapping possible (head = "PS")
				# for now a quick solution that only checks for "PS"
				#
				#if "PCP" in lipid.sub_names:
				#	print(rule_ion,t_rule,limit_name)
				#	print(lipid.sub_names_groups)
				#	print(lipid.sub_names)
				#	print('====')

				if rule_ion.split("]")[0]+"]" != stripped_ion:
					continue
				
				#print(rule_full,limit_name)
				#raw_input()


				#print(rule)
				#raw_input(" ;;;;;; ")
				if len(limit_name) != 0:
					viol_bool = []
					for lm in limit_name:
						if lipid.group_in_lipid(lm) == False:
							viol_bool.append(False)
						else:
							viol_bool.append(True)
					#if "PCP" in lipid.sub_names:
					#	print(viol_bool)
					if not any(viol_bool):
						continue

				start_ion = t_rule.split("]")[0]+"]"
				
				# TODO what about 2-? Not very generic code!
				if start_ion+"-" in lipid.ions_mw_neg.keys():
					start_mw = lipid.ions_mw_neg[start_ion+"-"]
					ion_mw = start_mw
					t_rule = t_rule[len(start_ion):]
				elif start_ion+"+" in lipid.ions_mw_pos.keys(): 
					start_mw = lipid.ions_mw_pos[start_ion+"+"]
					ion_mw = start_mw
					t_rule = t_rule[len(start_ion):]
				else: 
					start_mw = 0.0
				
				change_mw = 0.0
				group_name = ""
				start_reading_name = False
				subtract = True
				chem_form = ""

				#print(rule)
				#print(t_rule)
				#raw_input(" pfdlsp")

				for char_rule in t_rule:
					if char_rule == "M":
						start_mw = lipid.mol_wt
					elif char_rule == "-":
						if chem_form == "H" or chem_form == "H ": chem_form = "H1"
						try: 
							change_mw = mw_dict[chem_form]
						except KeyError:
							mw_dict[chem_form] = calculate_mass(chem_form)
							change_mw = mw_dict[chem_form]
						if subtract: start_mw -= change_mw
						else: start_mw += change_mw
						chem_form = ""
						change_mw = 0.0
						subtract = True
					elif char_rule == "+":
						if chem_form == "H" or chem_form == "H ": chem_form = "H1"
						try: 
							change_mw = mw_dict[chem_form]
						except KeyError:
							mw_dict[chem_form] = calculate_mass(chem_form)
							change_mw = mw_dict[chem_form]
						if subtract: start_mw -= change_mw
						else: start_mw += change_mw
						chem_form = ""
						change_mw = 0.0
						subtract = False
					elif char_rule == "}":
						if start_mw == 0.0:
							try:
								start_mw = smiles_to_mw[lipid.get_smiles_sub(group_name)]
								#get_smiles_sub(group_name)
								#start_mw = mw_dict[lipid.get_molecular_formula_sub(group_name)]
							except KeyError:
								smiles_to_mw[lipid.get_smiles_sub(group_name)] = calculate_mass(lipid.get_molecular_formula_sub(group_name))
								start_mw = smiles_to_mw[lipid.get_smiles_sub(group_name)]
							#print(lipid.name)
							#start_mw = calculate_mass(lipid.get_molecular_formula_sub(group_name))

							ion_mw = start_mw
							group_name = ""
							start_reading_name = False
							continue

						try:
							chem_form = smiles_to_mol_form[lipid.get_smiles_sub(group_name)]
						except KeyError:
							smiles_to_mol_form[lipid.get_smiles_sub(group_name)] = lipid.get_molecular_formula_sub(group_name)
							chem_form += smiles_to_mol_form[lipid.get_smiles_sub(group_name)]

						#chem_form = lipid.get_molecular_formula_sub(group_name)
						start_reading_name = False
					elif start_reading_name:
						group_name += char_rule
					elif char_rule == "{":
						start_reading_name = True
					elif char_rule == "[" or char_rule == "]":
						continue
					else:
						chem_form += char_rule

				if chem_form == "H" or chem_form == "H ": chem_form = "H1"

				try: 
					change_mw = mw_dict[chem_form]
				except KeyError:
					# TODO make a chemical formula sanitizer
					try: mw_dict[chem_form] = calculate_mass(chem_form)
					except: mw_dict[chem_form] = calculate_mass(chem_form.replace(" ",""))
					change_mw = mw_dict[chem_form]
				if subtract: start_mw -= change_mw
				else: start_mw += change_mw

				if start_mw < 0.0: continue
				peaks.append(start_mw)
				peak_names.append(t_rule)
				intensity.append(self.intensity)
			
			if len(peaks) == 0: continue
			msp_entry = self.get_msp_entry(lipid.name[0],peaks,peak_names,intensity,lipid.mol_wt,ion_mw,ion,lipid.chem_form)
			msp_entries.append(deepcopy(msp_entry))
		
		return(msp_entries)

	def get_msp_entry(self,name,peaks,peak_names,intensity,mol_wt,mol_wt_ion,name_ion,chem_form):
		ret_entry = ""
		ret_entry += "Name: %s\n" % (name)
		ret_entry += "Ion: %s\n" % (name_ion)
		ret_entry += "MW: %s\n" % (mol_wt)
		ret_entry += "PRECURSORMZ: %s\n" % (mol_wt_ion)
		ret_entry += "Num Peaks: %s\n" % (len(peaks))
		ret_entry += "Formula: %s\n" % (chem_form)
		ret_entry += "Formula Ion: %s\n" % (ion_chem_form(chem_form,name_ion))
		spectrum = zip(peaks,intensity,peak_names)
		spectrum = sorted(spectrum,key=operator.itemgetter(0))
		ret_entry += "\n".join([str(p)+" "+str(i)+" \""+pn+"\"" for p,i,pn in spectrum])
		ret_entry += "\n\n"
		return(ret_entry)

		
class Lipid():
	def __init__(self,smiles,sub_smiles,sub_names,sub_names_groups,calc_name=[],possible_ions=["[M-H]-","[M+H]+"],modifications={}):
		self.smiles = smiles
		self.sub_smiles = sub_smiles
		self.sub_names = sub_names
		self.sub_names_groups = sub_names_groups
		self.rdkit_mol = Chem.MolFromSmiles("".join(smiles))
		self.possible_ions = possible_ions
		self.calc_name = calc_name
		# todo quick and dirty solution for this... think of a better way.
		self.modifications = modifications

		self.chem_form = ""
		self.mol_wt = 0.0
		self.ions_mw_neg = {}
		self.ions_mw_pos = {}

		self.calculate_properties(self.possible_ions)
		self.name = self.get_lipid_name()
		#get name lipid... make function

	def calculate_properties(self,possible_ions):
		try:
			self.chem_form = Chem.rdMolDescriptors.CalcMolFormula(self.rdkit_mol)
		except:
			raw_input("stop, smiles incorrect")
			return()
		#print(self.chem_form)
		#self.mol_wt = calculate_mass(self.chem_form)
		
		try: 
			self.mol_wt = mw_dict[self.chem_form]
		except KeyError:
			try: mw_dict[self.chem_form] = calculate_mass(self.chem_form)
			except: mw_dict[self.chem_form] = calculate_mass(self.chem_form.replace("+","").replace("-",""))
			self.mol_wt = mw_dict[self.chem_form]
		
		for possible_ion in possible_ions:
			#print(possible_ion)
			possible_ion_temp = possible_ion.lstrip("[M")
			type_ion = possible_ion_temp[-1]
			possible_ion_temp = possible_ion_temp[:-1]
			
			subtract = True
			chem_form = ""
			ion_mw = self.mol_wt
			for ion_char in possible_ion_temp:
				if ion_char == "+" or ion_char == "-" or ion_char == "]":
					if ion_char == "-": subtract = True
					elif ion_char == "+": subtract = False
					if len(chem_form) > 0:
						if chem_form == "H":
							chem_form = "H1"
						if subtract: 
							try: 
								ion_mw -= mw_dict[chem_form]
							except KeyError:
								mw_dict[chem_form] = calculate_mass(chem_form)
								ion_mw -= mw_dict[chem_form]
						else:
							try: 
								ion_mw += mw_dict[chem_form]
							except KeyError:
								mw_dict[chem_form] = calculate_mass(chem_form)
								ion_mw += mw_dict[chem_form]
					chem_form = ""
					continue
				else:
					chem_form += ion_char
			#print("type_ion:")
			#print(type_ion)
			
			if type_ion == "+":	self.ions_mw_pos[possible_ion] = ion_mw
			elif type_ion == "-": self.ions_mw_neg[possible_ion] = ion_mw
			else: print("Not an ion type... %s" % (type_ion))

	def get_lipid_name(self):
		ret_names = []
		for rule in self.calc_name:
			if self.get_sub_name("head") not in rule[1]: continue
			start_reading = False
			start_reading_count = False
			partial_name = ""
			partial_name_count = ""
			temp_name = ""
			for name_char in rule[0]:

				if name_char == "}":
					#print(partial_name)
					
					if len(partial_name_count) == 0:
						temp_name += self.get_sub_name(partial_name)
					elif len(partial_name) == 0:
						if partial_name_count in self.modifications.keys():
							temp_name += str(self.get_count_smiles(self.modifications[partial_name_count]))
						else:
							temp_name += str(self.get_count_group(partial_name_count))
					else:
						temp_name += str(self.get_count_smiles_group(partial_name,partial_name_count))

					start_reading = False
					start_reading_count = False
					partial_name = ""
					partial_name_count = ""
					
				elif start_reading:
					if start_reading_count:
						partial_name_count += name_char
					elif name_char != "#":
						partial_name += name_char
					else:
						start_reading_count = True
				elif name_char == "{":
					start_reading = True
				else:
					temp_name += name_char
			#print(self.modifications)
			ret_names.append(temp_name)
		return(ret_names)

	def get_count_group(self,group_name_sub):
		#print(group_name_sub)
		#print(self.sub_names)
		return(len([sng for sng in self.sub_names if sng == group_name_sub]))

	def get_count_smiles(self,count_char):
		return(self.smiles.count(count_char))

	def get_count_smiles_group(self,group_name_sub,count_char):
		tot = 0
		name_count = [smiles for sng,smiles in zip(self.sub_names_groups,self.sub_smiles) if sng == group_name_sub]
		for nc in name_count:
			tot += nc.count(count_char)
		return(tot)

	def get_sub_name(self,group_name):
		return([sn for sng,sn in zip(self.sub_names_groups,self.sub_names) if sng == group_name][0])

	def get_smiles_sub(self,name):
		# todo change code to following:
		# [a for a,b,c in zip(l1,l2,l3) if name == b or name == c]
		smiles_indexes = [index_sn for index_sn,sn in enumerate(self.sub_names) if sn == name]
		smiles_indexes.extend([index_sn for index_sn,sn in enumerate(self.sub_names_groups) if sn == name])
		if smiles_indexes == []: return("")
		smiles_sub = "".join(operator.itemgetter(*smiles_indexes)(self.sub_smiles))
		return(smiles_sub)

	def get_molecular_formula_sub(self,name):
		# todo change code to following:
		# [a for a,b,c in zip(l1,l2,l3) if name == b or name == c]
		smiles_indexes = [index_sn for index_sn,sn in enumerate(self.sub_names) if sn == name]
		smiles_indexes.extend([index_sn for index_sn,sn in enumerate(self.sub_names_groups) if sn == name])
		if smiles_indexes == []: return("")
		smiles_sub = "".join(operator.itemgetter(*smiles_indexes)(self.sub_smiles))
		return(Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles_sub)))

	def group_in_lipid(self,name):
		smiles_indexes = [index_sn for index_sn,sn in enumerate(self.sub_names) if sn == name]
		smiles_indexes.extend([index_sn for index_sn,sn in enumerate(self.sub_names_groups) if sn == name])
		if len(smiles_indexes) > 0: return(True)
		else: return(False)


if __name__ == "__main__":
	db_in = ConfigDB(config_file="db_parse_noether.txt")
	db_in.create_db()