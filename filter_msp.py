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

def filter_msp(infile_n,outfile_n,outfile_formula_n,max_db=10,max_c=30,min_c=2,max_mods=3,even=True,ignore_mods=True,filter_class=False,accepted_classes=["PC","PG","PA","PS","PE","PI","TG"]):
	if filter_class:
		accepted_classes.extend(["L"+ac for ac in accepted_classes])

	infile = open(infile_n).readlines()
	outfile = open(outfile_n,"w")
	outfile_formula = open(outfile_formula_n,"w")

	ignore_entry = False
	analyzed = set()
	unique_form = set()
	form_to_names = dict()
	ion_mz_dict = dict()

	counter_analyzed = 0
	counter_accepted = 0

	for line in infile:
		line = line.strip()
		#if ":" in line:
		if line.startswith("Name:"):
			name = line.replace("Name: ","")

			#Todo: multiple fatty acid tails
			fa = line.split("(")[1].split(")")[0]
			mods = line.split("(")[2].split(")")[0]
			#if name.startswith("PC") and "|[M-H]-" in name: 
		#		ignore_entry = True
		#		continue
			# TODO check for FA3!
			try: fa1,fa2 = fa.split("_")
			except: fa1,fa2,fa3 = fa.split("_")
			fa1_c,fa1_db = map(int,fa1.split(":"))
			fa2_c,fa2_db = map(int,fa2.split(":"))
			mods = sum(map(int,list(mods)))

			counter_analyzed += 1
			print(name)
			if filter_class and len([ac for ac in accepted_classes if name.startswith(ac)]) == 0:
				print("filtered on class")
				ignore_entry = True
				continue
			if fa1_c > max_c or fa2_c > max_c:
				print("filtered on max length c")
				ignore_entry = True
				continue
			if fa1_db > max_db or fa2_db > max_db:
				print("filtered on max db")
				ignore_entry = True
				continue
			if mods > max_mods:
				print("filtered on max mods")
				ignore_entry = True
				continue
			if (fa1_c < min_c or fa2_c < min_c) and mods == 0:
				print("filtered on min c")
				ignore_entry = True
				continue
			if ignore_mods:
				if fa1_c % 2 != 0 and mods == 0 and even:
					print("filtered on even")
					ignore_entry = True
					continue
				if fa2_c % 2 != 0 and mods == 0 and even:
					print("filtered on even")
					ignore_entry = True
					continue
			if not ignore_mods:
				if fa1_c % 2 != 0 and even:
					ignore_entry = True
					print("filtered on even")
					continue
				if fa2_c % 2 != 0 and even:
					ignore_entry = True
					print("filtered on even")
					continue
			if name in analyzed:
				print("filtered on prev analyzed")
				ignore_entry = True
				continue
			ignore_entry = False
			if ignore_mods and mods == 0:
				analyzed.add(name.split("(")[0]+"("+fa2+"_"+fa1+")"+line.split(")")[1])
			elif not ignore_mods:
				analyzed.add(name.split("(")[0]+"("+fa2+"_"+fa1+")"+line.split(")")[1])

			counter_accepted += 1

			print(counter_analyzed,counter_accepted,len(unique_form))

		if not ignore_entry:
			
			if line.startswith("PRECURSORMZ"):
				ion_mz = float(line.split(" ")[1])
			if line.startswith("Ion:"):
				ion = line.split(": ")[1]
				if name.startswith("PC") and ion == "[M-H]-":
					ignore_entry = True
					continue
			if line.startswith("Formula:"):
				form = line.split(": ")[1]
				ion_form = ion_chem_form(form,ion)
				if ion_form not in unique_form:
					form_to_names[ion_form] = []
					unique_form.add(ion_form)
					print(len(unique_form))
				form_to_names[ion_form].append(name)
				ion_mz_dict[ion_form] = ion_mz

			outfile.write("%s\n" % (line))

	for k in form_to_names.keys():
		outfile_formula.write("%s,%s,%s\n" % (k,ion_mz_dict[k],form_to_names[k]))
	outfile_formula.close()
	outfile.close()

if __name__ == "__main__":
	filter_msp("db_noether.msp","filtered_db_noether.msp","lipidname_to_chemicalformula_noether.csv")