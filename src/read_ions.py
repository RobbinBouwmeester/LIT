infile = open("LipidBlast-pos.msp")
all_ions = []
for line in infile:
	if "Name: " not in line: continue
	line = line.strip()
	ion = line.split(";")[1].lstrip()
	if "[" not in ion: continue
	if not ion.endswith("+"): continue
	all_ions.append(ion)
print(list(set(all_ions)))
	#print(line)

infile = open("LipidBlast-neg.msp")
ms2_ions = []
ignore = False

for line in infile:
	if "Name: " in line:
		if "; GP" in line:
			ignore = False
		else:
			ignore = True
	if ignore: continue
	if ":" in line: continue

	line = line.strip()
	ion = line.split(" ")[-1].lstrip().replace("\"","")
	#if "[" not in ion: continue
	#if not ion.endswith("+"): continue
	print(ion)
	ms2_ions.append(ion)
	#raw_input()
print(list(set(ms2_ions)))
	#print(line)