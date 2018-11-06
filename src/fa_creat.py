from collections import Counter
import itertools
import exrex

unique = set()
unique_fa = []
fa_list = list(exrex.generate('(C){0,9}(C|C=C|(CC\(O\)C)|(CC\(=O\)C)|(CC\(OOH\)C)){1,6}(C){0,9}(C|C=O|C\(O\)=O)'))
print(len(fa_list))
for x in fa_list:
	if str(Counter(x)) not in unique:
		if Counter(x)["C"] < 5:
			print(x)
		unique.add(str(Counter(x)))
		unique_fa.append(x)
#print(unique_fa)
print(len(unique))

unique = set()
fa_list = list(exrex.generate('(C){5,10}(C|(C=C)){1,6}(C){4,9}'))
for x in fa_list:
	if str(Counter(x)) not in unique:
		unique.add(str(Counter(x)))
print(len(unique))

list(itertools.product(*l))