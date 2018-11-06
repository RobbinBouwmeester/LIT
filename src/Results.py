import matplotlib
matplotlib.use('TkAgg')

import pandas as pd
from scipy.stats import rankdata
import matplotlib.pyplot as plt
import numpy as np
from time import time



class IdentificationResults():
	def __init__(self):
		self.scan_res_ms = {}
		self.aligned_peaks_dict = {}

	def add_result(self,file_name,scan_identifier_ms2,scan_identifier_ms1,lipid_identifier,score_dict,aligned_peaks,db_entry,rt=0.0):
		key_res = "%s|%s|%s|%s" % (file_name,scan_identifier_ms1,scan_identifier_ms2,rt)
		
		# TODO clean this up... Just test if dictionary is there or not...
		try: self.lipid_identifier_to_entry.keys()
		except:	self.lipid_identifier_to_entry = {}
		self.lipid_identifier_to_entry[lipid_identifier] = db_entry

		if key_res in self.scan_res_ms.keys():
			self.scan_res_ms[key_res][lipid_identifier] = score_dict
			self.aligned_peaks_dict[key_res][lipid_identifier] = aligned_peaks
		else:
			self.scan_res_ms[key_res] = {}
			self.aligned_peaks_dict[key_res] = {}
			self.aligned_peaks_dict[key_res][lipid_identifier] = aligned_peaks
			self.scan_res_ms[key_res][lipid_identifier] = score_dict
			
	def get_aligned_peaks(self,spec_ident,lipid_ident):
		return(self.aligned_peaks_dict[spec_ident][lipid_ident])

	# TODO Make seperate function? Or at least connect with object
	def plot_aligned_spec(self,
						  theoretical_spectrum_mz,
						  theoretical_spectrum_intens,
						  experimental_spectrum_mz,
						  experimental_spectrum_intens,
						  fname,
						  annotations=[],
						  figwidth=20,
						  figheigth=15,
						  nomalize_exp=True):
		t0 = time()
		if nomalize_exp: experimental_spectrum_intens = [i/max(experimental_spectrum_intens) for i in experimental_spectrum_intens]
		print("time to normalize: %s" % (time()-t0))

		t0 = time()
		plt.figure(figsize=(figwidth, figheigth),dpi=250)
		plt.grid(True,ls="solid")
		plt.scatter(experimental_spectrum_mz,experimental_spectrum_intens,color="royalblue")
		plt.scatter(theoretical_spectrum_mz,[(esi/100.0)*-1 for esi in theoretical_spectrum_intens],color="Tomato")
		print("time to scatter: %s" % (time()-t0))

		t0 = time()
		plt.scatter([0],[-2.0],color="white")
		plt.yticks(np.arange(-1.0, 1.01, 0.2))
		print("time to add y-ticks: %s" % (time()-t0))

		t0 = time()
		max_mz = max([max(experimental_spectrum_mz),max(theoretical_spectrum_mz)])+50.0
		plt.plot((0.0,max_mz), (0.0,0.0), 'k-',color="black")

		ax = plt.gca()

		ax.vlines(experimental_spectrum_mz,[0.0]*len(experimental_spectrum_intens),experimental_spectrum_intens)
		ax.vlines(theoretical_spectrum_mz,[0.0]*len(theoretical_spectrum_intens),[(tsi/100.0)*-1 for tsi in theoretical_spectrum_intens])
		print("time to plot vertical lines: %s" % (time()-t0))
		
		#for x1,y1 in zip(,experimental_spectrum_intens):
		#	plt.plot((x1, x1), (0.0,y1), 'k-',color="black")
		#
		#for x1,y1 in zip(theoretical_spectrum_mz,theoretical_spectrum_intens):
		#	plt.plot((x1, x1), (0.0,(y1/100.0)*-1), 'k-',color="black")
		#

		t0 = time()
		if len(annotations) > 0:
			for mz,text in annotations:
				ax.annotate(text, xy=(mz, -1.05), xytext=((mz/(max_mz+max_mz*0.02))+(figwidth/2000),0.2),
						arrowprops=dict(facecolor='grey', shrink=0.01, width = 0.01,headwidth = 0),
						textcoords="axes fraction", rotation=90)
		print("time to add annotation: %s" % (time()-t0))
		plt.xlabel("m/z")
		plt.ylabel("Normalized intensity")

		#plt.axvline(experimental_spectrum_mz,color="blue") #,ymax=experimental_spectrum_intens
		#plt.axvline(theoretical_spectrum_mz,color="red") #,ymax=[(esi/100.0)*-1 for esi in theoretical_spectrum_intens]

		t0 = time()
		plt.savefig("%s" % (fname), format='pdf') #, bbox_inches='tight'
		plt.close()
		print("time to save and close: %s" % (time()-t0))
		
	def get_all_scores(self,score_name,unique_score_per_scan=True):
		ret_list = []
		for scan_id,scan_score_dict in self.scan_res_ms.items():
			if unique_score_per_scan: temp_scores = list(set([scan_score_dict[k][score_name] for k in scan_score_dict.keys()]))
			else: temp_scores = [scan_score_dict[k][score_name] for k in scan_score_dict.keys()]
			ret_list.extend(temp_scores)
		return(ret_list)

	def get_best_ranked_score(self,inverse_score=["avg_error_ms2","error_ms1"],ignore_score=[]):
		ret_list = []
		ranking_lipid_ident = {}
		for scan_id,scan_score_dict in self.scan_res_ms.items():
			possible_score_values = scan_score_dict[list(scan_score_dict.keys())[0]].keys()
			rows = []
			row_names = []
			for name,row in scan_score_dict.items():
				rows.append(row.values())
				row_names.append(name)
			res_df = pd.DataFrame(rows)
			res_df.columns = possible_score_values
			res_df.index = row_names

			rank_df = []
			for column in res_df:
				if column in ignore_score: continue
				if column in inverse_score: rank_df.append(rankdata(1.0/res_df[column],method="dense"))
				else: rank_df.append(rankdata(res_df[column],method="dense"))
			rank_df = pd.DataFrame(rank_df).transpose()
			rank_df.index = row_names
			ret_list.append([scan_id,(rank_df.sum(axis=1).idxmax())])
			ranking_lipid_ident[scan_id] = list(enumerate(rank_df.sum(axis=1).sort_values(ascending=False).index,1))
		#print(ranking_lipid_ident)
		return(ret_list,ranking_lipid_ident)
			#[scan_score_dict.keys()[0]]
		#
		

	def write_res_to_file(self,outfile_name="out.csv",best=False,order_score=[]):
		#"%s|%s|%s|%s" % (file_name,scan_identifier_ms1,scan_identifier_ms2,rt)
		colnames = ["exp_name","ms1_scan","ms2_scan","rt","lipid_ident","lipid_class","lipid_ion","score_rank"]
		colnames.extend(order_score)
		df_res = []

		res,ranked_lipids = self.get_best_ranked_score()

		if not best:
			for key_exp in ranked_lipids.keys():
				for lipid_rank in ranked_lipids[key_exp]:
					lipid_identifier = lipid_rank[1]
					temp_row = key_exp.split("|")
					temp_row.append(lipid_identifier.split("|")[0])
					temp_row.append(lipid_identifier.split("|")[0].split("(")[0])
					temp_row.append(lipid_identifier.split("|")[1])
					temp_row.append(lipid_rank[0])
					if len(order_score) == 0:
						order_score = list(self.scan_res_ms[key_exp][lipid_identifier].keys())
						colnames.extend(order_score)

					for score_id in order_score:
						temp_row.append(self.scan_res_ms[key_exp][lipid_identifier][score_id])
					df_res.append(temp_row)
		if best:
			for key_exp in ranked_lipids.keys():
				lipid_identifier = ranked_lipids[key_exp][0][1]
				temp_row = key_exp.split("|")
				temp_row.append(lipid_identifier.split("|")[0])
				temp_row.append(lipid_identifier.split("|")[0].split("(")[0])
				temp_row.append(lipid_identifier.split("|")[1])
				temp_row.append(1)
				if len(order_score) == 0:
					order_score = list(self.scan_res_ms[key_exp][lipid_identifier].keys())
					colnames.extend(order_score)

				for score_id in order_score:
					temp_row.append(self.scan_res_ms[key_exp][lipid_identifier][score_id])

				df_res.append(temp_row)

		df_res = pd.DataFrame(df_res,columns=colnames)
		df_res.to_csv(outfile_name)