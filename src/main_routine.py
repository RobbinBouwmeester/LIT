import matplotlib
import datetime
matplotlib.use('TkAgg')
import pickle
import os
import sys
from experiment import Experiment
from lpb import LipidBLAST
from dbs import DBSelf
from dbs import determine_chunks_prec
from isotope import Isotopic_dists,isotope_filters
import pandas as pd
import logging
from search_engine import PrecursorFilter
from search_engine import MS2Filter
from search_engine import IsotopeAlignment
from create_acyls import get_acyl_dict
from collections import Counter
from scoring_function import ScoringFunction
# TODO remove capital
from Results import IdentificationResults
import matplotlib.pyplot as plt
from os import listdir
from numpy import median
import numpy as np
from scipy.signal import find_peaks_cwt
from joblib import Parallel, delayed
from multiprocessing import Pool
from time import time
import pickle
import os.path

# Make a config file
# https://wiki.python.org/moin/ConfigParserExamples

def plot_aligned_spec(theoretical_spectrum_mz,
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
    #print("time to normalize: %s" % (time()-t0))

    #t0 = time()
    plt.figure(figsize=(figwidth, figheigth),dpi=250)
    plt.grid(True,ls="solid")
    plt.scatter(experimental_spectrum_mz,experimental_spectrum_intens,color="royalblue")
    plt.scatter(theoretical_spectrum_mz,[(esi/100.0)*-1 for esi in theoretical_spectrum_intens],color="Tomato")
    #print("time to scatter: %s" % (time()-t0))

    #t0 = time()
    plt.scatter([0],[-2.0],color="white")
    plt.yticks(np.arange(-1.0, 1.01, 0.2))
    #print("time to add y-ticks: %s" % (time()-t0))

    #t0 = time()
    max_mz = max([max(experimental_spectrum_mz),max(theoretical_spectrum_mz)])+50.0
    plt.plot((0.0,max_mz), (0.0,0.0), 'k-',color="black")

    ax = plt.gca()

    ax.vlines(experimental_spectrum_mz,[0.0]*len(experimental_spectrum_intens),experimental_spectrum_intens)
    ax.vlines(theoretical_spectrum_mz,[0.0]*len(theoretical_spectrum_intens),[(tsi/100.0)*-1 for tsi in theoretical_spectrum_intens])
    #print("time to plot vertical lines: %s" % (time()-t0))
    
    #for x1,y1 in zip(,experimental_spectrum_intens):
    #    plt.plot((x1, x1), (0.0,y1), 'k-',color="black")
    #
    #for x1,y1 in zip(theoretical_spectrum_mz,theoretical_spectrum_intens):
    #    plt.plot((x1, x1), (0.0,(y1/100.0)*-1), 'k-',color="black")
    #

    #t0 = time()
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

    #t0 = time()
    plt.savefig("%s" % (fname), format='pdf') #, bbox_inches='tight'
    plt.close()
    print("time to save and close: %s" % (time()-t0))

def moving_average(a, n=2):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def get_peaks(xic_list,peak_thres=0.1,peak_min_dist=1):
    peak_indexes = [0]
	#peak_indexes = peakutils.indexes([i for t,i in xic_list], thres=peak_thres, min_dist=peak_min_dist)
    #print([peak_indexes])
    ret_list_peak = [xic_list[idx] for idx in peak_indexes]
    #ret_list_peak = [[coord[0],coord[1]] for index,coord in enumerate(xic_list) if index in peak_indexes]
    return(ret_list_peak)
        
def read_lpb(lpb_files=["LipidBlast-pos.msp","LipidBlast-neg.msp"]):
    logging.info("Reading the LPB database ...")
    lpb = LipidBLAST(aggregate_acyls=False)
    logging.info("Done reading the LPB database ...")
    logging.info(lpb)
    return(lpb)

def read_dbs(lpb_files=["filtered_db_noether_backup.msp"],pos_ions=True,neg_ions=True,lower_limit_mz=0.0,upper_limit_mz=10000.0):
    logging.info("Reading the database ...")
    dbs = DBSelf(f_names=lpb_files,aggregate_acyls=False,pos_ions=pos_ions,neg_ions=neg_ions,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz)
    logging.info("Done reading the database ...")
    logging.info(dbs)
    print(dbs)
    return(dbs)

def read_rt_file(infile):
    infile = open(infile)
    rt_dict = {}
    for line in infile:

        if line.startswith("IDENTIFIER"): continue
        split_line = line.strip().split(",")
        rt_dict[split_line[0]] = float(split_line[1])
    return(rt_dict)

def get_script_path():
    return(os.path.dirname(os.path.realpath(sys.argv[0])))
    
def perform_search(mzml_file,exp,lower_limit_mz,upper_limit_mz,scoring_funct,
                    ms1_error=5,ms2_error=5,ms_tol_ppm=True,fa_spec=[],head_spec=[],
                    score_ms1_error=True,score_ms2_error=True,score_hypergeom=True,
                    score_intensity_expl=True,score_pred_tr=True,
                    filter_fa_spec=True,filter_head_spec=True,filter_intensity_explained=True,
                    filter_hypergeom=False,min_intensity_explained=0.1,min_hypergeom=0.0,
                    search_negative_ions=True,search_positive_ions=True,db_file="filtered_db_noether_backup.msp"):

    print("=======================")
    print(lower_limit_mz,upper_limit_mz)
    print("=======================")

    # TODO save them in a temp folder in the execution folder
    if search_positive_ions:
        try:
            lpb_pos = pickle.load(open(os.path.join(get_script_path(),"temp/lpb_pos_%s_%s.pickle" % (lower_limit_mz,upper_limit_mz)),"rb"))
        except IOError:
            lpb_pos = read_dbs(lpb_files=[db_file],pos_ions=True,neg_ions=False,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz)
            pickle.dump(lpb_pos, open(os.path.join(get_script_path(),"temp/lpb_pos_%s_%s.pickle" % (lower_limit_mz,upper_limit_mz)), "wb" ))

    if search_negative_ions:
        try:
            lpb_neg = pickle.load(open(os.path.join(get_script_path(),"temp/lpb_neg_%s_%s.pickle" % (lower_limit_mz,upper_limit_mz)),"rb"))
        except IOError:
            lpb_neg = read_dbs(lpb_files=[db_file],pos_ions=False,neg_ions=True,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz)
            pickle.dump(lpb_neg, open(os.path.join(get_script_path(),"temp/lpb_neg_%s_%s.pickle" % (lower_limit_mz,upper_limit_mz)),"wb"))
    
    if ms_tol_ppm:
        if search_positive_ions: precf_pos = PrecursorFilter(lpb_pos,ppm=ms1_error)
        else: precf_pos = []
        if search_negative_ions: precf_neg = PrecursorFilter(lpb_neg,ppm=ms1_error)
        else: precf_neg = []
    else:
        if search_positive_ions: precf_pos = PrecursorFilter(lpb_pos,dalton=ms1_error)
        else: precf_pos = []
        if search_negative_ions: precf_neg = PrecursorFilter(lpb_neg,dalton=ms1_error)
        else: precf_neg = []
        

    fa_dict,fa_form_dict = get_acyl_dict()
    ms2_filt = MS2Filter(fa_dict,ppm=ms2_error)

    results = []

    for scan_num in exp.ms2:
        try:
            exp.ms2_to_ms1[scan_num]
            exp.scan_to_spectrum[exp.ms2_to_ms1[scan_num]]
        except:
            return()

        scan_num_to_precf_res = {}
        scan_num_to_precf_res[scan_num] = []
        if exp.scan_to_spectrum[scan_num].positive_scan:
            scan_num_to_precf_res[scan_num].extend(precf_pos.retrieve_entry_pre_c_mass(exp.scan_to_spectrum[scan_num].prec_mass,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz))
        
        if exp.scan_to_spectrum[scan_num].negative_scan:
            scan_num_to_precf_res[scan_num].extend(precf_neg.retrieve_entry_pre_c_mass(exp.scan_to_spectrum[scan_num].prec_mass,lower_limit_mz=lower_limit_mz,upper_limit_mz=upper_limit_mz))
        
        chem_forms = list(set([lpb_entr.chem_form for name,delta,lpb_entr in scan_num_to_precf_res[scan_num]]))

        ms2_res = ms2_filt.align_db(list(zip(exp.scan_to_spectrum[scan_num].mz_array,
                                      exp.scan_to_spectrum[scan_num].intensity_array)),
                                      scan_num_to_precf_res[scan_num])
        
        for candidate in ms2_res:
            aligned_peaks = [[b[0],a[-1]] for a,b,c in candidate[1]]
            score_dict = {}

            # TODO division in bins based on ppm
            if score_hypergeom:
                score_dict["hyper_geom_score"] = scoring_funct.calc_hypergeom(candidate[1],candidate[2],
                                                                exp.scan_to_spectrum[scan_num].mz_array,
                                                                int(exp.scan_to_spectrum[scan_num].prec_mass/0.007))
            if score_intensity_expl:
                score_dict["intensity_explained"] = scoring_funct.calc_intensity_explained(candidate[1],candidate[2],
                                                                exp.scan_to_spectrum[scan_num].intensity_array)
            if score_ms1_error: score_dict["error_ms1"] = candidate[0][1]
            if score_ms2_error: score_dict["avg_error_ms2"] = sum([ap[2] for ap in candidate[1]])/len(candidate[1])
            
            #print([k[-1] for k,l,m in candidate[1]])
            if filter_head_spec and len([k[-1] for k,l,m in candidate[1] if k[-1] in head_spec]) == 0: continue
            #print("Head_spec")
            if filter_fa_spec and len([k[-1] for k,l,m in candidate[1] if k[-1] in fa_spec]) == 0: continue
            #print("FA_spec")
            if filter_intensity_explained and score_dict["intensity_explained"] < min_intensity_explained: continue
            #print("int")
            if filter_hypergeom and score_dict["hyper_geom_score"] < min_hypergeom: continue
            #print("hyper")

            results.append([mzml_file.split("/")[-1],
                           scan_num,exp.ms2_to_ms1[scan_num],
                           candidate[0][0],
                           score_dict,
                           aligned_peaks,
                           candidate[2],
                           exp.scan_to_spectrum[exp.ms2_to_ms1[scan_num]].scan_start_time])
    
    print("++++++++++++++++++++++")
    print(lower_limit_mz,upper_limit_mz)
    print("++++++++++++++++++++++")

    return(results)

def main_routine_func(ms1_error = 5, ms2_error = 20, ms_tol_ppm = True,
        plot_XIC=False,plot_ms2=False,aggregate_results=True,
        search_negative_ions = True,search_positive_ions = True,
        score_ms1_error=True,score_ms2_error=True,score_hypergeom=True,
        score_intensity_expl=True,score_pred_tr=True,
        min_intensity_explained=5.0,min_hypergeom=2.0,
        filter_head_spec=True,filter_fa_spec=True,
        rt_file = "rt_pred/MASSTRPLAN_preds_l3_116.csv",mzml_files_loc="",
        db_file = "filtered_db_noether_backup.msp",
        n_chunks = 32, n_cores = 8, gui=None):

    if os.path.isfile(mzml_files_loc):
        mzml_files = [mzml_files_loc]
    else:
        mzml_files = [os.path.join(mzml_files_loc,file) for file in os.listdir(mzml_files_loc) if file.lower().endswith(".mzml")]
        mgf_files = [os.path.join(mzml_files_loc,file) for file in os.listdir(mzml_files_loc) if file.lower().endswith(".mgf")]
        mzml_files.extend(mgf_files)
        
    tot_chunk_count = 0

    logging.basicConfig(filename="prec_filter.log",
                        level=logging.DEBUG,
                        filemode="w",
                        format="%(levelname)s:%(created)f:%(asctime)s:%(message)s")

    chunk_ranges = determine_chunks_prec(db_file,chunks=n_chunks)
    
    head_spec = [hs.rstrip() for hs in open("head_spec.txt")]
    fa_spec = [fs.rstrip() for fs in open("fa_spec.txt").readlines()]
        
    param_dict = {"ms1_error":ms1_error,"ms2_error":ms2_error,"ms_tol_ppm":ms_tol_ppm,"search_negative_ions":search_negative_ions,
                  "search_positive_ions":search_positive_ions,"score_ms1_error":score_ms1_error,"score_ms2_error":score_ms2_error,
                  "score_hypergeom":score_hypergeom,"score_intensity_expl":score_intensity_expl,"score_pred_tr":score_pred_tr,
                  "min_intensity_explained":min_intensity_explained,"min_hypergeom":min_hypergeom,"filter_head_spec":filter_head_spec,
                  "filter_fa_spec":filter_fa_spec,"head_spec":head_spec,"fa_spec":fa_spec,"db_file":db_file}

    iso_dists = Isotopic_dists()
    iso_filt = isotope_filters()
    iso_align = IsotopeAlignment()
    scoring_funct = ScoringFunction()
    outfile_high_confid = open("high_confident_ident.csv","w")
    scan_num_to_precf_res = {}
    aggregated_results = IdentificationResults()

    xic_movingavg_length = 4
    
    for mzml_file in mzml_files:
        base_path = os.path.join(os.path.abspath(__file__),".".join(mzml_file.split(".")[:-1]))

        results = IdentificationResults()
        if not os.path.exists(base_path):
            os.makedirs(base_path)
        if not os.path.exists(os.path.join(base_path,"MS2")):
            os.makedirs(os.path.join(base_path,"MS2"))
        if not os.path.exists(os.path.join(base_path,"XIC")):
            os.makedirs(os.path.join(base_path,"XIC"))

        if mzml_file.endswith(".mzML"): exp = Experiment(mzml_file)
        elif mzml_file.endswith(".mgf"): exp = Experiment(mzml_file,mgf=True)

        pool = Pool(processes=8)
        chunk_results = {}
        
        for lower_limit_mz,upper_limit_mz in chunk_ranges:
            chunk_results[lower_limit_mz] = pool.apply_async(perform_search, args = (mzml_file,exp,lower_limit_mz,upper_limit_mz,scoring_funct),kwds=param_dict)
            
        pool.close()
        #pool.join()

        for lower_limit_mz,cr in chunk_results.items():
            tot_chunk_count += 1
            if gui: gui.progress_search.setProperty("value", ((tot_chunk_count/(float(n_chunks)*len(mzml_files)))*100))
            for chunked_res in cr.get():
                results.add_result(*chunked_res)
                
                if aggregate_results:
                    aggregated_results.add_result(*chunked_res)

        # TODO change name "res" is stupid; change to res_ranked....
        res,ranked_lipids = results.get_best_ranked_score()

        results.write_res_to_file(outfile_name=os.path.join(base_path,"ranked_res.csv"))
        results.write_res_to_file(outfile_name=os.path.join(base_path,"best_res.csv"),best=True)
        
        if plot_XIC:
            prepared_plotting_vars = []
            for n_exp,n_lip in res:
                n_lip_entry = results.lipid_identifier_to_entry[n_lip]
                mz_list_entry = [mz for mz,intens,name in n_lip_entry]
                intens_list_entry = [intens for mz,intens,name in n_lip_entry]
                
                scan_num = n_exp.split("|")[2]

                aligned_peaks = results.get_aligned_peaks(n_exp,n_lip)

                xic_list = exp.get_XIC(exp.scan_to_spectrum[scan_num].prec_mass,
                                       positive_mode=False,
                                       negative_mode=exp.scan_to_spectrum[scan_num].negative_scan)

                xic_movingavg = moving_average([y[1] for x,y in enumerate(xic_list)],n=xic_movingavg_length)

                xic_movingavg_time = [y[0] for x,y in enumerate(xic_list)][:len(xic_movingavg)]

                xic_peaks = get_peaks(list(zip(xic_movingavg_time,xic_movingavg)))

                try: closest_peak = min(xic_peaks,key=lambda x: abs(float(x[0]) - exp.scan_to_spectrum[scan_num].scan_start_time))
                except ValueError: closest_peak = exp.scan_to_spectrum[scan_num].scan_start_time

                plt.plot(np.array(xic_movingavg_time),xic_movingavg)

                try: plt.axvline(closest_peak[0],color="green")
                except: pass

                plt.axvline(exp.scan_to_spectrum[scan_num].scan_start_time,color="red")
                plt.axvline(exp.scan_to_spectrum[exp.ms2_to_ms1[scan_num]].scan_start_time,color="blue")

                # Pipes or double points not allowed by windows...
                plt.savefig(os.path.join(base_path,"XIC/%s_%s.png" % (n_exp.replace("|","+"),n_lip.split("|")[0].replace(":","-"))))
                plt.close()

                prepared_plotting_vars.append([mz_list_entry,intens_list_entry,scan_num,n_exp,n_lip,aligned_peaks])
                
        if plot_ms2:
            pool = Pool(processes=8)
            chunk_results = {}
            for mz_list_entry,intens_list_entry,scan_num,n_exp,n_lip,aligned_peaks in prepared_plotting_vars:
                chunk_results[n_exp+n_lip] = pool.apply_async(plot_aligned_spec, args = (mz_list_entry,intens_list_entry,exp.scan_to_spectrum[scan_num].mz_array,exp.scan_to_spectrum[scan_num].intensity_array,os.path.join(base_path,"MS2/%s_%s.pdf" % (n_exp.replace("|","_").replace(":","_"),n_lip.replace("|","_").replace(":","_"))),aligned_peaks))
            pool.close()
            pool.join()
        
            for lower_limit_mz,cr in chunk_results.items():
                cr.get()

    if aggregate_results:
        aggregated_results.write_res_to_file(outfile_name=os.path.join(base_path,"aggr_ranked_res.csv"))
        aggregated_results.write_res_to_file(outfile_name=os.path.join(base_path,"aggr_best_res.csv"),best=True)

        res,ranked_lipids = aggregated_results.get_best_ranked_score()
        poss_lipids = list(set([n_lip for n_exp,n_lip in res]))
        for lip_ident in poss_lipids:
            get_rt_list = [float(n_exp.split("|")[-1]) for n_exp,n_lip in res if n_lip == lip_ident]
            if len(get_rt_list) > 3: #and np.std(get_rt_list) < 2.5:
                print(median(get_rt_list))
                rt_lip_median = median([float(n_exp.split("|")[-1]) for n_exp,n_lip in res if n_lip == lip_ident])
                outfile_high_confid.write("%s\t%s\n" % (lip_ident,rt_lip_median))
        outfile_high_confid.flush()

if __name__ == "__main__":
    main(mzml_files_loc = "C:/Users/asus/Documents/mzml/")
