#!/usr/bin/python3

import sys
import math
import glob
import os.path
import decimal

import pandas as pd

# function to read result files from directory
# Parameters:
#  - _file 				file path to read
# Returns:
#  - instance			list of instance parameters
#  - df.mean()		DF with means of key values
def read_file_text(_file):
	
	instance = _file.rpartition("/")[-1].rpartition(".txt")[0].split("_")
	instance[0] = True if instance[0] == "true" else False
	instance[1] = int(instance[1])
	instance[4] = True if instance[4] == "true" else False
	instance[5] = float(instance[5]) if instance[5] != '' else -100
	instance[6] = float(instance[6]) if instance[6] != '' else -100
	df = pd.read_csv(_file, sep=";")
	
	return instance, df.mean()

def read_file_param(_inst,_wait,_ameth,_rmeth,_bundle,_alpha,_beta):
	
	fstem = "../results"
	inst = "true" if _inst else "false"
	bundle = "true" if _bundle else "false"
	alpha = "{:.2f}".format(_alpha) if _alpha > -99 else ''
	beta = "{:.2f}".format(_beta) if _beta > -99 else ''
	fname = inst+"_"+str(_wait)+"_"+_ameth+"_"+_rmeth+"_"+bundle+"_"+alpha+"_"+beta
	instance = [_inst,_wait,_ameth,_rmeth,_bundle,_alpha,_beta]
	df = pd.read_csv(fstem+'/'+fname+'.txt', sep=";")
	
	return instance, df

def read_all(res_path="../results"):
	
	for ifile in glob.glob(res_path+"/"+"*.txt",recursive=False):
		instance, df = read_file_text(ifile)
		print(df,instance)
		#inst,dff = read_file_param(*instance)
		#print(df.mean(),instance)
	
	return




# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================
# =================================================================================




# main ==============================================
if __name__ == '__main__':
	if ( (len(sys.argv) < 1) or (len(sys.argv) > 3) ):
		print("Usage: read_result.py file [-v]")
		print("\t- file\t\t\t instance file")
		print("\t- [-v]\t\t\t verbose mode\t\t (optional)")
		sys.exit()
	
	read_all()
