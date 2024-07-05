#!/usr/bin/python3

import sys
import math
import glob
import os.path
import decimal

import pandas as pd
import matplotlib.pyplot as plt

# generate different root colours
def get_hex_col(n):
	mid = int(n/2)
	# blue to green to red
	hsv_tuples = [(0.0, x * 2.0 / n, 1- x * 2.0 / n ) for x in range(mid)]
	hsv_tuples = hsv_tuples + [(x * 2.0 / n, 1- x * 2.0 / n, 0.0 ) for x in range(n-1-mid)]
	hsv_tuples.append((1.0,0.0,0.0))
	# yellow to blue
	#HSV_tuples = [(1.0 - x * 1.0 / n , 1 - x * 1.0 / n, x * 1.0 / n ) for x in range(n-1)]
	#HSV_tuples.append((0.0,0.0,1.0))
	hex_out = []
	for rgb in hsv_tuples:
		rgb = tuple(int(x * 255) for x in rgb)
		hex_out.append('#%02x%02x%02x' % (rgb))
	return hex_out

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
	instance[7] = float(instance[7]) if instance[7] != '' else 1.0
	df = pd.read_csv(_file, sep=";")
	
	return instance, df.mean(), df.sem()

def read_file_param(_inst,_wait,_ameth,_rmeth,_bundle,_alpha,_beta,_scale):
	
	fstem = "./results"
	inst = "true" if _inst else "false"
	bundle = "true" if _bundle else "false"
	alpha = "{:.2f}".format(_alpha) if _alpha > -99 else ''
	beta = "{:.2f}".format(_beta) if _beta > -99 else ''
	scale = "{:.2f}".format(_scale)
	fname = inst+"_"+str(_wait)+"_"+_ameth+"_"+_rmeth+"_"+bundle+"_"+alpha+"_"+beta+"_"+scale
	instance = [_inst,_wait,_ameth,_rmeth,_bundle,_alpha,_beta,_scale]
	df = pd.read_csv(fstem+'/'+fname+'.txt', sep=";")
	
	return instance, df

def read_all(res_path="./results",force_read=True):
	if (not os.path.isfile("./eval/data_dafd.csv") or force_read):
		param_labels = ["instance", "max_wait", "a_meth", "r_meth", "bundling", "alpha", "beta", "scaling_factor"]
		with open(res_path+"/"+os.listdir(res_path)[0], "r") as file:
			print("Result \".csv\" not found -> read data from \".txt\"")
			data_labels = file.readline().strip().split(';')
		
		df_all = pd.DataFrame(columns=param_labels+data_labels)
		for ifile in glob.glob(res_path+"/"+"*.txt",recursive=False):
			instance, df_mean, _ = read_file_text(ifile)
			df_all.loc[len(df_all)] = instance+df_mean.values.tolist()
		
		#df_all = df_all.sort_values(by=["instance","a_meth","r_meth"],ascending=True)
		df_all = df_all.sort_values(by=["instance","max_wait","scaling_factor",data_labels[1]],ascending=True)
		df_all.reset_index(drop=True,inplace=True)
		df_all.to_csv("./eval/data_dafd.csv", encoding='utf-8', index=False)
	else:
		print("Result \".csv\" found -> opened")
		df_all = pd.read_csv("./eval/data_dafd.csv")
	
	return df_all

def plot_ameth_over_alpha (df,ins=False,wait=1500,bun=True,scale=1.0,yval="averageDelay"):
	
	df_tmp = df.loc[ (df["instance"] == ins) & (df["max_wait"] == wait) & (df["scaling_factor"] == scale) & (df["a_meth"] == 'w') & (df["bundling"] == bun) & (df["alpha"] < 0.999),["r_meth","alpha","beta",yval]]
	
	methods = list(df_tmp.loc[:,"r_meth"].unique())
	alphas = sorted(list(df_tmp.loc[:,"alpha"].unique()))
	betas = list(df_tmp.loc[:,"beta"].unique())
	#print(methods,alphas,betas)
	
	df_tmp_comp = df.loc[ (df["instance"] == ins) & (df["max_wait"] == wait) & (df["scaling_factor"] == scale) & (df["a_meth"] == 'n') & (df["r_meth"] == 's') & (df["bundling"] == bun),[yval]]
	static_val = df_tmp_comp[yval].iloc[0]
	
	para = []
	y_list = []
	for meth in methods:
		if meth == 'l':
			for bet in betas:
				if bet == -100:
					continue
				y_list.append(df_tmp.loc[ (df_tmp["r_meth"] == meth) & (df_tmp["beta"] == bet),["alpha",yval]].sort_values(by=["alpha"],ascending=True)[yval].tolist())
				para.append((meth,bet))
		else:
			y_list.append(df_tmp.loc[ (df_tmp["r_meth"] == meth),["alpha",yval]].sort_values(by=["alpha"],ascending=True)[yval].tolist())
			para.append((meth,-100))
	
	line_colours = get_hex_col(len(para))
	fig = plt.figure(figsize=(8, 4))
	ax = fig.add_subplot(111)
	plt.title("{} per alpha (inst={}, wait={}, bundle={}, scale={:.2f})".format(yval,ins,wait,bun,scale))
	plt.xlabel("alpha")
	plt.ylabel(yval)
	for i in range(len(para)):
		if para[i][0] == 'n': continue
		plt.plot(alphas, y_list[i], marker='x', color=line_colours[i],label=para[i])
	
	plt.hlines(y=static_val, xmin=alphas[0], xmax=alphas[-1], colors='aqua', linestyles='--', lw=2, label='static_assignment')
	
	# add legend, grid and make the figure smaller
	plt.legend(loc='best',fontsize="small")
	plt.grid(alpha=0.5, linestyle='-', linewidth=0.5,color='#cccccc')
	fig.tight_layout()#(rect=[-0.02, 0.03, 1, 0.95])
	
	#plt.show()
	plt.savefig("./eval/plots/{}/plot_{}_{}_{}_{:.2f}_{}.png".format(yval,ins,wait,bun,scale,yval))
	print("saved plot to \"./eval/plots/{}/plot_{}_{}_{}_{:.2f}_{}.png\"".format(yval,ins,wait,bun,scale,yval))
	
	plt.close()
	#sys.exit(-1)
	
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
	
	df = read_all(force_read=False)
	print(df.head())
	#print(df.dtypes)
	
	ins = [False,True]
	wai = [1200,1500,1800]
	bun = [False,True]
	scale = [0.90,0.95,1.0,1.05,1.10,1.15,1.20,1.50,1.75,2.00]
	val = ["averageDelay","percDelayed","percBundled"]
	#val = ["percDelayed"]
	#val = ["percBundled"]
	
	for inn in ins:
		for waa in wai:
			for buu in bun:
				for vaa in val:
					if not buu and vaa == "percBundled": continue
					for sca in scale:
						plot_ameth_over_alpha(df,ins=inn,wait=waa,bun=buu,scale=sca,yval=vaa)
