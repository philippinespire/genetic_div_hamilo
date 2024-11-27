import momi
import numpy as np
from autograd.numpy import log

sfsfile="/home/rdclark/PIRE/gazza_minuta/momi2/gmi.all.A.nohighhet.Ham.sfs.gz"
sfs = momi.Sfs.load(sfsfile)
NeConstant=1e4

#specify model
model_inf_2changeexpg_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2changeexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2changeexpg_contemp.add_size_param("n_hist")
model_inf_2changeexpg_contemp.add_size_param("n_alb")
model_inf_2changeexpg_contemp.add_size_param("n_cont")
model_inf_2changeexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_contemp.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_contemp.set_size("CBat",g=lambda params: log(params.n_alb/params.n_hist)/params.t_exp, t= "t_bot")
model_inf_2changeexpg_contemp.set_size("CBat", g=0, t="t_exp")
#run model
model_inf_2changeexpg_contemp.optimize(method="TNC")

n_bootstraps=100
model_inf_2changeexpg_copy = model_inf_2changeexpg_contemp.copy() 
bootstrap_2changeexpg_contemp = [] 
for i in range(n_bootstraps):
	# resample the data
	resampled_sfs = sfs.resample()
	print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for 2change (expg) + contemp")
	# tell model to use the new dataset
	model_inf_2changeexpg_copy.set_data(resampled_sfs, length=467359)
	# choose new random parameters for submodel, optimize
	model_inf_2changeexpg_copy.set_params(randomize=True)
	model_inf_2changeexpg_copy.optimize()
	# append results
	bootstrap_2changeexpg_contemp.append(model_inf_2changeexpg_copy.get_params()) 

boot="/home/rdclark/PIRE/gazza_minuta/momi2/Gmi_contemp2changeexpg_bootstraps.csv" 
f = open(boot,"a") 
f.write("n_hist") 
for i in range(len(bootstrap_2changeexpg_contemp)):
	f.write(",")
	f.write(str(bootstrap_2changeexpg_contemp[i].get('n_hist'))) 
f.write('\n') 
f.write("n_alb") 
for i in range(len(bootstrap_2changeexpg_contemp)):
	f.write(',')
	f.write(str(bootstrap_2changeexpg_contemp[i].get('n_alb'))) 
f.write('\n') 
f.write("n_cont") 
for i in range(len(bootstrap_2changeexpg_contemp)):
	f.write(',')
	f.write(str(bootstrap_2changeexpg_contemp[i].get('n_cont'))) 
f.write('\n') 
f.write("t_exp") 
for i in range(len(bootstrap_2changeexpg_contemp)):
	f.write(',')
	f.write(str(bootstrap_2changeexpg_contemp[i].get('t_exp'))) 
f.write('\n') 
f.write("t_bot") 
for i in range(len(bootstrap_2changeexpg_contemp)):
	f.write(',')
	f.write(str(bootstrap_2changeexpg_contemp[i].get('t_bot')))
f.close()
