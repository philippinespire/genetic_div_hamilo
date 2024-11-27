import momi
import numpy as np
from autograd.numpy import log

sfsfile="/home/rdclark/PIRE/gazza_minuta/momi2/gmi.all.A.nohighhet.Ham.sfs.gz"
sfs = momi.Sfs.load(sfsfile)
NeConstant=1e4

#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_3change_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, pre-alb size, historic size, times of three size changes
model_inf_3change_temporal.add_size_param("n_hist")
model_inf_3change_temporal.add_size_param("n_rec")
model_inf_3change_temporal.add_size_param("n_alb")
model_inf_3change_temporal.add_size_param("n_cont")
model_inf_3change_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3change_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_3change_temporal.add_leaf("CBat",N="n_cont")
model_inf_3change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_3change_temporal.set_size("CBat", N="n_alb", t="t_bot")
model_inf_3change_temporal.set_size("AHam",N="n_rec",t="t_rec")
model_inf_3change_temporal.set_size("AHam",N="n_hist",t="t_exp")
model_inf_3change_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_3change_temporal.optimize(method="TNC")

n_bootstraps=100
model_inf_3change_copy = model_inf_3change_temporal.copy() 
bootstrap_3change_temporal = [] 
for i in range(n_bootstraps):
	# resample the data
	resampled_sfs = sfs.resample()
	print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for 3change + temporal")
	# tell model to use the new dataset
	model_inf_3change_copy.set_data(resampled_sfs, length=467359)
	# choose new random parameters for submodel, optimize
	model_inf_3change_copy.set_params(randomize=True)
	model_inf_3change_copy.optimize()
	# append results
	bootstrap_3change_temporal.append(model_inf_3change_copy.get_params()) 

boot="/home/rdclark/PIRE/gazza_minuta/momi2/Gmi_temporal3change_bootstraps.csv" 
f = open(boot,"a") 
f.write("n_hist") 
for i in range(len(bootstrap_3change_temporal)):
	f.write(",")
	f.write(str(bootstrap_3change_temporal[i].get('n_hist'))) 
f.write('\n') 
f.write("n_rec")
for i in range(len(bootstrap_3change_temporal)):
        f.write(",")
        f.write(str(bootstrap_3change_temporal[i].get('n_rec')))
f.write('\n')
f.write("n_alb") 
for i in range(len(bootstrap_3change_temporal)):
	f.write(',')
	f.write(str(bootstrap_3change_temporal[i].get('n_alb'))) 
f.write('\n') 
f.write("n_cont") 
for i in range(len(bootstrap_3change_temporal)):
	f.write(',')
	f.write(str(bootstrap_3change_temporal[i].get('n_cont'))) 
f.write('\n') 
f.write("t_exp") 
for i in range(len(bootstrap_3change_temporal)):
	f.write(',')
	f.write(str(bootstrap_3change_temporal[i].get('t_exp'))) 
f.write('\n') 
f.write("t_rec")
for i in range(len(bootstrap_3change_temporal)):
        f.write(',')
        f.write(str(bootstrap_3change_temporal[i].get('t_rec')))
f.write('\n')
f.write("t_bot") 
for i in range(len(bootstrap_3change_temporal)):
	f.write(',')
	f.write(str(bootstrap_3change_temporal[i].get('t_bot')))
f.close()
