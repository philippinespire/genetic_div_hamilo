import momi
import numpy as np
from autograd.numpy import log

sfsfile="/home/rdclark/PIRE/leiognathus_leuciscus/momi2/lle.all.ela.nohighhet.sfs.gz"
sfs = momi.Sfs.load(sfsfile)
NeConstant=1e4

#specify model
model_inf_2histchange_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2histchange_temporal.add_size_param("n_hist")
model_inf_2histchange_temporal.add_size_param("n_rec")
model_inf_2histchange_temporal.add_size_param("n_cont")
model_inf_2histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_2histchange_temporal.add_leaf("AHam",N="n_cont",t=109)
model_inf_2histchange_temporal.move_lineages("CNas","AHam",t=110)
model_inf_2histchange_temporal.set_size("AHam", N="n_hist", t="t_exp")
model_inf_2histchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_2histchange_temporal.optimize(method="TNC")

n_bootstraps=100
model_inf_2histchange_copy = model_inf_2histchange_temporal.copy() 
bootstrap_2histchange_temporal = [] 
for i in range(n_bootstraps):
	# resample the data
	resampled_sfs = sfs.resample()
	print(f"Fitting {i+1}-th bootstrap out of {n_bootstraps} for 2histchange + temporal")
	# tell model to use the new dataset
	model_inf_2histchange_copy.set_data(resampled_sfs, length=485532)
	# choose new random parameters for submodel, optimize
	model_inf_2histchange_copy.set_params(randomize=True)
	model_inf_2histchange_copy.optimize()
	# append results
	bootstrap_2histchange_temporal.append(model_inf_2histchange_copy.get_params()) 

boot="/home/rdclark/PIRE/leiognathus_leuciscus/momi2/Ela_temp2histchange_bootstraps.csv" 
f = open(boot,"a") 
f.write("n_hist") 
for i in range(len(bootstrap_2histchange_temporal)):
	f.write(",")
	f.write(str(bootstrap_2histchange_temporal[i].get('n_hist'))) 
f.write('\n') 
f.write("n_rec") 
for i in range(len(bootstrap_2histchange_temporal)):
	f.write(',')
	f.write(str(bootstrap_2histchange_temporal[i].get('n_rec'))) 
f.write('\n') 
f.write("n_cont") 
for i in range(len(bootstrap_2histchange_temporal)):
	f.write(',')
	f.write(str(bootstrap_2histchange_temporal[i].get('n_cont'))) 
f.write('\n') 
f.write("t_exp") 
for i in range(len(bootstrap_2histchange_temporal)):
	f.write(',')
	f.write(str(bootstrap_2histchange_temporal[i].get('t_exp'))) 
f.write('\n') 
f.write("t_rec") 
for i in range(len(bootstrap_2histchange_temporal)):
	f.write(',')
	f.write(str(bootstrap_2histchange_temporal[i].get('t_rec')))
f.close()
