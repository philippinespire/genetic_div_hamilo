# Running MOMI2 on Wahab

Code for *Gazza minuta*.

## Create input files for momi2

```sh
cd ~/PIRE/gazza_minuta/momi2
salloc
module load container_env python3
bash

#make allele counts file to make SFS from
crun.python3 -p ~/.conda/envs/momi-py36 python -m momi.read_vcf --no_aa ../gmi.all.A.nohighhet.Ham.vcf.gz ind2pop.txt gmi.all.A.nohighhet.Ham.allelecounts

#make SFS
crun.python3 -p ~/.conda/envs/momi-py36 python -m momi.extract_sfs gmi.all.A.nohighhet.Ham.sfs.gz 100 gmi.all.A.nohighhet.Ham.allelecounts

#read SFS into momi and run models
crun.python3 -p ~/.conda/envs/momi-py36 python
import momi
from autograd.numpy import log

sfsfile="gmi.all.A.nohighhet.Ham.sfs.gz"
sfs = momi.Sfs.load(sfsfile)
NeConstant=1e4

#check file
print("populations", sfs.populations)
populations ('AHam', 'CBat')
print("percent missing data per population", sfs.p_missing)
percent missing data per population [0.08277804 0.03494044]
```

## Constant pop size, contemp only

```bash
#specify model
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CBat",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.2506224208032938
            jac: array([1.56363903e-14])
  kl_divergence: 0.2506224208032938
 log_likelihood: -18372.833443486383
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_constant': 16165.186753328191})
         status: 1
        success: True
              x: array([9.69061524])
```

## Constant pop size, temporal only

```bash
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_temponly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.16220847010228867
            jac: array([1.34123456e-11])
  kl_divergence: 0.16220847010228867
 log_likelihood: -9636.232159114246
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 10
            nit: 4
     parameters: ParamsDict({'n_constant': 16023.936072205768})
         status: 1
        success: True
              x: array([9.68183889])
```

## Constant pop size, temp & contemp

```bash
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=467359)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CBat",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)#adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CBat","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.7702807318652237
            jac: array([4.70280014e-13])
  kl_divergence: 0.7702807318652237
 log_likelihood: -31102.53323245625
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 4
     parameters: ParamsDict({'n_constant': 15797.57240966082})
         status: 1
        success: True
              x: array([9.66761156])
```

## Recent size change, contemp only

```bash
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_change_contemp.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_change_contemp.add_size_param("n_alb")
model_inf_change_contemp.add_size_param("n_bot")
model_inf_change_contemp.add_time_param("t_bot",upper=1e2) #force bot to be in last 100 years (gens?)
model_inf_change_contemp.add_leaf("CBat",N="n_bot")
model_inf_change_contemp.set_size("CBat", N="n_alb", t="t_bot") #says CBat pop changes from n_alb to n_bot at t_bot
#run model
model_inf_change_contemp.optimize(method="TNC")

            fun: 0.24450080626167342
            jac: array([ 6.72152988e-10, -7.53868258e-09, -4.39117672e-07])
  kl_divergence: 0.24450080626167342
 log_likelihood: -18330.190276589456
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 11
     parameters: ParamsDict({'n_alb': 15303.182266042042, 'n_bot': 10000000000.000004, 't_bot': 99.99179043069756})
         status: 1
        success: True
              x: array([ 9.63581608, 23.02585093,  9.4075429 ])
```

## Recent size change, temp and contemp

```bash
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_change_temporal.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_change_temporal.add_size_param("n_alb")
model_inf_change_temporal.add_size_param("n_bot")
model_inf_change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_change_temporal.add_leaf("CBat",N="n_bot")
model_inf_change_temporal.set_size("CBat", N="n_alb", t="t_bot")
model_inf_change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_change_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_change_temporal.optimize(method="TNC")

            fun: 0.7699491955271237
            jac: array([ 6.85826172e-06, -2.21925528e-06, -1.72950426e-04])
  kl_divergence: 0.7699491955271237
 log_likelihood: -31100.05731908332
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 17
            nit: 6
     parameters: ParamsDict({'n_alb': 15686.130012185804, 'n_bot': 1372427.516387697, 't_bot': 24.199900882073898})
         status: 1
        success: True
              x: array([ 9.66053216, 14.13209164, -1.14175106])
```

## Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_contemp.add_size_param("n_rec")
model_inf_recchange_contemp.add_size_param("n_cont")
model_inf_recchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_contemp.add_leaf("CBat",N="n_cont")
model_inf_recchange_contemp.set_size("CBat", N="n_rec", t="t_rec")
#run model
model_inf_recchange_contemp.optimize(method="TNC")

            fun: 0.23184684488892954
            jac: array([-1.76006065e-07, -2.87967495e-07, -4.10732147e-08])
  kl_divergence: 0.23184684488892954
 log_likelihood: -18242.04278166692
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 51
            nit: 18
     parameters: ParamsDict({'n_rec': 14017.756483375073, 'n_cont': 65244.97555925329, 't_rec': 499.99841934225117})
         status: 1
        success: True
              x: array([ 9.54808013, 11.08590432, 12.4134895 ])
```

## Pre-Albatross size change, temporal only

```bash
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temponly.add_size_param("n_rec")
model_inf_recchange_temponly.add_size_param("n_cont")
model_inf_recchange_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_recchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_recchange_temponly.optimize(method="TNC")

            fun: 0.1338797011026867
            jac: array([ 1.64004486e-08, -3.50646031e-08, -1.10798173e-05])
  kl_divergence: 0.1338797011026867
 log_likelihood: -9504.95664357009
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 11
     parameters: ParamsDict({'n_rec': 13195.300312442134, 'n_cont': 10000000000.000004, 't_rec': 499.7739173188368})
         status: 1
        success: True
              x: array([ 9.48761601, 23.02585093,  7.44985249])
```

## Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temporal.add_size_param("n_rec")
model_inf_recchange_temporal.add_size_param("n_cont")
model_inf_recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temporal.add_leaf("CBat",N="n_cont")
model_inf_recchange_temporal.add_leaf("AHam", N="n_cont", t=109)
model_inf_recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_recchange_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_recchange_temporal.optimize(method="TNC")

            fun: 0.7639513728921632
            jac: array([-2.62128510e-07, -4.55627005e-07, -5.17893676e-06])
  kl_divergence: 0.7639513728921632
 log_likelihood: -31055.265579645435
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 21
            nit: 9
     parameters: ParamsDict({'n_rec': 14347.544833174676, 'n_cont': 26138.37828687529, 't_rec': 499.71566147393264})
         status: 1
        success: True
              x: array([ 9.57133411, 10.17115995,  7.22043789])
```

## Historic size change, contemp only

```bash
#specify model
model_inf_histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_histchange_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_contemp.add_size_param("n_hist")
model_inf_histchange_contemp.add_size_param("n_cont")
model_inf_histchange_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_contemp.add_leaf("CBat",N="n_cont")
model_inf_histchange_contemp.set_size("CBat", N="n_hist", t="t_exp")
#run model
model_inf_histchange_contemp.optimize(method="TNC")
			    
            fun: 0.07152124356863589
            jac: array([-1.44818379e-07,  1.37754548e-07, -1.97577878e-07])
  kl_divergence: 0.07152124356863589
 log_likelihood: -17125.214642869756
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 36
            nit: 12
     parameters: ParamsDict({'n_hist': 7428.225014085398, 'n_cont': 58953.51133302024, 't_exp': 30517.659565273352})
         status: 1
        success: True
              x: array([ 8.91304221, 10.98450447, -3.8554763 ])
```

## Historic size change, temporal only

```bash
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_histchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temponly.add_size_param("n_hist")
model_inf_histchange_temponly.add_size_param("n_cont")
model_inf_histchange_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_histchange_temponly.optimize(method="TNC")
			    
            fun: 0.015552430033766038
            jac: array([ 6.06585359e-07, -2.60913932e-06,  1.96970249e-06])
  kl_divergence: 0.015552430033766038
 log_likelihood: -8956.628069436712
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 50
            nit: 16
     parameters: ParamsDict({'n_hist': 8350.7226662616, 'n_cont': 53614.480017999405, 't_exp': 27746.5696715168})
         status: 1
        success: True
              x: array([ 9.03010336, 10.88957446, -4.00342426])
```

## Historic size change, temp and contemp		  

```bash
#specify model
model_inf_histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_histchange_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temporal.add_size_param("n_hist")
model_inf_histchange_temporal.add_size_param("n_cont")
model_inf_histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_temporal.add_leaf("CBat",N="n_cont")
model_inf_histchange_temporal.add_leaf("AHam",N="n_cont",t=109)
model_inf_histchange_temporal.set_size("AHam", N="n_hist", t="t_exp") #says, at some time (t_exp) in past, Alb was at n_hist size
model_inf_histchange_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_histchange_temporal.optimize(method="TNC")
			    
            fun: 0.6160067286777614
            jac: array([1.38551981e-07, 1.18583448e-07, 7.93490083e-06])
  kl_divergence: 0.6160067286777614
 log_likelihood: -29950.414976652282
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 13
     parameters: ParamsDict({'n_hist': 7038.525502277384, 'n_cont': 53151.65545163074, 't_exp': 32581.843532502986})
         status: 1
        success: True
              x: array([ 8.85915398, 10.88090453, -3.75748463])
```
			  
## Recent and Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recchange_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recchange_contemp.add_size_param("n_rec")
model_inf_2recchange_contemp.add_size_param("n_alb")
model_inf_2recchange_contemp.add_size_param("n_cont")
model_inf_2recchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recchange_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recchange_contemp.add_leaf("CBat",N="n_cont")
model_inf_2recchange_contemp.set_size("CBat", N="n_rec", t="t_rec")
model_inf_2recchange_contemp.set_size("CBat", N="n_alb", t="t_bot")
#run model
model_inf_2recchange_contemp.optimize(method="TNC")

            fun: 0.22808909474166966
            jac: array([ 3.33333058e-09, -1.43665061e-08,  2.22180589e-09, -1.04920279e-29,
        3.93935841e-07])
  kl_divergence: 0.22808909474166966
 log_likelihood: -18215.86629414111
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 86
            nit: 19
     parameters: ParamsDict({'n_rec': 12472.287900375804, 'n_alb': 10000000000.000004, 'n_cont': 5.152213620183855, 't_rec': 500.0, 't_bot': 0.05301016451940768})
         status: 1
        success: True
              x: array([ 9.43126449, 23.02585093,  1.63942645, 62.70585525, -7.54191154])
```

## Recent and Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recchange_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recchange_temporal.add_size_param("n_rec")
model_inf_2recchange_temporal.add_size_param("n_alb")
model_inf_2recchange_temporal.add_size_param("n_cont")
model_inf_2recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recchange_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recchange_temporal.add_leaf("CBat",N="n_cont")
model_inf_2recchange_temporal.set_size("CBat", N="n_alb", t="t_bot")
model_inf_2recchange_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_2recchange_temporal.move_lineages("CBat","AHam",t=110)
model_inf_2recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_2recchange_temporal.optimize(method="TNC")

            fun: 0.761772124860525
            jac: array([ 6.81900333e-06, -3.03293166e-06,  6.23756628e-08, -1.22802079e-23,
       -1.18934371e-08])
  kl_divergence: 0.761772124860525
 log_likelihood: -31038.99095534516
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 74
            nit: 19
     parameters: ParamsDict({'n_rec': 14123.790160607798, 'n_alb': 42204.47463679947, 'n_cont': 1.0, 't_rec': 500.0, 't_bot': 0.00757433491549252})
         status: 1
        success: True
              x: array([ 9.5556159 , 10.65028153,  0.        , 47.97644655, -9.48808417])
```
			  
## Recent and historic size change, contemp only

```bash
#specify model
model_inf_2change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2change_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_contemp.add_size_param("n_hist")
model_inf_2change_contemp.add_size_param("n_alb")
model_inf_2change_contemp.add_size_param("n_cont")
model_inf_2change_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2change_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2change_contemp.add_leaf("CBat",N="n_cont")
model_inf_2change_contemp.set_size("CBat", N="n_alb", t="t_bot")
model_inf_2change_contemp.set_size("CBat", N="n_hist", t="t_exp")
#run model
model_inf_2change_contemp.optimize(method="TNC")

            fun: 0.07155496869664156
            jac: array([ 1.56511466e-05, -1.02625671e-06,  5.05112656e-07,  4.32911879e-05,
        3.31095153e-05])
  kl_divergence: 0.07155496869664156
 log_likelihood: -17125.449572111444
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 87
            nit: 19
     parameters: ParamsDict({'n_hist': 7430.315413670844, 'n_alb': 59002.98378383065, 'n_cont': 3992782.33971139, 't_exp': 30571.856558195395, 't_bot': 1.6198173716462398})
         status: 1
        success: True
              x: array([ 8.91332359, 10.98534329, 15.19999887, -3.8527824 , -4.10652598])
```

## Recent and historic size change, temp and contemp

```bash
#specify model
model_inf_2change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2change_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_temporal.add_size_param("n_hist")
model_inf_2change_temporal.add_size_param("n_alb"
model_inf_2change_temporal.add_size_param("n_cont")
model_inf_2change_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2change_temporal.add_leaf("CBat",N="n_cont")
model_inf_2change_temporal.set_size("CBat", N="n_alb", t="t_bot")
model_inf_2change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_2change_temporal.move_lineages("CBat","AHam",t=110)
model_inf_2change_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2change_temporal.optimize(method="TNC")

            fun: 0.6098680176073412
            jac: array([ 7.53811060e-06, -1.96090588e-05,  1.96743666e-07,  2.30811152e-06,
        4.44461282e-08])
  kl_divergence: 0.6098680176073412
 log_likelihood: -29904.571082378385
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 66
            nit: 17
     parameters: ParamsDict({'n_hist': 7470.72651001182, 'n_alb': 57600.18092916944, 'n_cont': 7.058370834976229, 't_exp': 31183.35057984658, 't_bot': 0.07369439268594485})
         status: 1
        success: True
              x: array([ 8.91874753, 10.96128099,  1.95421426, -3.82285985, -7.21226154])
```

## Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_contemp.add_size_param("n_hist")
model_inf_2histchange_contemp.add_size_param("n_rec")
model_inf_2histchange_contemp.add_size_param("n_cont")
model_inf_2histchange_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_contemp.add_leaf("CBat",N="n_cont")
model_inf_2histchange_contemp.set_size("CBat", N="n_rec", t="t_rec")
model_inf_2histchange_contemp.set_size("CBat", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_contemp.optimize(method="TNC")

            fun: 0.04903809782126269
            jac: array([-2.76091953e-08, -2.50565389e-07,  3.94514050e-08, -1.02442931e-06,
        4.95677806e-11])
  kl_divergence: 0.04903809782126269
 log_likelihood: -16968.597049593554
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 25
            nit: 6
     parameters: ParamsDict({'n_hist': 9440.17385285424, 'n_rec': 395397.4431814586, 'n_cont': 1126.2912424573221, 't_exp': 21086.151226208312, 't_rec': 111.00006610029787})
         status: 1
        success: True
              x: array([  9.15272968,  12.88764672,   7.02668543,  -4.48074695,
       -15.58791648])
```

## Pre-Albatross and historic size change, temporal only

```bash
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_temponly.add_size_param("n_hist")
model_inf_2histchange_temponly.add_size_param("n_rec")
model_inf_2histchange_temponly.add_size_param("n_cont")
model_inf_2histchange_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_2histchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
model_inf_2histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_temponly.optimize(method="TNC")

            fun: 0.014848650092945205
            jac: array([ 6.94627836e-08,  1.04348965e-07, -2.26065191e-07, -4.13375733e-07,
        1.69782909e-06])
  kl_divergence: 0.014848650092945205
 log_likelihood: -8953.366753190949
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 44
            nit: 11
     parameters: ParamsDict({'n_hist': 8908.866377958599, 'n_rec': 66054.05161115069, 'n_cont': 7174.523044522323, 't_exp': 24952.429841077937, 't_rec': 182.22538086298738})
         status: 1
        success: True
              x: array([ 9.09480228, 11.09822865,  8.87829156, -4.17761244, -1.49549316])
```

## Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_temporal.add_size_param("n_hist")
model_inf_2histchange_temporal.add_size_param("n_rec")
model_inf_2histchange_temporal.add_size_param("n_cont")
model_inf_2histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_temporal.add_leaf("CBat",N="n_cont")
model_inf_2histchange_temporal.add_leaf("AHam", N="n_cont",t=109)
model_inf_2histchange_temporal.move_lineages("CBat","AHam",t=110)
model_inf_2histchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_2histchange_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_temporal.optimize(method="TNC")

            fun: 0.5878142707574121
            jac: array([ 3.50810121e-07, -1.96391560e-06,  3.73177667e-06, -2.11949530e-08,
        6.26755333e-08])
  kl_divergence: 0.5878142707574121
 log_likelihood: -29739.873700903114
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 79
            nit: 19
     parameters: ParamsDict({'n_hist': 8933.336648711847, 'n_rec': 98376.30650448188, 'n_cont': 9080.811583025314, 't_rec': 499.998530931626, 't_exp': 25212.725105915804})
         status: 1
        success: True
              x: array([ 9.09754525, 11.49655527,  9.11391885, 12.48670241, -4.16008701])
```

## Recent, Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_3change_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3change_contemp.add_size_param("n_hist")
model_inf_3change_contemp.add_size_param("n_rec")
model_inf_3change_contemp.add_size_param("n_alb")
model_inf_3change_contemp.add_size_param("n_cont")
model_inf_3change_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3change_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3change_contemp.add_time_param("t_bot",upper=1e2)
model_inf_3change_contemp.add_leaf("CBat",N="n_cont")
model_inf_3change_contemp.set_size("CBat", N="n_alb", t="t_bot")
model_inf_3change_contemp.set_size("CBat", N="n_rec", t="t_rec")
model_inf_3change_contemp.set_size("CBat", N="n_hist", t="t_exp")
#run model
model_inf_3change_contemp.optimize(method="TNC")

            fun: 0.04897154741356404
            jac: array([ 5.05592248e-06,  7.89478324e-06, -1.01691424e-10,  2.31889170e-05,
       -1.55037042e-06, -1.64438288e-10, -1.56741448e-06])
  kl_divergence: 0.04897154741356404
 log_likelihood: -16968.133459453526
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 19
            nit: 4
     parameters: ParamsDict({'n_hist': 9408.296548626979, 'n_rec': 304586.0299384168, 'n_alb': 5023167.796634158, 'n_cont': 285.03675539903486, 't_exp': 21246.924405844784, 't_rec': 112.5303930572937, 't_bot': 26.066053801840262})
         status: 1
        success: True
              x: array([ 9.14934719, 12.62670885, 15.42957133,  5.65261814, -4.46618468,
       -5.53411281, -1.04253823])
```

## Recent, Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_3change_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3change_temporal.add_size_param("n_hist")
model_inf_3change_temporal.add_size_param("n_rec")
model_inf_3change_temporal.add_size_param("n_alb")
model_inf_3change_temporal.add_size_param("n_cont")
model_inf_3change_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3change_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_3change_temporal.add_leaf("CBat",N="n_cont")
model_inf_3change_temporal.add_leaf("AHam", N="n_alb",t=109)
model_inf_3change_temporal.move_lineages("CBat","AHam",t=110)
model_inf_3change_temporal.set_size("CBat", N="n_alb", t="t_bot")
model_inf_3change_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_3change_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_3change_temporal.optimize(method="TNC")

            fun: 0.5805378963305511
            jac: array([ 3.69789803e-06, -4.61075157e-06,  1.36891349e-06, -1.08278281e-07,
        4.46514099e-06,  3.57252805e-06, -1.33569786e-08])
  kl_divergence: 0.5805378963305511
 log_likelihood: -29685.533736683316
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 46
            nit: 9
     parameters: ParamsDict({'n_hist': 9553.067461841363, 'n_rec': 229800.72793804668, 'n_alb': 2825.663566342768, 'n_cont': 190245720.69705817, 't_exp': 21434.724883204864, 't_rec': 356.56399042119176, 't_bot': 99.99983417065636})
         status: 1
        success: True
              x: array([ 9.16461758, 12.34496781,  7.94649851, 19.06382706, -4.44943273,
        0.53766856, 13.30971988])
```

## Recent exponential change, contemp only 

```bash
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
#specify model
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_expg_contemp.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_expg_contemp.add_size_param("n_alb")
model_inf_expg_contemp.add_size_param("n_bot")
model_inf_expg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_expg_contemp.add_leaf("CBat",N="n_bot",g=lambda params: log(params.n_bot/params.n_alb)/params.t_bot) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_expg_contemp.set_size("CBat",g=0, t="t_bot")
#run model
model_inf_expg_contemp.optimize(method="TNC")

            fun: 0.24490886991318164
            jac: array([ 1.08661246e-08, -3.01576620e-05, -6.75259157e-06])
  kl_divergence: 0.24490886991318164
 log_likelihood: -18333.03284798586
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 12
     parameters: ParamsDict({'n_alb': 15345.839245197996, 'n_bot': 10000000000.000004, 't_bot': 99.86655348753256})
         status: 1
        success: True
              x: array([ 9.63859966, 23.02585093,  6.61788937])
```

## Recent exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data
model_inf_expg_temporal.set_data(sfs,length=467359)
#set parameterss to infer - contemp size, alb size, time of bottleneck
model_inf_expg_temporal.add_size_param("n_alb")
model_inf_expg_temporal.add_size_param("n_bot")
model_inf_expg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_expg_temporal.add_leaf("CBat",N="n_bot",g=lambda params: log(params.n_bot/params.n_alb)/params.t_bot)
model_inf_expg_temporal.set_size("CBat",g=0, t="t_bot")
model_inf_expg_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_expg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_expg_temporal.optimize(method="TNC")

            fun: 0.7702795401388605
            jac: array([-1.31307495e-05, -2.43945987e-06, -1.18905632e-06])
  kl_divergence: 0.7702795401388605
 log_likelihood: -31102.52433264377
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 18
            nit: 7
     parameters: ParamsDict({'n_alb': 15794.559120744812, 'n_bot': 24144.858091275506, 't_bot': 0.3552324987966256})
         status: 1
        success: True
              x: array([ 9.6674208 , 10.09182672, -5.63659431])
```
		  
## Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_contemp.add_size_param("n_rec")
model_inf_recexpg_contemp.add_size_param("n_cont")
model_inf_recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_contemp.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_contemp.set_size("CBat", g=0, t="t_rec")
#run model
model_inf_recexpg_contemp.optimize(method="TNC")

            fun: 0.2341330822113826
            jac: array([ 5.16159524e-11, -5.07095183e-06, -2.11773060e-08])
  kl_divergence: 0.2341330822113826
 log_likelihood: -18257.96871085513
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 47
            nit: 16
     parameters: ParamsDict({'n_rec': 15373.383935052241, 'n_cont': 10000000000.000004, 't_rec': 499.99738323812005})
         status: 1
        success: True
              x: array([ 9.64039298, 23.02585093, 11.90939027])
```

## Pre-Albatross exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recexpg_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temponly.add_size_param("n_rec")
model_inf_recexpg_temponly.add_size_param("n_cont")
model_inf_recexpg_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temponly.set_size("AHam", g=0, t="t_rec")
#run model
model_inf_recexpg_temponly.optimize(method="TNC")

            fun: 0.13571769348901108
            jac: array([ 3.76707350e-10, -1.39520679e-04, -4.36861448e-07])
  kl_divergence: 0.13571769348901108
 log_likelihood: -9513.473900288318
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 38
            nit: 14
     parameters: ParamsDict({'n_rec': 13322.38780568862, 'n_cont': 10000000000.000004, 't_rec': 499.9906308651865})
         status: 1
        success: True
              x: array([ 9.49720119, 23.02585093, 10.63388978])
```

## Pre-Albatross size change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_recexpg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temporal.add_size_param("n_rec")
model_inf_recexpg_temporal.add_size_param("n_cont")
model_inf_recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_temporal.add_leaf("CBat",N="n_cont")
model_inf_recexpg_temporal.add_leaf("AHam", N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_recexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_recexpg_temporal.optimize(method="TNC")

            fun: 0.7662874960762641
            jac: array([ 2.40866197e-06, -9.80232072e-07, -2.85806499e-07])
  kl_divergence: 0.7662874960762641
 log_likelihood: -31072.7117475843
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 16
     parameters: ParamsDict({'n_rec': 15119.580090804253, 'n_cont': 30785.114992962892, 't_rec': 499.9690238296591})
         status: 1
        success: True
              x: array([ 9.62374588, 10.33478657,  9.43803678])
```
			  
## Historic exponential change, contemp only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_histexpg_contemp.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_contemp.add_size_param("n_hist")
model_inf_histexpg_contemp.add_size_param("n_cont")
model_inf_histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_contemp.add_leaf("CBat",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_contemp.set_size("CBat",g=0, t="t_exp")
#run model
model_inf_histexpg_contemp.optimize(method="TNC")

            fun: 0.08419708272069557
            jac: array([-2.01759000e-07,  1.52776059e-06, -3.34660793e-06])
  kl_divergence: 0.08419708272069557
 log_likelihood: -17213.514538403004
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 45
            nit: 17
     parameters: ParamsDict({'n_hist': 1105.8097431892802, 'n_cont': 67757.4473181167, 't_exp': 103490.52764110279})
         status: 1
        success: True
              x: array([ 7.00833314, 11.12368966, -2.26064874])
```
## Historic exponential change, temporal only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_histexpg_temponly.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temponly.add_size_param("n_hist")
model_inf_histexpg_temponly.add_size_param("n_cont")
model_inf_histexpg_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_temponly.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_temponly.set_size("AHam",g=0, t="t_exp")
#run model
model_inf_histexpg_temponly.optimize(method="TNC")

            fun: 0.017447474709382457
            jac: array([-3.99418586e-07, -1.60079045e-06, -1.50039822e-05])
  kl_divergence: 0.017447474709382457
 log_likelihood: -8965.409706463519
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 17
     parameters: ParamsDict({'n_hist': 1592.379310506907, 'n_cont': 69235.34856627553, 't_exp': 91942.05640388276})
         status: 1
        success: True
              x: array([ 7.3729846 , 11.14526683, -2.40529582])
 ```

## Historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data
model_inf_histexpg_temporal.set_data(sfs,length=467359)
#set parameterss to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temporal.add_size_param("n_hist")
model_inf_histexpg_temporal.add_size_param("n_cont")
model_inf_histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_temporal.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp)
model_inf_histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_histexpg_temporal.add_leaf("CBat",N="n_cont")
model_inf_histexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_histexpg_temporal.optimize(method="TNC")

            fun: 0.6322492835526239
            jac: array([ 1.35578631e-06, -4.80761182e-06, -1.28004297e-05])
  kl_divergence: 0.6322492835526239
 log_likelihood: -30071.714376457756
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 55
            nit: 19
     parameters: ParamsDict({'n_hist': 1532.492452047867, 'n_cont': 57392.801170395054, 't_exp': 104906.50840292552})
         status: 1
        success: True
              x: array([ 7.33465074, 10.95767416, -2.24403589])
```

## Recent and Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_contemp.add_size_param("n_rec")
model_inf_2recexpg_contemp.add_size_param("n_alb")
model_inf_2recexpg_contemp.add_size_param("n_cont")
model_inf_2recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_contemp.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_contemp.set_size("CBat", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec, t= "t_bot")
model_inf_2recexpg_contemp.set_size("CBat", g=0, t="t_rec")
#run model
model_inf_2recexpg_contemp.optimize(method="TNC")

            fun: 0.2341283124129927
            jac: array([-1.83920302e-06, -5.48357279e-06,  2.34146855e-09, -1.24230025e-06,
       -5.18479330e-06])
  kl_divergence: 0.2341283124129927
 log_likelihood: -18257.935484439546
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 97
            nit: 17
     parameters: ParamsDict({'n_rec': 4606.3271955840955, 'n_alb': 10000000000.000004, 'n_cont': 9705700932.846136, 't_rec': 499.84110817558485, 't_bot': 41.32839890274925})
         status: 1
        success: True
              x: array([ 8.43518611, 23.02585093, 22.99597927,  7.80270246, -0.35040592])
```

## Recent and Pre-Albatross exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recexpg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_temporal.add_size_param("n_rec")
model_inf_2recexpg_temporal.add_size_param("n_alb")
model_inf_2recexpg_temporal.add_size_param("n_cont")
model_inf_2recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec)
model_inf_2recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_2recexpg_temporal.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_temporal.set_size("CBat", g=0, t="t_bot")
model_inf_2recexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2recexpg_temporal.optimize(method="TNC")

            fun: 0.7662873084887306
            jac: array([ 2.97013083e-10, -1.97368501e-09,  1.02790733e-09, -1.36186960e-08,
        8.47039950e-08])
  kl_divergence: 0.7662873084887306
 log_likelihood: -31072.7103466806
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 49
            nit: 12
     parameters: ParamsDict({'n_rec': 15119.279427638154, 'n_alb': 30787.532142634274, 'n_cont': 449486024.11043656, 't_rec': 499.9985241585976, 't_bot': 0.01244351808789977})
         status: 1
        success: True
              x: array([ 9.62372599, 10.33486509, 19.92361532, 12.48210256, -8.99160117])
```

## Recent and historic exponential size change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
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
model_inf_2changeexpg_contemp.set_size("CBat",g=0,t="t_exp")
#run model
model_inf_2changeexpg_contemp.optimize(method="TNC")
            
            fun: 0.04870326458602142
            jac: array([-1.08037568e-06, -2.39126804e-06,  1.78473769e-07,  2.02500697e-06,
        1.94357373e-06])
  kl_divergence: 0.04870326458602142
 log_likelihood: -16966.264601276864
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 60
            nit: 17
     parameters: ParamsDict({'n_hist': 8800.361351926438, 'n_alb': 2797631.8179397783, 'n_cont': 79.01453077733078, 't_exp': 26724.918087998245, 't_bot': 80.86368910488369})
         status: 1
        success: True
              x: array([ 9.08254806, 14.84428384,  4.36963177, -4.06376705,  1.44117726])
```

## Recent and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data
model_inf_2changeexpg_temporal.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_2changeexpg_temporal.add_size_param("n_alb")
model_inf_2changeexpg_temporal.add_size_param("n_hist")
model_inf_2changeexpg_temporal.add_size_param("n_cont")
model_inf_2changeexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_hist)/params.t_exp)
model_inf_2changeexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2changeexpg_temporal.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_temporal.set_size("CBat",g=0, t="t_bot")
model_inf_2changeexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2changeexpg_temporal.optimize(method="TNC")

            fun: 0.6236917697992787
            jac: array([2.50421056e-06, 6.43433732e-09, 4.77152186e-07, 4.80104663e-08,
       1.58085691e-06])
  kl_divergence: 0.6236917697992787
 log_likelihood: -30007.806863747774
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 49
            nit: 13
     parameters: ParamsDict({'n_alb': 66818.54173230584, 'n_hist': 1.5242387583720107, 'n_cont': 19.957993552171164, 't_exp': 283607.9583250005, 't_bot': 2.0028977039488485})
         status: 1
        success: True
              x: array([11.10973589,  0.42149511,  2.99362974, -0.96253129, -3.89034293])
```

## Pre-Albatross and historic exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_contemp.add_size_param("n_hist")
model_inf_2histexpg_contemp.add_size_param("n_rec")
model_inf_2histexpg_contemp.add_size_param("n_cont")
model_inf_2histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_contemp.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_contemp.set_size("CBat",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_2histexpg_contemp.set_size("CBat", g=0, t="t_exp")
#run model
model_inf_2histexpg_contemp.optimize(method="TNC")

            fun: 0.04875312527258749
            jac: array([ 2.55170485e-06,  2.41148186e-06,  1.33422366e-06, -3.25887700e-06,
        8.82532455e-06])
  kl_divergence: 0.04875312527258749
 log_likelihood: -16966.611930819483
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 57
            nit: 15
     parameters: ParamsDict({'n_hist': 8170.201048709391, 'n_rec': 3270213.957763446, 'n_cont': 484.29301428905535, 't_exp': 26381.285981877038, 't_rec': 427.43709199799594})
         status: 1
        success: True
              x: array([ 9.0082488 , 15.00036597,  6.18269012, -4.08488018,  1.47267058])
```

## Pre-Albatross and historic exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temponly.add_size_param("n_hist")
model_inf_2histexpg_temponly.add_size_param("n_rec")
model_inf_2histexpg_temponly.add_size_param("n_cont")
model_inf_2histexpg_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temponly.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_2histexpg_temponly.set_size("AHam", g=0, t="t_exp")
#run model
model_inf_2histexpg_temponly.optimize(method="TNC")

            fun: 0.014855546453198538
            jac: array([ 8.63650788e-08,  1.03075734e-07, -3.23052083e-08, -1.82079037e-07,
       -1.58294612e-07])
  kl_divergence: 0.014855546453198538
 log_likelihood: -8953.398710924363
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 44
            nit: 13
     parameters: ParamsDict({'n_hist': 7674.72553618606, 'n_rec': 145347.22033216464, 'n_cont': 1815.4414903846935, 't_exp': 37874.53502194981, 't_rec': 411.28348551966764})
         status: 1
        success: True
              x: array([ 8.94568781, 11.88688078,  7.50408396, -3.54143131,  1.21928092])
```

## Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temporal.add_size_param("n_hist")
model_inf_2histexpg_temporal.add_size_param("n_rec")
model_inf_2histexpg_temporal.add_size_param("n_cont")
model_inf_2histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_temporal.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_2histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2histexpg_temporal.add_leaf("CBat", N="n_cont",t=109)
model_inf_2histexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2histexpg_temporal.optimize(method="TNC")

            fun: 0.5929774794424969
            jac: array([ 5.94083936e-07, -4.12913877e-07, -1.75634577e-06,  2.53313474e-06,
       -1.15814202e-09])
  kl_divergence: 0.5929774794424969
 log_likelihood: -29778.432543363328
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 60
            nit: 17
     parameters: ParamsDict({'n_hist': 471.49520161520456, 'n_rec': 86806.07892093407, 'n_cont': 2073.9298803753272, 't_exp': 106990.35375540002, 't_rec': 499.99997954233174})
         status: 1
        success: True
              x: array([ 6.15590893, 11.37143193,  7.63720058, -2.21998586, 16.76073206])
```

## Recent, Pre-Albatross and historical exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_3expg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3expg_contemp.add_size_param("n_hist")
model_inf_3expg_contemp.add_size_param("n_rec")
model_inf_3expg_contemp.add_size_param("n_alb")
model_inf_3expg_contemp.add_size_param("n_cont")
model_inf_3expg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3expg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3expg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_3expg_contemp.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_contemp.set_size("CBat", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec, t= "t_bot")
model_inf_3expg_contemp.set_size("CBat",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_3expg_contemp.set_size("CBat", g=0, t="t_exp")
#run model
model_inf_3expg_contemp.optimize(method="TNC")

            fun: 0.048700818618524266
            jac: array([ 2.61423254e-07, -2.13248281e-07, -6.05016562e-07,  4.38828835e-06,
        2.24379117e-06,  1.12729465e-07,  1.91367079e-06])
  kl_divergence: 0.048700818618524266
 log_likelihood: -16966.24756266728
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 36
            nit: 9
     parameters: ParamsDict({'n_hist': 16721.642016569283, 'n_rec': 5333933.921869299, 'n_alb': 27461.498991393968, 'n_cont': 57.46557302774521, 't_exp': 26716.255486020556, 't_rec': 263.7220099155648, 't_bot': 34.09573734933886})
         status: 1
        success: True
              x: array([ 9.72445909, 15.4895996 , 10.22054027,  4.05118604, -4.06429403,
       -0.4363897 , -0.65903075])
```

## Recent, Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.393, muts_per_gen=2.5e-8)
#add data to model
model_inf_3expg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3expg_temporal.add_size_param("n_hist")
model_inf_3expg_temporal.add_size_param("n_rec")
model_inf_3expg_temporal.add_size_param("n_alb")
model_inf_3expg_temporal.add_size_param("n_cont")
model_inf_3expg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3expg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3expg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_3expg_temporal.add_leaf("AHam",N="n_alb", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec)
model_inf_3expg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_3expg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_3expg_temporal.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_temporal.set_size("CBat", g=0, t="t_bot")
model_inf_3expg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_3expg_temporal.optimize(method="TNC")

            fun: 0.5962904924878343
            jac: array([ 1.51083588e-06, -3.11087504e-05, -5.11584109e-05, -1.13178972e-05,
        1.15284536e-05, -9.59756616e-09, -3.79868961e-15])
  kl_divergence: 0.5962904924878343
 log_likelihood: -29803.174124785906
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 25
            nit: 5
     parameters: ParamsDict({'n_hist': 3.8194796502711603, 'n_rec': 89138.70920662006, 'n_cont': 10000000000.000004, 'n_alb': 2446.3417067262058, 't_exp': 211682.4178762637, 't_rec': 499.99981425294624, 't_bot': 99.99999999996633})
         status: 1
        success: True
              x: array([ 1.3401142 , 11.39794897, 23.02585093,  7.80234901, -1.36320676,
       14.5547036 , 28.71993117])
```
