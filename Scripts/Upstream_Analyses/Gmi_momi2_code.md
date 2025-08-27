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
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CBat",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.2506224208032932
            jac: array([1.22467999e-12])
  kl_divergence: 0.2506224208032932
 log_likelihood: -18372.83344348638
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_constant': 115465.61966796315})
         status: 1
        success: True
              x: array([11.6567281])
```

## Constant pop size, temporal only

```bash
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_temponly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.1622084701022879
            jac: array([2.16901899e-10])
  kl_divergence: 0.1622084701022879
 log_likelihood: -9636.232159114243
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 10
            nit: 4
     parameters: ParamsDict({'n_constant': 114456.6863871662})
         status: 1
        success: True
              x: array([11.64795175])
```

## Constant pop size, temp & contemp

```bash
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=467359)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CBat",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)#adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CBat","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.7727099118273987
            jac: array([8.48773469e-13])
  kl_divergence: 0.7727099118273987
 log_likelihood: -31120.674348413773
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_constant': 114546.46410280639})
         status: 1
        success: True
              x: array([11.64873582])
```

## Recent size change, contemp only

```bash
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.2496786701477289
            jac: array([ 9.27704030e-09, -1.05856874e-08, -3.26120352e-06])
  kl_divergence: 0.2496786701477289
 log_likelihood: -18366.259276419718
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 30
            nit: 10
     parameters: ParamsDict({'n_alb': 114404.28263921668, 'n_bot': 10000000000.000004, 't_bot': 99.64920514269306})
         status: 1
        success: True
              x: array([11.64749379, 23.02585093,  5.64920975])
```

## Recent size change, temp and contemp

```bash
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.7727004710484383
            jac: array([-4.68489606e-09, -2.88283881e-07, -1.73732681e-06])
  kl_divergence: 0.7727004710484383
 log_likelihood: -31120.603844676498
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 26
            nit: 9
     parameters: ParamsDict({'n_alb': 114512.94212236632, 'n_bot': 591805.699880096, 't_bot': 50.59623795485525})
         status: 1
        success: True
              x: array([11.64844313, 13.29093365,  0.02385065])
```

## Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.24617876588321644
            jac: array([-3.70166483e-06, -4.32131431e-08, -1.15510167e-08])
  kl_divergence: 0.24617876588321644
 log_likelihood: -18341.878943313124
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 7
            nit: 3
     parameters: ParamsDict({'n_rec': 110790.70287119939, 'n_cont': 10000000000.000004, 't_rec': 499.9985814164964})
         status: 1
        success: True
              x: array([11.61539814, 23.02585093, 12.52167213])
```

## Pre-Albatross size change, temporal only

```bash
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temponly.add_size_param("n_rec")
model_inf_recchange_temponly.add_size_param("n_cont")
model_inf_recchange_temponly.add_time_param("t_rec",upper=390)
model_inf_recchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_recchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_recchange_temponly.optimize(method="TNC")

            fun: 0.15881121697536954
            jac: array([ 2.73230875e-09, -3.81136792e-08, -9.13351813e-07])
  kl_divergence: 0.15881121697536954
 log_likelihood: -9620.489288124103
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 12
     parameters: ParamsDict({'n_rec': 111715.92671574, 'n_cont': 10000000000.000004, 't_rec': 389.8942490569808})
         status: 1
        success: True
              x: array([11.62371456, 23.02585093,  8.21254409])
```

## Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.7716331859824356
            jac: array([-1.26238647e-08, -1.54076999e-10, -5.55011237e-10])
  kl_divergence: 0.7716331859824356
 log_likelihood: -31112.63335980359
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 1
            nit: 0
     parameters: ParamsDict({'n_rec': 113739.82851710713, 'n_cont': 10000000000.000004, 't_rec': 499.9996961495581})
         status: 0
        success: True
              x: array([11.64166891, 23.02585093, 14.06255351])
```

## Historic size change, contemp only

```bash
#specify model
model_inf_histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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
			    
            fun: 0.07152124360759715
            jac: array([-1.14056398e-06, -1.07484737e-06, -2.21745557e-07])
  kl_divergence: 0.07152124360759715
 log_likelihood: -17125.21464314116
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 36
            nit: 15
     parameters: ParamsDict({'n_hist': 53056.70516575036, 'n_cont': 421085.3801397299, 't_exp': 218138.44082369638})
         status: 1
        success: True
              x: array([10.87911653, 12.9505909 , -1.32347425])
```

## Historic size change, temporal only

```bash
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temponly.add_size_param("n_hist")
model_inf_histchange_temponly.add_size_param("n_cont")
model_inf_histchange_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_histchange_temponly.optimize(method="TNC")
			    
            fun: 0.024759104503991176
            jac: array([ 1.62035024e-06,  5.36691939e-06, -7.19535575e-05])
  kl_divergence: 0.024759104503991176
 log_likelihood: -8999.291798931736
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 41
            nit: 13
     parameters: ParamsDict({'n_hist': 67779.54943459321, 'n_cont': 345339.5162909038, 't_exp': 99603.89423840042})
         status: 1
        success: True
              x: array([11.1240158 , 12.75228332,  5.74801939])
```

## Historic size change, temp and contemp		  

```bash
#specify model
model_inf_histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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
			    
            fun: 0.6160235824007526
            jac: array([ 3.30775334e-05, -3.78825530e-06,  2.06343907e-05])
  kl_divergence: 0.6160235824007526
 log_likelihood: -29950.54084025558
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 52
            nit: 17
     parameters: ParamsDict({'n_hist': 50371.26481853803, 'n_cont': 381062.94312186056, 't_exp': 231947.83485693182})
         status: 1
        success: True
              x: array([10.82717615, 12.85071985, -1.24141528])
```
			  
## Recent and Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.24617875722994387
            jac: array([ 2.14620578e-08, -3.80383358e-08, -5.17589726e-09, -2.96861553e-09,
        0.00000000e+00])
  kl_divergence: 0.24617875722994387
 log_likelihood: -18341.878883034427
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 1
            nit: 0
     parameters: ParamsDict({'n_rec': 110794.9441262405, 'n_alb': 10000000000.000004, 'n_cont': 10000000000.000004, 't_rec': 499.9996354128366, 't_bot': 62.28220155243827})
         status: 0
        success: True
              x: array([11.61543642, 23.02585093, 23.02585093, 13.88032331,  0.50154361])
```

## Recent and Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.7712386038041957
            jac: array([ 1.09050380e-04, -1.11965328e-08, -1.27941951e-07, -7.03307574e-07,
        7.54765617e-06])
  kl_divergence: 0.7712386038041957
 log_likelihood: -31109.686620096494
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 60
            nit: 12
     parameters: ParamsDict({'n_rec': 113006.59306005535, 'n_alb': 10000000000.000004, 'n_cont': 18475.630920407875, 't_rec': 499.7902665581774, 't_bot': 44.0783407065318})
         status: 1
        success: True
              x: array([11.63520144, 23.02585093,  9.8242079 ,  7.52495792, -0.23798325])
```
			  
## Recent and historic size change, contemp only

```bash
#specify model
model_inf_2change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.04896093024202654
            jac: array([-5.48842062e-08, -6.84947153e-07,  1.05059294e-07,  4.01767850e-06,
        7.91531735e-07])
  kl_divergence: 0.04896093024202654
 log_likelihood: -16968.059500236595
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 92
            nit: 21
     parameters: ParamsDict({'n_hist': 67143.03367478559, 'n_alb': 2096045.5761089486, 'n_cont': 1026.3834430724612, 't_exp': 152106.19760542616, 't_bot': 92.4503795278016})
         status: 1
        success: True
              x: array([11.11458045, 14.55556307,  6.93379668, -1.78618075,  2.50517477])
```

## Recent and historic size change, temp and contemp

```bash
#specify model
model_inf_2change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2change_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_temporal.add_size_param("n_hist")
model_inf_2change_temporal.add_size_param("n_alb")
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

            fun: 0.6091397426634759
            jac: array([ 7.24556871e-07, -7.10798522e-07, -4.83258251e-08, -2.67394472e-09,
        2.60058618e-06])
  kl_divergence: 0.6091397426634759
 log_likelihood: -29899.1323250976
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_hist': 53139.39413036798, 'n_alb': 410043.5438496062, 'n_cont': 522.5157166814631, 't_exp': 220098.69802963772, 't_bot': 5.741563505409801})
         status: 1
        success: True
              x: array([10.88067382, 12.92401864,  6.25865506, -1.31158996, -2.79830877])
```

## Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.04898645101097275
            jac: array([ 1.37334521e-06,  6.55672430e-07, -1.58840831e-06, -8.64885108e-06,
        1.12646777e-05])
  kl_divergence: 0.04898645101097275
 log_likelihood: -16968.237277913075
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 73
            nit: 17
     parameters: ParamsDict({'n_hist': 67243.23553907468, 'n_rec': 2282421.280210224, 'n_cont': 3428.3783929658534, 't_exp': 151602.77455614725, 't_rec': 318.01728712903326})
         status: 1
        success: True
              x: array([11.11607171, 14.6407474 ,  8.13984266, -1.79032318,  0.1288906 ])
```

## Pre-Albatross and historic size change, temporal only

```bash
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histchange_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_temponly.add_size_param("n_hist")
model_inf_2histchange_temponly.add_size_param("n_rec")
model_inf_2histchange_temponly.add_size_param("n_cont")
model_inf_2histchange_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_2histchange_temponly.add_time_param("t_rec",upper=390)
model_inf_2histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_2histchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
model_inf_2histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_temponly.optimize(method="TNC")

            fun: 0.019144311746246142
            jac: array([ 2.57345575e-05, -4.34254433e-06,  1.32807685e-05, -4.71523646e-07,
       -1.30776474e-05])
  kl_divergence: 0.019144311746246142
 log_likelihood: -8973.272849292345
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 25
            nit: 8
     parameters: ParamsDict({'n_hist': 71208.77533550606, 'n_rec': 123885246.02000551, 'n_cont': 85.99741865592048, 't_exp': 99886.92099995934, 't_rec': 9.948428522273133})
         status: 1
        success: True
              x: array([11.17337134, 18.63486626,  4.45431728, 10.28292586, -3.64289236])
```

## Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.5859013553276375
            jac: array([-1.35539123e-06,  8.13476987e-07, -1.77171398e-06,  1.92939629e-06,
       -4.73202528e-07])
  kl_divergence: 0.5859013553276375
 log_likelihood: -29725.588048473557
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 61
            nit: 15
     parameters: ParamsDict({'n_hist': 62673.54523930567, 'n_rec': 669433.5271880438, 'n_cont': 10701.89618315035, 't_exp': 176233.13603213706, 't_rec': 499.97034728817414})
         status: 1
        success: True
              x: array([11.04569471, 13.41418715,  9.27817622, -1.60049632,  9.48170481])
```

## Recent, Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.0489505510098634
            jac: array([ 1.60989456e-05, -9.27940545e-06, -3.20745845e-06, -3.28184212e-05,
       -3.83419455e-06,  1.13477384e-06,  3.26619287e-05])
  kl_divergence: 0.0489505510098634
 log_likelihood: -16967.987198505347
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 54
            nit: 10
     parameters: ParamsDict({'n_hist': 67142.00408071902, 'n_rec': 2081679.86124527, 'n_alb': 26303.182165432525, 'n_cont': 8.608880417890582, 't_exp': 152083.28459306044, 't_rec': 201.78744686009026, 't_bot': 0.7078033429324033})
         status: 1
        success: True
              x: array([11.11456512, 14.54868575, 10.17744521,  2.15279428, -1.78636902,
       -1.18928547, -4.94365597])
```

## Recent, Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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
8
            fun: 0.5797897015755427
            jac: array([-7.16400414e-09, -9.92673337e-08, -1.20640497e-07, -1.79714280e-09,
       -6.42354786e-06, -1.51827282e-07, -5.47511552e-10])
  kl_divergence: 0.5797897015755427
 log_likelihood: -29679.946218252913
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 3
     parameters: ParamsDict({'n_hist': 67528.11024947988, 'n_rec': 1325909.870736475, 'n_alb': 5082.579646472889, 'n_cont': 10000000000.000004, 't_exp': 153112.14915551632, 't_rec': 499.2826752239294, 't_bot': 99.99998471125822})
         status: 1
        success: True
              x: array([11.12029924, 14.09760948,  8.53357422, 23.02585093, -1.7779397 ,
        6.29396019, 15.69356387])
```

## Recent exponential change, contemp only 

```bash
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
#specify model
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.24975779833261616
            jac: array([ 1.59422249e-09, -7.21632946e-06, -3.45982389e-07])
  kl_divergence: 0.24975779833261616
 log_likelihood: -18366.810483355643
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 34
            nit: 12
     parameters: ParamsDict({'n_alb': 114484.00206821255, 'n_bot': 10000000000.000004, 't_bot': 99.95943650913726})
         status: 1
        success: True
              x: array([11.64819037, 23.02585093,  7.80965133])
```

## Recent exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data
model_inf_expg_temporal.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_expg_temporal.add_size_param("n_alb")
model_inf_expg_temporal.add_size_param("n_bot")
model_inf_expg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_expg_temporal.add_leaf("CBat",N="n_bot",g=lambda params: log(params.n_bot/params.n_alb)/params.t_bot)
model_inf_expg_temporal.set_size("CBat",g=0, t="t_bot")
model_inf_expg_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_expg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_expg_temporal.optimize(method="TNC")

            fun: 0.7727006312549426
            jac: array([1.80468594e-10, 8.32080628e-09, 1.70390732e-07])
  kl_divergence: 0.7727006312549426
 log_likelihood: -31120.60504109867
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 10
     parameters: ParamsDict({'n_alb': 114511.2559611455, 'n_bot': 9999999123.055517, 't_bot': 55.12053743450615})
         status: 1
        success: True
              x: array([11.6484284 , 23.02585084,  0.20554209])
```
		  
## Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.24653743697673552
            jac: array([ 9.76008271e-10, -3.15852267e-05, -1.71105675e-07])
  kl_divergence: 0.24653743697673552
 log_likelihood: -18344.377446150578
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 30
            nit: 10
     parameters: ParamsDict({'n_rec': 111098.14277418869, 'n_cont': 10000000000.000004, 't_rec': 499.9773610400062})
         status: 1
        success: True
              x: array([11.61816926, 23.02585093,  9.75160411])
```

## Pre-Albatross exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recexpg_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temponly.add_size_param("n_rec")
model_inf_recexpg_temponly.add_size_param("n_cont")
model_inf_recexpg_temponly.add_time_param("t_rec",upper=390)
model_inf_recexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temponly.set_size("AHam", g=0, t="t_rec")
#run model
model_inf_recexpg_temponly.optimize(method="TNC")


            fun: 0.15915863230543387
            jac: array([-4.27399713e-05, -2.55970061e-05, -5.13879751e-05])
  kl_divergence: 0.15915863230543387
 log_likelihood: -9622.099210763621
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 35
            nit: 14
     parameters: ParamsDict({'n_rec': 111940.44445830795, 'n_cont': 10000000000.000004, 't_rec': 383.3811952569464})
         status: 1
        success: True
              x: array([11.62572226, 23.02585093,  4.05911498])
```

## Pre-Albatross size change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.7719207813483154
            jac: array([ 7.22277326e-07, -6.69208619e-06, -6.05272799e-05])
  kl_divergence: 0.7719207813483154
 log_likelihood: -31114.78112199598
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 29
            nit: 9
     parameters: ParamsDict({'n_rec': 113950.11387403413, 'n_cont': 10000000000.000004, 't_rec': 465.1054682033294})
         status: 1
        success: True
              x: array([11.64351603, 23.02585093,  2.31726467])
```
			  
## Historic exponential change, contemp only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.0841970087672007
            jac: array([ 2.58292167e-06, -1.30941533e-05,  2.35141728e-06])
  kl_divergence: 0.0841970087672007
 log_likelihood: -17213.51402324296
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 23
            nit: 9
     parameters: ParamsDict({'n_hist': 6066.383813958035, 'n_cont': 483919.7640729334, 't_exp': 787288.8017318724})
         status: 1
        success: True
              x: array([ 8.71051796, 13.0896744 ,  1.2958766 ])
```
## Historic exponential change, temporal only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histexpg_temponly.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temponly.add_size_param("n_hist")
model_inf_histexpg_temponly.add_size_param("n_cont")
model_inf_histexpg_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_histexpg_temponly.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_temponly.set_size("AHam",g=0, t="t_exp")
#run model
model_inf_histexpg_temponly.optimize(method="TNC")

            fun: 0.03905122609425718
            jac: array([-1.57997589e-05, -1.26738145e-05, -5.07905176e-06])
  kl_divergence: 0.03905122609425718
 log_likelihood: -9065.521490381028
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 77
            nit: 23
     parameters: ParamsDict({'n_hist': 72364.67312977242, 'n_cont': 556865.1433996052, 't_exp': 99872.5461735106})
         status: 1
        success: True
              x: array([11.18947352, 13.23007838,  8.54781209])
 ```

## Historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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

            fun: 0.6318302819570523
            jac: array([ 1.87374789e-06, -1.85154421e-06, -4.71430686e-06])
  kl_divergence: 0.6318302819570523
 log_likelihood: -30068.585272542026
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 41
            nit: 12
     parameters: ParamsDict({'n_hist': 11086.773084443048, 'n_cont': 410101.0256964739, 't_exp': 736304.334230244})
         status: 1
        success: True
              x: array([ 9.31350806, 12.92415881,  1.01317346])
```

## Recent and Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_contemp.add_size_param("n_rec")
model_inf_2recexpg_contemp.add_size_param("n_alb")
model_inf_2recexpg_contemp.add_size_param("n_cont")
model_inf_2recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_contemp.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_contemp.set_size("CBat", g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot), t= "t_bot")
model_inf_2recexpg_contemp.set_size("CBat", g=0, t="t_rec")
#run model
model_inf_2recexpg_contemp.optimize(method="TNC")

            fun: 0.24646527975199975
            jac: array([-5.67854581e-09, -2.52295969e-05,  5.15476213e-09, -4.17435562e-09,
       -7.84696124e-10])
  kl_divergence: 0.24646527975199975
 log_likelihood: -18343.87479892307
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 51
            nit: 13
     parameters: ParamsDict({'n_rec': 111036.85472847187, 'n_alb': 10000000000.000004, 'n_cont': 9999948126.38765, 't_rec': 499.9994458122683, 't_bot': 99.99890887594451})
         status: 1
        success: True
              x: array([11.61761745, 23.02585093, 23.02584574, 13.46158498, 11.42570615])
```

## Recent and Pre-Albatross exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recexpg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_temporal.add_size_param("n_rec")
model_inf_2recexpg_temporal.add_size_param("n_alb")
model_inf_2recexpg_temporal.add_size_param("n_cont")
model_inf_2recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot))
model_inf_2recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_2recexpg_temporal.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_temporal.set_size("CBat", g=0, t="t_bot")
model_inf_2recexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2recexpg_temporal.optimize(method="TNC")

            fun: 0.7700743446475339
            jac: array([ 5.38945052e-06,  5.41867481e-06, -4.15498593e-05,  5.09533058e-07,
       -5.82121299e-10])
  kl_divergence: 0.7700743446475339
 log_likelihood: -31100.991932714544
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 100
            nit: 26
     parameters: ParamsDict({'n_rec': 19277.733276057173, 'n_alb': 6030.961532891244, 'n_cont': 10000000000.000004, 't_rec': 165.3119103444535, 't_bot': 99.99999266956658})
         status: 1
        success: True
              x: array([ 9.86670599,  8.70466174, 23.02585093, -1.81845548, 16.42864603])
```

## Recent and historic exponential size change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2changeexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2changeexpg_contemp.add_size_param("n_hist")
model_inf_2changeexpg_contemp.add_size_param("n_alb")
model_inf_2changeexpg_contemp.add_size_param("n_cont")
model_inf_2changeexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_contemp.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_contemp.set_size("CBat",g=lambda params: log(params.n_alb/params.n_hist)/(params.t_exp-params.t_bot), t= "t_bot")
model_inf_2changeexpg_contemp.set_size("CBat",g=0,t="t_exp")
#run model
model_inf_2changeexpg_contemp.optimize(method="TNC")
            
            fun: 0.048694876398121374
            jac: array([-2.32567451e-07, -3.06028332e-07, -1.33168216e-07,  8.98405151e-06,
        3.10093545e-07])
  kl_divergence: 0.048694876398121374
 log_likelihood: -16966.206169159952
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 71
            nit: 18
     parameters: ParamsDict({'n_hist': 63954.62473358225, 'n_alb': 20046711.84161918, 'n_cont': 60.349407868520636, 't_exp': 191147.55090317182, 't_bot': 74.85833552215723})
         status: 1
        success: True
              x: array([11.06592912, 16.8135757 ,  4.10015114, -1.49630462,  1.09107106])
```

## Recent and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data
model_inf_2changeexpg_temporal.set_data(sfs,length=467359)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_2changeexpg_temporal.add_size_param("n_alb")
model_inf_2changeexpg_temporal.add_size_param("n_hist")
model_inf_2changeexpg_temporal.add_size_param("n_cont")
model_inf_2changeexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_hist)/(params.t_exp-params.t_bot))
model_inf_2changeexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2changeexpg_temporal.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_temporal.set_size("CBat",g=0, t="t_bot")
model_inf_2changeexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2changeexpg_temporal.optimize(method="TNC")

            fun: 0.6217755789577732
            jac: array([-2.00821275e-05,  9.81952519e-08, -1.93176344e-06, -1.60339408e-07,
        3.94780016e-06])
  kl_divergence: 0.6217755789577732
 log_likelihood: -29993.49675054341
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 69
            nit: 16
     parameters: ParamsDict({'n_alb': 470372.4943323863, 'n_hist': 7257.289618249424, 'n_cont': 147.87316061668506, 't_exp': 771081.3531621123, 't_bot': 15.735037650861958})
         status: 1
        success: True
              x: array([13.0612802 ,  8.88976171,  4.99635488,  1.20137357, -1.67807622])
```

## Pre-Albatross and historic exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_contemp.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_contemp.add_size_param("n_hist")
model_inf_2histexpg_contemp.add_size_param("n_rec")
model_inf_2histexpg_contemp.add_size_param("n_cont")
model_inf_2histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_contemp.add_leaf("CBat",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_contemp.set_size("CBat",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_2histexpg_contemp.set_size("CBat", g=0, t="t_exp")
#run model
model_inf_2histexpg_contemp.optimize(method="TNC")

            fun: 0.04870143866575696
            jac: array([-9.29925129e-07, -2.36697185e-07, -6.18213239e-07,  1.46720477e-05,
        4.29786659e-07])
  kl_divergence: 0.04870143866575696
 log_likelihood: -16966.2518819163
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 65
            nit: 17
     parameters: ParamsDict({'n_hist': 63979.41977182353, 'n_rec': 20155433.525938973, 'n_cont': 454.7736253818734, 't_exp': 191004.6615651123, 't_rec': 475.85504341004344})
         status: 1
        success: True
              x: array([11.06631674, 16.81898446,  6.11979977, -1.49727037,  2.71542461])
```

## Pre-Albatross and historic exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_temponly.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temponly.add_size_param("n_hist")
model_inf_2histexpg_temponly.add_size_param("n_rec")
model_inf_2histexpg_temponly.add_size_param("n_cont")
model_inf_2histexpg_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_2histexpg_temponly.add_time_param("t_rec",upper=390)
model_inf_2histexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temponly.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_2histexpg_temponly.set_size("AHam", g=0, t="t_exp")
#run model
model_inf_2histexpg_temponly.optimize(method="TNC")

            fun: 0.020592961298606535
            jac: array([ 5.72636343e-08, -1.38335032e-04,  1.99268276e-08, -3.72074061e-10,
       -9.59399623e-07])
  kl_divergence: 0.020592961298606535
 log_likelihood: -8979.985891317983
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 81
            nit: 18
     parameters: ParamsDict({'n_hist': 70695.41820799882, 'n_rec': 10000000000.000004, 'n_cont': 34.51280108351472, 't_exp': 99889.99805861108, 't_rec': 71.57199534455103})
         status: 1
        success: True
              x: array([11.16613604, 23.02585093,  3.5413303 , 17.65191655, -1.49269253])
```

## Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temporal.add_size_param("n_hist")
model_inf_2histexpg_temporal.add_size_param("n_rec")
model_inf_2histexpg_temporal.add_size_param("n_cont")
model_inf_2histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_temporal.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_2histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2histexpg_temporal.add_leaf("CBat", N="n_cont",t=109)
model_inf_2histexpg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_2histexpg_temporal.optimize(method="TNC")

            fun: 0.6003664812378268
            jac: array([-9.17736451e-07,  1.96622921e-06,  1.22769871e-06, -1.14166640e-06,
       -2.36938484e-08])
  kl_divergence: 0.6003664812378268
 log_likelihood: -29833.61360877085
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 84
            nit: 22
     parameters: ParamsDict({'n_hist': 5223.091373114251, 'n_rec': 516131.06659131875, 'n_cont': 1841.7575233247223, 't_exp': 735165.6937540155, 't_rec': 499.9996193283845})
         status: 1
        success: True
              x: array([ 8.56084472, 13.15411602,  7.51847557,  1.0072958 , 13.83715182])
```

## Recent, Pre-Albatross and historical exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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
model_inf_3expg_contemp.set_size("CBat", g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot), t= "t_bot")
model_inf_3expg_contemp.set_size("CBat",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_3expg_contemp.set_size("CBat", g=0, t="t_exp")
#run model
model_inf_3expg_contemp.optimize(method="TNC")

            fun: 0.04869988177562233
            jac: array([ 1.83978431e-05, -5.71947792e-07, -7.03227126e-06,  1.45843692e-07,
       -8.86449549e-06,  2.68862964e-06, -9.50281716e-08])
  kl_divergence: 0.04869988177562233
 log_likelihood: -16966.241036619624
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 94
            nit: 18
     parameters: ParamsDict({'n_hist': 64040.493932178906, 'n_rec': 20745065.087824628, 'n_alb': 377.93222842305914, 'n_cont': 71.69352379985685, 't_exp': 190558.39691226542, 't_rec': 398.08121661903743, 't_bot': 0.22290089826751072})
         status: 1
        success: True
              x: array([11.06727088, 16.84781895,  5.93471489,  4.27240042, -1.50029038,
        1.03558891, -6.1039667 ])
```

## Recent, Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
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
model_inf_3expg_temporal.add_leaf("AHam",N="n_alb", g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot))
model_inf_3expg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_3expg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_3expg_temporal.add_leaf("CBat", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_temporal.set_size("CBat", g=0, t="t_bot")
model_inf_3expg_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_3expg_temporal.optimize(method="TNC")

            fun: 0.600541243654024
            jac: array([-1.37095976e-08, -8.02489930e-10, -1.51420198e-09, -3.85820453e-05,
       -7.63768345e-10, -7.53716537e-07, -1.26733213e-07])
  kl_divergence: 0.600541243654024
 log_likelihood: -29834.91873449501
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 68
            nit: 16
     parameters: ParamsDict({'n_hist': 342.01145256119304, 'n_rec': 180441.35968887853, 'n_alb': 2223.337517525294, 'n_cont': 10000000000.000004, 't_exp': 991116.0484110452, 't_rec': 499.98526480586435, 't_bot': 99.99850854022061})
         status: 1
        success: True
              x: array([ 5.83484422, 12.10316113,  7.70676473, 23.02585093,  4.70444429,
       10.18105795, 11.11315519])
```
