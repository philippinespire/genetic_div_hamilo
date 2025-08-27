# Running MOMI2 on Wahab

Code for *Equulites laterofenestra*.

## Create input files for momi2

```sh
cd ~/PIRE/leiognathus_leuciscus/momi2
salloc
module load container_env python3
bash

#make allele counts file to make SFS from
crun.python3 -p ~/.conda/envs/momi-py36 python -m momi.read_vcf --no_aa ../lle.all.ela.nohighhet.vcf.gz ind2pop.txt lle.all.ela.nohighhet.allelecounts

#make SFS
crun.python3 -p ~/.conda/envs/momi-py36 python -m momi.extract_sfs lle.all.ela.nohighhet.sfs.gz 100 lle.all.ela.nohighhet.allelecounts

#read SFS into momi and run models
crun.python3 -p ~/.conda/envs/momi-py36 python
import momi
from autograd.numpy import log

sfsfile="lle.all.ela.nohighhet.sfs.gz"
sfs = momi.Sfs.load(sfsfile)
NeConstant=1e4

#check file
print("populations", sfs.populations)
populations ('AHam', 'CNas')
print("percent missing data per population", sfs.p_missing)
percent missing data per population [0.04156797 0.00903024]
```

## Constant pop size, contemp only

```bash
#specify model
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CNas",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.15475256184257447
            jac: array([1.84050544e-09])
  kl_divergence: 0.15475256184257447
 log_likelihood: -164614.5161962893
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_constant': 862440.5056874502})
         status: 1
        success: True
              x: array([13.66752145])
```

## Constant pop size, temporal only

```bash
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_temponly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.06514603518730128
            jac: array([5.58987528e-10])
  kl_divergence: 0.06514603518730128
 log_likelihood: -49024.40797016061
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 4
     parameters: ParamsDict({'n_constant': 914717.8014384819})
         status: 1
        success: True
              x: array([13.72637088])
```

## Constant pop size, temp & contemp

```bash
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=485532)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CNas",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109) #adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CNas","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.3798758029717348
            jac: array([2.38630056e-08])
  kl_divergence: 0.3798758029717348
 log_likelihood: -232942.24344177224
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 13
            nit: 3
     parameters: ParamsDict({'n_constant': 866468.6894291132})
         status: 0
        success: True
              x: array([13.67218125])
```

## Recent size change, contemp only

```bash
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_change_contemp.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_change_contemp.add_size_param("n_alb")
model_inf_change_contemp.add_size_param("n_bot")
model_inf_change_contemp.add_time_param("t_bot",upper=1e2) #force bot to be in last 100 years (gens?)
model_inf_change_contemp.add_leaf("CNas",N="n_bot")
model_inf_change_contemp.set_size("CNas", N="n_alb", t="t_bot") #says CBat pop changes from n_alb to n_bot at t_bot
#run model
model_inf_change_contemp.optimize(method="TNC")

            fun: 0.13793470537879832
            jac: array([-6.85756250e-07, -1.54665496e-06,  1.86091039e-05])
  kl_divergence: 0.13793470537879832
 log_likelihood: -163853.8950020021
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 16
            nit: 6
     parameters: ParamsDict({'n_alb': 867423.0264118328, 'n_bot': 4390.475653793292, 't_bot': 65.86648475737071})
         status: 1
        success: True
              x: array([13.67328206,  8.38719285,  0.65734998])
```

## Recent size change, temp and contemp

```bash
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_change_temporal.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_change_temporal.add_size_param("n_alb")
model_inf_change_temporal.add_size_param("n_bot")
model_inf_change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_change_temporal.add_leaf("CNas",N="n_bot")
model_inf_change_temporal.set_size("CNas", N="n_alb", t="t_bot")
model_inf_change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_change_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_change_temporal.optimize(method="TNC")

            fun: 0.33953234699666857
            jac: array([1.30063625e-06, 3.71170521e-08, 2.28405293e-05])
  kl_divergence: 0.33953234699666857
 log_likelihood: -231053.36283301964
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 20
            nit: 5
     parameters: ParamsDict({'n_alb': 874134.5803649353, 'n_bot': 2594.5955820918985, 't_bot': 58.93605618250327})
         status: 1
        success: True
              x: array([13.68098962,  7.86118594,  0.36132261])
```

## Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recchange_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_contemp.add_size_param("n_rec")
model_inf_recchange_contemp.add_size_param("n_cont")
model_inf_recchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_contemp.add_leaf("CNas",N="n_cont")
model_inf_recchange_contemp.set_size("CNas", N="n_rec", t="t_rec")
#run model
model_inf_recchange_contemp.optimize(method="TNC")

            fun: 0.13798188887641527
            jac: array([-1.20001387e-05, -3.01403838e-06,  1.08463745e-05])
  kl_divergence: 0.13798188887641527
 log_likelihood: -163856.02897004882
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 93
            nit: 19
     parameters: ParamsDict({'n_rec': 867663.9576673126, 'n_cont': 8267.823530152658, 't_rec': 124.60187087535796})
         status: 1
        success: True
              x: array([13.67355977,  9.02012658, -3.31777979])
```
		  
## Pre-Albatross size change, temporal only

```bash
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recchange_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temponly.add_size_param("n_rec")
model_inf_recchange_temponly.add_size_param("n_cont")
model_inf_recchange_temponly.add_time_param("t_rec",upper=390)
model_inf_recchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_recchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_recchange_temponly.optimize(method="TNC")

            fun: 0.0650232338623682
            jac: array([-3.46715318e-07, -5.04355891e-07, -4.33015265e-05])
  kl_divergence: 0.0650232338623682
 log_likelihood: -49021.11063178483
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 29
            nit: 8
     parameters: ParamsDict({'n_rec': 914118.4122795431, 'n_cont': 223473983.03856885, 't_rec': 252.37295579371514})
         status: 1
        success: True
              x: array([13.7257154 , 19.22480556,  0.60636053])
```

## Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temporal.add_size_param("n_rec")
model_inf_recchange_temporal.add_size_param("n_cont")
model_inf_recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_recchange_temporal.add_leaf("AHam", N="n_cont", t=109)
model_inf_recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_recchange_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_recchange_temporal.optimize(method="TNC")

            fun: 0.3395737705050992
            jac: array([-7.54404939e-06, -2.35561465e-06,  4.71440676e-07])
  kl_divergence: 0.3395737705050992
 log_likelihood: -231055.30228168436
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 38
            nit: 9
     parameters: ParamsDict({'n_rec': 874444.592717893, 'n_cont': 4863.395352587883, 't_rec': 111.07086986381833})
         status: 1
        success: True
              x: array([13.68134421,  8.48949211, -8.61030713])
```

## Historic size change, contemp only

```bash
#specify model
model_inf_histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histchange_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_contemp.add_size_param("n_hist")
model_inf_histchange_contemp.add_size_param("n_cont")
model_inf_histchange_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_contemp.add_leaf("CNas",N="n_cont")
model_inf_histchange_contemp.set_size("CNas", N="n_hist", t="t_exp")
#run model
model_inf_histchange_contemp.optimize(method="TNC")
			    
            fun: 0.1444152927716583
            jac: array([-5.08357383e-06,  1.48007815e-06,  3.77751139e-06])
  kl_divergence: 0.1444152927716583
 log_likelihood: -164146.99252801898
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 57
            nit: 19
     parameters: ParamsDict({'n_hist': 902076.6531832581, 'n_cont': 391538.4835135454, 't_exp': 10007.251529369676})
         status: 1
        success: True
              x: array([ 13.71245478,  12.87783909, -11.8242405 ])
```
			  
## Historic size change, temporal only

```bash
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histchange_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temponly.add_size_param("n_hist")
model_inf_histchange_temponly.add_size_param("n_cont")
model_inf_histchange_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_histchange_temponly.optimize(method="TNC")
			    
            fun: 0.03225338059763309
            jac: array([ 7.60862192e-10, -1.42396786e-06, -1.50153344e-07])
  kl_divergence: 0.03225338059763309
 log_likelihood: -48141.20730177343
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 45
            nit: 15
     parameters: ParamsDict({'n_hist': 804748.9061613939, 'n_cont': 10000000000.000004, 't_exp': 99889.16618168529})
         status: 1
        success: True
              x: array([13.59828559, 23.02585093, 11.58929543])
```

## Historic size change, temp and contemp		  

```bash
#specify model
model_inf_histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temporal.add_size_param("n_hist")
model_inf_histchange_temporal.add_size_param("n_cont")
model_inf_histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_histchange_temporal.add_leaf("AHam",N="n_cont",t=109)
model_inf_histchange_temporal.set_size("AHam", N="n_hist", t="t_exp") #says, at some time (t_exp) in past, Alb was at n_hist size
model_inf_histchange_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_histchange_temporal.optimize(method="TNC")
			    
            fun: 0.3463323207427854
            jac: array([-2.79508126e-06,  2.65692526e-05, -1.62583730e-05])
  kl_divergence: 0.3463323207427854
 log_likelihood: -231371.73760381283
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 15
     parameters: ParamsDict({'n_hist': 636388.1775659714, 'n_cont': 1397631.4985445628, 't_exp': 998068.8672182811})
         status: 1
        success: True
              x: array([13.363564  , 14.15028958,  6.23764563])
```

## Recent and Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recchange_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recchange_contemp.add_size_param("n_rec")
model_inf_2recchange_contemp.add_size_param("n_alb")
model_inf_2recchange_contemp.add_size_param("n_cont")
model_inf_2recchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recchange_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recchange_contemp.add_leaf("CNas",N="n_cont")
model_inf_2recchange_contemp.set_size("CNas", N="n_rec", t="t_rec")
model_inf_2recchange_contemp.set_size("CNas", N="n_alb", t="t_bot")
#run model
model_inf_2recchange_contemp.optimize(method="TNC")

            fun: 0.13788119994797515
            jac: array([-4.41502911e-12, -2.16252087e-08,  5.13635130e-07, -9.37255317e-12,
       -5.01395773e-07])
  kl_divergence: 0.13788119994797515
 log_likelihood: -163851.47511188226
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 90
            nit: 19
     parameters: ParamsDict({'n_rec': 867063.3323495127, 'n_alb': 19142998.62555278, 'n_cont': 1.0, 't_rec': 111.0011360759977, 't_bot': 0.015046582347038942})
         status: 1
        success: True
              x: array([ 13.6728673 ,  16.7674476 ,   0.        , -12.74375148,
        -8.80162411])
```

## Recent and Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recchange_temporal.add_size_param("n_rec")
model_inf_2recchange_temporal.add_size_param("n_alb")
model_inf_2recchange_temporal.add_size_param("n_cont")
model_inf_2recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recchange_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_2recchange_temporal.set_size("CNas", N="n_alb", t="t_bot")
model_inf_2recchange_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_2recchange_temporal.move_lineages("CNas","AHam",t=110)
model_inf_2recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_2recchange_temporal.optimize(method="TNC")

            fun: 0.33952125400651906
            jac: array([-1.74226324e-05, -3.56170318e-10,  3.63078063e-07, -8.91346296e-08,
        2.34309245e-05])
  kl_divergence: 0.33952125400651906
 log_likelihood: -231052.84345922084
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 7
            nit: 2
     parameters: ParamsDict({'n_rec': 873922.2483635749, 'n_alb': 671283929.7387217, 'n_cont': 2084.859624470853, 't_rec': 115.33830335459164, 't_bot': 47.483419650785585})
         status: 1
        success: True
              x: array([13.68074669, 20.32470275,  7.64245681, -4.4848809 , -0.10074835])
```

## Recent and historic size change, contemp only

```bash
#specify model
model_inf_2change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2change_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_contemp.add_size_param("n_hist")
model_inf_2change_contemp.add_size_param("n_alb")
model_inf_2change_contemp.add_size_param("n_cont")
model_inf_2change_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2change_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2change_contemp.add_leaf("CNas",N="n_cont")
model_inf_2change_contemp.set_size("CNas", N="n_alb", t="t_bot")
model_inf_2change_contemp.set_size("CNas", N="n_hist", t="t_exp")
#run model
model_inf_2change_contemp.optimize(method="TNC")

            fun: 0.05269709756830118
            jac: array([ 3.00827749e-05, -1.29066105e-05, -9.61149275e-05,  1.55330390e-05,
        9.52705501e-05])
  kl_divergence: 0.05269709756830118
 log_likelihood: -159998.85371355675
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 73
            nit: 22
     parameters: ParamsDict({'n_hist': 690706.4478946883, 'n_alb': 24443399.520275723, 'n_cont': 8.982156900748475, 't_exp': 603471.9700564917, 't_bot': 1.0596441254339273})
         status: 1
        success: True
              x: array([13.44547019, 17.01187078,  2.19524004,  0.40324325, -4.53658408])
```

## Recent and historic size change, temp and contemp

```bash
#specify model
model_inf_2change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2change_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_temporal.add_size_param("n_hist")
model_inf_2change_temporal.add_size_param("n_alb")
model_inf_2change_temporal.add_size_param("n_cont")
model_inf_2change_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2change_temporal.add_leaf("CNas",N="n_cont")
model_inf_2change_temporal.set_size("CNas", N="n_alb", t="t_bot")
model_inf_2change_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_2change_temporal.move_lineages("CNas","AHam",t=110)
model_inf_2change_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2change_temporal.optimize(method="TNC")

            fun: 0.2693536470914207
            jac: array([-1.67952164e-06, -2.52073122e-07, -1.48836563e-06, -1.46251150e-07,
        1.56900029e-05])
  kl_divergence: 0.2693536470914207
 log_likelihood: -227767.59610345593
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 69
            nit: 18
     parameters: ParamsDict({'n_hist': 617548.8746870809, 'n_alb': 2034699.4121091221, 'n_cont': 957.947834266565, 't_exp': 999967.0579046658, 't_bot': 34.316344083098386})
         status: 1
        success: True
              x: array([13.33351349, 14.52585866,  6.86479332, 10.31067561, -0.64922838])
```

## Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histchange_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_contemp.add_size_param("n_hist")
model_inf_2histchange_contemp.add_size_param("n_rec")
model_inf_2histchange_contemp.add_size_param("n_cont")
model_inf_2histchange_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_contemp.add_leaf("CNas",N="n_cont")
model_inf_2histchange_contemp.set_size("CNas", N="n_rec", t="t_rec")
model_inf_2histchange_contemp.set_size("CNas", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_contemp.optimize(method="TNC")

            fun: 0.05272937501586763
            jac: array([ 7.32013215e-08, -3.52776623e-06, -7.49094223e-05,  1.09639690e-05,
        3.88947990e-05])
  kl_divergence: 0.05272937501586763
 log_likelihood: -160000.31352567783
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 22
            nit: 6
     parameters: ParamsDict({'n_hist': 691132.2528781184, 'n_rec': 75778487.14696927, 'n_cont': 1744.8335226242136, 't_exp': 596022.3332892809, 't_rec': 219.57162393814193})
         status: 1
        success: True
              x: array([13.44608648, 18.143325  ,  7.46441443,  0.3719983 , -0.94890826])
```

## Pre-Albatross and historic size change, temporal only

```bash
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histchange_temponly.set_data(sfs, length=485532)
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

            fun: 0.03225325426234555
            jac: array([ 3.06377117e-08, -1.41947854e-06, -4.47808159e-09, -2.38185766e-08,
       -1.79782906e-21])
  kl_divergence: 0.03225325426234555
 log_likelihood: -48141.20390954462
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 36
            nit: 10
     parameters: ParamsDict({'n_hist': 804748.9837076663, 'n_rec': 10000000000.000004, 'n_cont': 10000000000.000004, 't_exp': 99889.86773251485, 't_rec': 319.4488321831875})
         status: 0
        success: True
              x: array([13.59828569, 23.02585093, 23.02585093, 13.43049248,  1.51025888])
```

## Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_temporal.add_size_param("n_hist")
model_inf_2histchange_temporal.add_size_param("n_rec")
model_inf_2histchange_temporal.add_size_param("n_cont")
model_inf_2histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_2histchange_temporal.add_leaf("AHam",N="n_cont",t=109)
model_inf_2histchange_temporal.move_lineages("CNas","AHam",t=110)
model_inf_2histchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_2histchange_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_temporal.optimize(method="TNC")

            fun: 0.32858232040370355
            jac: array([-3.33622494e-07, -1.05364277e-06,  5.49082611e-06, -5.07947669e-06,
        4.10106063e-05])
  kl_divergence: 0.32858232040370355
 log_likelihood: -230540.68258793702
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 64
            nit: 18
     parameters: ParamsDict({'n_hist': 1.2015813028828977, 'n_rec': 655632.0201630667, 'n_cont': 4463.671626557806, 't_exp': 999940.7457080488, 't_rec': 201.79175498447682})
         status: 1
        success: True
              x: array([ 0.18363844, 13.39335497,  8.40372694,  9.72356215, -1.18922357])
```

## Recent, Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_3change_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3change_contemp.add_size_param("n_hist")
model_inf_3change_contemp.add_size_param("n_rec")
model_inf_3change_contemp.add_size_param("n_alb")
model_inf_3change_contemp.add_size_param("n_cont")
model_inf_3change_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3change_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3change_contemp.add_time_param("t_bot",upper=1e2)
model_inf_3change_contemp.add_leaf("CNas",N="n_cont")
model_inf_3change_contemp.set_size("CNas", N="n_alb", t="t_bot")
model_inf_3change_contemp.set_size("CNas", N="n_rec", t="t_rec")
model_inf_3change_contemp.set_size("CNas", N="n_hist", t="t_exp")
#run model
model_inf_3change_contemp.optimize(method="TNC")

            fun: 0.05270685135529537
            jac: array([ 1.48726364e-05, -3.22194574e-06, -1.48617254e-05, -7.81175982e-05,
       -1.01859004e-05,  7.57903841e-06,  2.87890071e-05])
  kl_divergence: 0.05270685135529537
 log_likelihood: -159999.29484808113
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 100
            nit: 17
     parameters: ParamsDict({'n_hist': 691060.4541214148, 'n_rec': 85318658.74539937, 'n_alb': 7666.880844694758, 'n_cont': 602.6370463331224, 't_exp': 594621.7664494055, 't_rec': 202.31070194110282, 't_bot': 65.2534170804026})
         status: 1
        success: True
              x: array([13.44598259, 18.26190373,  8.94466514,  6.4013151 ,  0.36614454,
       -1.18178232,  0.63019718])
```

## Recent, Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_3change_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3change_temporal.add_size_param("n_hist")
model_inf_3change_temporal.add_size_param("n_rec")
model_inf_3change_temporal.add_size_param("n_alb")
model_inf_3change_temporal.add_size_param("n_cont")
model_inf_3change_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3change_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3change_temporal.add_time_param("t_bot",upper=1e2)
model_inf_3change_temporal.add_leaf("CNas",N="n_cont")
model_inf_3change_temporal.add_leaf("AHam", N="n_alb",t=109)
model_inf_3change_temporal.move_lineages("CNas","AHam",t=110)
model_inf_3change_temporal.set_size("CNas", N="n_alb", t="t_bot")
model_inf_3change_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_3change_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_3change_temporal.optimize(method="TNC")

            fun: 0.24799423581098984
            jac: array([ 2.50382083e-06, -1.24842515e-06,  1.92658730e-06, -5.33851247e-06,
        3.68872649e-05, -1.14700119e-05, -1.89127854e-05])
  kl_divergence: 0.24799423581098984
 log_likelihood: -226767.54846730616
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 26
            nit: 6
     parameters: ParamsDict({'n_hist': 688429.7273753138, 'n_rec': 24749783.440175988, 'n_alb': 2257.771871529254, 'n_cont': 7455.377347092129, 't_exp': 633521.4702046391, 't_rec': 338.80641620857307, 't_bot': 79.04876765940242})
         status: 1
        success: True
              x: array([13.44216853, 17.0243273 ,  7.72213371,  8.91669084,  0.53144326,
        0.34589019,  1.32786751])
```

## Recent exponential change, contemp only 

```bash
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
#specify model
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_expg_contemp.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_expg_contemp.add_size_param("n_alb")
model_inf_expg_contemp.add_size_param("n_bot")
model_inf_expg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_expg_contemp.add_leaf("CNas",N="n_bot",g=lambda params: log(params.n_bot/params.n_alb)/params.t_bot) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_expg_contemp.set_size("CNas",g=0, t="t_bot")
#run model
model_inf_expg_contemp.optimize(method="TNC")

            fun: 0.13788519624466078
            fun: 0.13790137131752625
            jac: array([-5.34087518e-07, -5.49075742e-08,  1.08214497e-06])
  kl_divergence: 0.13790137131752625
 log_likelihood: -163852.38740241295
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 20
            nit: 8
     parameters: ParamsDict({'n_alb': 867220.4679051142, 'n_bot': 924.3978661939286, 't_bot': 95.2426112984719})
         status: 1
        success: True
              x: array([13.67304851,  6.82914257,  2.99672851])
```

## Recent exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data
model_inf_expg_temporal.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_expg_temporal.add_size_param("n_alb")
model_inf_expg_temporal.add_size_param("n_bot")
model_inf_expg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_expg_temporal.add_leaf("CNas",N="n_bot",g=lambda params: log(params.n_bot/params.n_alb)/params.t_bot)
model_inf_expg_temporal.set_size("CNas",g=0, t="t_bot")
model_inf_expg_temporal.add_leaf("AHam",N="n_alb",t=109)
model_inf_expg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_expg_temporal.optimize(method="TNC")

            fun: 0.33948021907659653
            jac: array([-6.59322041e-06, -1.02167783e-05,  1.23587792e-05])
  kl_divergence: 0.33948021907659653
 log_likelihood: -231050.92220380186
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 9
     parameters: ParamsDict({'n_alb': 873724.6520386358, 'n_bot': 97.93530937432844, 't_bot': 20.19633316538813})
         status: 1
        success: True
              x: array([13.68052056,  4.58430715, -1.37406839])
```

## Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_contemp.add_size_param("n_rec")
model_inf_recexpg_contemp.add_size_param("n_cont")
model_inf_recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_contemp.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_contemp.set_size("CNas", g=0, t="t_rec")
#run model
model_inf_recexpg_contemp.optimize(method="TNC")

            fun: 0.13790648046294465
            jac: array([ 5.27491429e-08, -5.60472639e-05,  3.95937956e-06])
  kl_divergence: 0.13790648046294465
 log_likelihood: -163852.6184737328
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 80
            nit: 22
     parameters: ParamsDict({'n_rec': 867339.7711163707, 'n_cont': 1158.2863620374997, 't_rec': 115.95582534731766})
         status: 1
        success: True
              x: array([13.67318607,  7.05469692, -4.35019386])
```

## Pre-Albatross exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recexpg_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temponly.add_size_param("n_rec")
model_inf_recexpg_temponly.add_size_param("n_cont")
model_inf_recexpg_temponly.add_time_param("t_rec",upper=390)
model_inf_recexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temponly.set_size("AHam", g=0, t="t_rec")
#run model
model_inf_recexpg_temponly.optimize(method="TNC")

            fun: 0.0649762381434033
            jac: array([-3.81066318e-10, -2.19564617e-06, -2.11457605e-07])
  kl_divergence: 0.0649762381434033
 log_likelihood: -49019.848749734905
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 33
            nit: 11
     parameters: ParamsDict({'n_rec': 913885.002275881, 'n_cont': 10000000000.000004, 't_rec': 389.5137840082854})
         status: 1
        success: True
              x: array([13.72546002, 23.02585093,  6.68600158])
```

## Pre-Albatross size change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_recexpg_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temporal.add_size_param("n_rec")
model_inf_recexpg_temporal.add_size_param("n_cont")
model_inf_recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_temporal.add_leaf("CNas",N="n_cont")
model_inf_recexpg_temporal.add_leaf("AHam", N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_recexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_recexpg_temporal.optimize(method="TNC")

            fun: 0.33950675719904866
            jac: array([-5.00161944e-06, -5.74868486e-04,  1.09525306e-06])
  kl_divergence: 0.33950675719904866
 log_likelihood: -231052.16471869507
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 11
     parameters: ParamsDict({'n_rec': 874649.8191328878, 'n_cont': 4899.514190557941, 't_rec': 112.08697567326021})
         status: 1
        success: True
              x: array([13.68157888,  8.49689133, -5.87738192])
```
		  
## Historic exponential change, contemp only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histexpg_contemp.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_contemp.add_size_param("n_hist")
model_inf_histexpg_contemp.add_size_param("n_cont")
model_inf_histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_contemp.add_leaf("CNas",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_contemp.set_size("CNas",g=0, t="t_exp")
#run model
model_inf_histexpg_contemp.optimize(method="TNC")

            fun: 0.14165794698984685
            jac: array([-7.0643568e-07,  1.6634691e-06,  1.5737754e-07])
  kl_divergence: 0.14165794698984685
 log_likelihood: -164022.286050345
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 54
            nit: 17
     parameters: ParamsDict({'n_hist': 891802.091434614, 'n_cont': 193220.0186846684, 't_exp': 10000.431629427841})
         status: 1
        success: True
              x: array([ 13.70099952,  12.17158481, -14.64564765])
 ```

## Historic exponential change, temporal only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_histexpg_temponly.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temponly.add_size_param("n_hist")
model_inf_histexpg_temponly.add_size_param("n_cont")
model_inf_histexpg_temponly.add_time_param("t_exp",lower=9890,upper=99890)
model_inf_histexpg_temponly.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_temponly.set_size("AHam",g=0, t="t_exp")
#run model
model_inf_histexpg_temponly.optimize(method="TNC")

            fun: 0.034323417903482555
            jac: array([ 8.08293625e-09, -2.38313413e-04, -4.69734190e-06])
  kl_divergence: 0.034323417903482555
 log_likelihood: -48196.78987347279
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 27
            nit: 10
     parameters: ParamsDict({'n_hist': 804827.6226805545, 'n_cont': 10000000000.000004, 't_exp': 99864.85032957778})
         status: 1
        success: True
              x: array([13.5983834 , 23.02585093,  8.18244068])
 ```

## Historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data
model_inf_histexpg_temporal.set_data(sfs,length=485532)
#set parameterss to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temporal.add_size_param("n_hist")
model_inf_histexpg_temporal.add_size_param("n_cont")
model_inf_histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_temporal.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp)
model_inf_histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_histexpg_temporal.add_leaf("CNas",N="n_cont")
model_inf_histexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_histexpg_temporal.optimize(method="TNC")

            fun: 0.3637995294167631
            jac: array([ 1.82953242e-06, -4.01778600e-05,  3.19318723e-06])
  kl_divergence: 0.3637995294167631
 log_likelihood: -232189.55231392846
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 39
            nit: 15
     parameters: ParamsDict({'n_hist': 900457.6242542629, 'n_cont': 179766.4094874783, 't_exp': 10005.768926499019})
         status: 1
        success: True
              x: array([ 13.71065838,  12.09941356, -12.05296838])
```

## Recent and Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_contemp.add_size_param("n_rec")
model_inf_2recexpg_contemp.add_size_param("n_alb")
model_inf_2recexpg_contemp.add_size_param("n_cont")
model_inf_2recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_contemp.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_contemp.set_size("CNas", g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot), t= "t_bot")
model_inf_2recexpg_contemp.set_size("CNas", g=0, t="t_rec")
#run model
model_inf_2recexpg_contemp.optimize(method="TNC")

            fun: 0.13788317625904126
            jac: array([ 1.82988272e-09, -2.57469257e-07,  1.53611913e-08, -1.49683103e-06,
        2.07970247e-06])
  kl_divergence: 0.13788317625904126
 log_likelihood: -163851.56449450285
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 71
            nit: 16
     parameters: ParamsDict({'n_rec': 867064.9452697256, 'n_alb': 10000000000.000004, 'n_cont': 167.40551645989237, 't_rec': 257.672585814103, 't_bot': 45.544838329219814})
         status: 1
        success: True
              x: array([13.67286916, 23.02585093,  5.12041911, -0.50208697, -0.17868034])
```

## Recent and Pre-Albatross exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2recexpg_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_temporal.add_size_param("n_rec")
model_inf_2recexpg_temporal.add_size_param("n_alb")
model_inf_2recexpg_temporal.add_size_param("n_cont")
model_inf_2recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot))
model_inf_2recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_2recexpg_temporal.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_temporal.set_size("CNas", g=0, t="t_bot")
model_inf_2recexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2recexpg_temporal.optimize(method="TNC")

            fun: 0.33944764272301975
            jac: array([-7.84194798e-07,  1.39140248e-05,  1.00121585e-05,  1.17098044e-05,
        3.42218171e-06])
  kl_divergence: 0.33944764272301975
 log_likelihood: -231049.3969789274
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 57
            nit: 13
     parameters: ParamsDict({'n_rec': 295117.96516554896, 'n_alb': 27226.251524356143, 'n_cont': 1253.197818725793, 't_rec': 275.5856658563322, 't_bot': 86.26160254978956})
         status: 1
        success: True
              x: array([12.59513044, 10.21193692,  7.13345382, -0.31006285,  1.83718992])
```

## Recent and historic size change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2changeexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2changeexpg_contemp.add_size_param("n_hist")
model_inf_2changeexpg_contemp.add_size_param("n_alb")
model_inf_2changeexpg_contemp.add_size_param("n_cont")
model_inf_2changeexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_contemp.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_contemp.set_size("CNas",g=lambda params: log(params.n_alb/params.n_hist)/(params.t_exp-params.t_bot), t= "t_bot")
model_inf_2changeexpg_contemp.set_size("CNas",g=0,t="t_exp")
#run model
model_inf_2changeexpg_contemp.optimize(method="TNC")
            
            fun: 0.05261924840020995
            jac: array([ 7.13901280e-06, -6.10583990e-06,  9.83508152e-06, -1.70583138e-05,
       -8.53739246e-06])
  kl_divergence: 0.05261924840020995
 log_likelihood: -159995.33282923148
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 8
     parameters: ParamsDict({'n_hist': 678667.7659156057, 'n_alb': 42360343.53783952, 'n_cont': 10.333042154447817, 't_exp': 804608.0920460667, 't_bot': 16.80318729771402})
         status: 1
        success: True
              x: array([13.42788699, 17.56172319,  2.33534674,  1.4028417 , -1.59964045])
```
			  
## Recent and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data
model_inf_2changeexpg_temporal.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_2changeexpg_temporal.add_size_param("n_alb")
model_inf_2changeexpg_temporal.add_size_param("n_hist")
model_inf_2changeexpg_temporal.add_size_param("n_cont")
model_inf_2changeexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_hist)/(params.t_exp-params.t_bot))
model_inf_2changeexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2changeexpg_temporal.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_temporal.set_size("CNas",g=0, t="t_bot")
model_inf_2changeexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2changeexpg_temporal.optimize(method="TNC")

            fun: 0.280467730713091
            jac: array([-3.93144603e-08, -1.43296073e-07, -2.01662618e-07, -1.06672785e-06,
        5.80776765e-07])
  kl_divergence: 0.280467730713091
 log_likelihood: -228287.95749862253
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 20
            nit: 6
     parameters: ParamsDict({'n_alb': 2156291.9163691313, 'n_hist': 643367.683678645, 'n_cont': 288.7594899925959, 't_exp': 999922.1624270558, 't_bot': 96.18977892634287})
         status: 1
        success: True
              x: array([14.5839006 , 13.37447167,  5.66559413,  9.45075734,  3.22863589])
```
			  
## Pre-Albatross and historic exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_contemp.add_size_param("n_hist")
model_inf_2histexpg_contemp.add_size_param("n_rec")
model_inf_2histexpg_contemp.add_size_param("n_cont")
model_inf_2histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_contemp.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_contemp.set_size("CNas",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_2histexpg_contemp.set_size("CNas", g=0, t="t_exp")
#run model
model_inf_2histexpg_contemp.optimize(method="TNC")

            fun: 0.052687657765304576
            jac: array([-1.20816253e-06,  3.38375771e-06, -2.19908550e-06,  6.61280977e-07,
        1.27179056e-06])
  kl_divergence: 0.052687657765304576
 log_likelihood: -159998.42677958662
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 37
            nit: 11
     parameters: ParamsDict({'n_hist': 689918.6912297576, 'n_rec': 9823562459.40747, 'n_cont': 177.30060493487022, 't_exp': 661972.0504549508, 't_rec': 395.69338771625814})
         status: 1
        success: True
              x: array([13.44432903, 23.00804967,  5.17784663,  0.65687311,  1.00407801])
```

## Pre-Albatross and historic exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_temponly.set_data(sfs, length=485532)
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

            fun: 0.03431687880717364
            jac: array([-5.12397241e-09, -2.38063055e-04, -3.02194453e-08, -7.79088118e-08,
       -1.49108460e-06])
  kl_divergence: 0.03431687880717364
 log_likelihood: -48196.6142921978
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 34
            nit: 9
     parameters: ParamsDict({'n_hist': 804821.6057821573, 'n_rec': 10000000000.000004, 'n_cont': 89912886.58251509, 't_exp': 99889.58281953415, 't_rec': 86.94753946701071})
         status: 1
        success: True
              x: array([13.59837592, 23.02585093, 18.31435183, 12.28179669, -1.24860099])
```

## Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_2histexpg_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temporal.add_size_param("n_hist")
model_inf_2histexpg_temporal.add_size_param("n_rec")
model_inf_2histexpg_temporal.add_size_param("n_cont")
model_inf_2histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_temporal.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_2histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2histexpg_temporal.add_leaf("CNas", N="n_cont",t=109)
model_inf_2histexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2histexpg_temporal.optimize(method="TNC")

            fun: 0.3584850070975374
            jac: array([ 6.06395279e-07,  4.53382796e-06, -1.93541919e-06,  7.28628923e-09,
       -9.24545092e-06])
  kl_divergence: 0.3584850070975374
 log_likelihood: -231940.72637894232
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 14
     parameters: ParamsDict({'n_hist': 876411.8203478971, 'n_rec': 491348.0788098974, 'n_cont': 3166.3386374873926, 't_exp': 10000.069458256643, 't_rec': 498.4517336832406})
         status: 1
        success: True
              x: array([ 13.68359137,  13.10490807,   8.0603312 , -16.47248948,
         5.52245548])
```

## Recent, Pre-Albatross and historical exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_3expg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_3expg_contemp.add_size_param("n_hist")
model_inf_3expg_contemp.add_size_param("n_rec")
model_inf_3expg_contemp.add_size_param("n_alb")
model_inf_3expg_contemp.add_size_param("n_cont")
model_inf_3expg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_3expg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_3expg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_3expg_contemp.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_contemp.set_size("CNas", g=lambda params: log(params.n_alb/params.n_rec)/(params.t_rec-params.t_bot), t= "t_bot")
model_inf_3expg_contemp.set_size("CNas",g=lambda params: log(params.n_rec/params.n_hist)/(params.t_exp-params.t_rec), t= "t_rec")
model_inf_3expg_contemp.set_size("CNas", g=0, t="t_exp")
#run model
model_inf_3expg_contemp.optimize(method="TNC")

            fun: 0.05262800949817243
            jac: array([ 1.17478214e-06,  1.04000565e-06, -3.27499840e-06,  3.54847470e-06,
        5.04901650e-08,  3.38095411e-07,  9.68253318e-07])
  kl_divergence: 0.05262800949817243
 log_likelihood: -159995.72906740903
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 45
            nit: 10
     parameters: ParamsDict({'n_hist': 678954.1873187566, 'n_rec': 43213079.372557364, 'n_alb': 2895.57661948658, 'n_cont': 384.70784992859, 't_exp': 805153.8753861334, 't_rec': 357.1717520216322, 't_bot': 87.23005313299866})
         status: 1
        success: True
              x: array([13.42830893, 17.58165377,  7.97093955,  5.95248421,  1.40632551,
        0.54438663,  1.92145441])
```

## Recent, Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=3.5e-9)
#add data to model
model_inf_3expg_temporal.set_data(sfs, length=485532)
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
model_inf_3expg_temporal.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_temporal.set_size("CNas", g=0, t="t_bot")
model_inf_3expg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_3expg_temporal.optimize(method="TNC")

            fun: 0.2763345425490948
            jac: array([ 1.84668836e-06, -3.96561468e-06, -4.05901318e-07,  1.71244010e-06,
       -4.02642346e-09, -4.63970392e-08,  5.31072978e-07])
  kl_divergence: 0.2763345425490948
 log_likelihood: -228094.44162878423
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 21
            nit: 5
     parameters: ParamsDict({'n_hist': 231172.5455384937, 'n_rec': 829692.6321847352, 'n_alb': 3519.2649024134785, 'n_cont': 2910.692947394513, 't_exp': 999999.6904749826, 't_rec': 499.9980926757859, 't_bot': 77.76897672450497})
         status: 1
        success: True
              x: array([12.35091966, 13.62881059,  8.16600741,  7.97614646, 14.97817627,
       12.22562839,  1.25225384])
```
