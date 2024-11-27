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
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CNas",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.15475256184257447
            jac: array([3.33339239e-11])
  kl_divergence: 0.15475256184257447
 log_likelihood: -164614.5161962893
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 4
     parameters: ParamsDict({'n_constant': 120741.66909753464})
         status: 1
        success: True
              x: array([11.70140858])
```

## Constant pop size, temporal only

```bash
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_tempnly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.06514603518730182
            jac: array([6.78378909e-10])
  kl_divergence: 0.06514603518730182
 log_likelihood: -49024.407970160624
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 4
     parameters: ParamsDict({'n_constant': 128060.49227027716})
         status: 1
        success: True
              x: array([11.76025803])
```

## Constant pop size, temp & contemp

```bash
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=485532)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CNas",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)#adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CNas","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.37876440437973063
            jac: array([2.83952086e-08])
  kl_divergence: 0.37876440437973063
 log_likelihood: -232890.2077596946
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 11
            nit: 4
     parameters: ParamsDict({'n_constant': 120268.3739416653})
         status: 0
        success: True
              x: array([11.69748097])
```

## Recent size change, contemp only

```bash
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.13834618646793775
            jac: array([-7.06194006e-06, -7.53990287e-06,  1.92967261e-04])
  kl_divergence: 0.13834618646793775
 log_likelihood: -163872.5050572206
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 16
            nit: 6
     parameters: ParamsDict({'n_alb': 121816.29190115027, 'n_bot': 5234.815346592085, 't_bot': 58.5777481139065})
         status: 1
        success: True
              x: array([11.71026938,  8.56308685,  0.34653668])
```

## Recent size change, temp and contemp

```bash
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

			fun: 0.33944634418784714
            jac: array([-4.29685607e-09,  1.70280363e-08,  1.34576895e-07])
  kl_divergence: 0.33944634418784714
 log_likelihood: -231049.33618151062
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 94
            nit: 29
     parameters: ParamsDict({'n_alb': 122296.81527686241, 'n_bot': 1.0, 't_bot': 0.016042659439470232})
         status: 1
        success: True
              x: array([11.71420628,  0.        , -8.73751364])
```

## Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.1387590290133308
            jac: array([-9.42460229e-06, -1.45114077e-05,  7.18907503e-06])
  kl_divergence: 0.1387590290133308
 log_likelihood: -163891.1766870211
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 17
            nit: 5
     parameters: ParamsDict({'n_rec': 122188.84217740857, 'n_cont': 9628.140488885167, 't_rec': 111.92764307827672})
         status: 1
        success: True
              x: array([11.71332301,  9.17244539, -6.03630004])
```
		  
## Pre-Albatross size change, temporal only

```bash
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temponly.add_size_param("n_rec")
model_inf_recchange_temponly.add_size_param("n_cont")
model_inf_recchange_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_recchange_temponly.set_size("AHam", N="n_rec", t="t_rec")
#run model
model_inf_recchange_temponly.optimize(method="TNC")

            fun: 0.06279349918227647
            jac: array([-1.12941177e-06, -2.97007682e-08, -4.28005277e-05])
  kl_divergence: 0.06279349918227647
 log_likelihood: -48961.24002588969
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 29
            nit: 11
     parameters: ParamsDict({'n_rec': 126478.3020143099, 'n_cont': 10000000000.000004, 't_rec': 490.7083765892394})
         status: 1
        success: True
              x: array([11.74782605, 23.02585093,  3.71029024])
```

## Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temporal.add_size_param("n_rec")
model_inf_recchange_temporal.add_size_param("n_cont")
model_inf_recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_recchange_temporal.add_leaf("AHam", N=n_cont", t=109)
model_inf_recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_recchange_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_recchange_temporal.optimize(method="TNC")

            fun: 0.340463873210023
            jac: array([ 3.72698323e-07, -4.00865503e-08,  9.83665375e-07])
  kl_divergence: 0.340463873210023
 log_likelihood: -231096.9768903289
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 55
            nit: 20
     parameters: ParamsDict({'n_rec': 123296.72812146386, 'n_cont': 6503.049586885535, 't_rec': 111.11718063431152})
         status: 1
        success: True
              x: array([11.72234915,  8.78002651, -8.10731671])
```

## Historic size change, contemp only

```bash
#specify model
model_inf_histchange_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
			    
            fun: 0.11916073211424216
            jac: array([ 7.52692140e-06, -1.17853663e-06,  1.42667698e-04])
  kl_divergence: 0.11916073211424216
 log_likelihood: -163004.80451316602
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 96
            nit: 32
     parameters: ParamsDict({'n_hist': 4234.8370388839985, 'n_cont': 203417.63202119607, 't_exp': 358806.73577822506})
         status: 1
        success: True
              x: array([ 8.35110013, 12.22301645, -0.60881291])
```
			  
## Historic size change, temporal only

```bash
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_histchange_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_histchange_temponly.add_size_param("n_hist")
model_inf_histchange_temponly.add_size_param("n_cont")
model_inf_histchange_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histchange_temponly.add_leaf("AHam",N="n_cont")
model_inf_histchange_temponly.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_histchange_temponly.optimize(method="TNC")
			    
            fun: 0.0035947942906817926
            jac: array([ 9.63844580e-05, -2.95107730e-05,  3.56397723e-07])
  kl_divergence: 0.0035947942906817926
 log_likelihood: -47371.69560084548
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 38
            nit: 14
     parameters: ParamsDict({'n_hist': 90843.35730911845, 'n_cont': 321961.02430462907, 't_exp': 113189.88910275811})
         status: 1
        success: True
              x: array([11.41689195, 12.68218577, -2.15106   ])
```

## Historic size change, temp and contemp		  

```bash
#specify model
model_inf_histchange_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
			    
            fun: 0.3425292009776663
            jac: array([ 4.94254458e-06,  6.48868030e-06, -2.55145390e-05])
  kl_divergence: 0.3425292009776663
 log_likelihood: -231193.67553640995
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 96
            nit: 27
     parameters: ParamsDict({'n_hist': 10337.591255370244, 'n_cont': 205026.3998527552, 't_exp': 343759.27705411235})
         status: 1
        success: True
              x: array([ 9.24354217, 12.23089403, -0.67610767])
```

## Recent and Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.13713251528566153
            jac: array([ 5.78017087e-11, -7.96128631e-09,  1.25980340e-07, -1.49110499e-07,
       -4.40626442e-09])
  kl_divergence: 0.13713251528566153
 log_likelihood: -163817.6143506598
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 82
            nit: 20
     parameters: ParamsDict({'n_rec': 120102.70060607034, 'n_alb': 10000000000.000004, 'n_cont': 1.0, 't_rec': 499.9437589553411, 't_bot': 0.014292578658204703})
         status: 1
        success: True
              x: array([11.69610249, 23.02585093,  0.        ,  8.84154321, -8.8530421 ])
```

## Recent and Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.3392501746093374
            jac: array([-1.20645146e-07,  3.14980559e-06,  8.49957971e-07, -4.62854652e-63,
       -3.09167279e-07])
  kl_divergence: 0.3392501746093374
 log_likelihood: -231040.1515218448
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 4
            nit: 1
     parameters: ParamsDict({'n_rec': 121920.05911118022, 'n_alb': 281958.065453314, 'n_cont': 3.2332026152582074, 't_rec': 500.0, 't_bot': 0.05671336363635188})
         status: 1
        success: True
              x: array([ 11.71112086,  12.54951363,   1.17347317, 135.33817936,
        -7.4743483 ])
```

## Recent and historic size change, contemp only

```bash
#specify model
model_inf_2change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.05272793742996014
            jac: array([-1.25417348e-05, -1.68472335e-06, -2.39271515e-05,  1.41805554e-05,
        4.47449209e-05])
  kl_divergence: 0.05272793742996014
 log_likelihood: -160000.24850798
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 51
            nit: 12
     parameters: ParamsDict({'n_hist': 96757.89120876072, 'n_alb': 12229345.10983416, 'n_cont': 234.68914587945088, 't_exp': 59774.94336713778, 't_bot': 21.262829837321732})
         status: 1
        success: True
              x: array([11.47996717, 16.31934896,  5.45826186, -2.93860756, -1.30915488])
```

## Recent and historic size change, temp and contemp

```bash
#specify model
model_inf_2change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2change_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2change_temporal.add_size_param("n_hist")
model_inf_2change_temporal.add_size_param("n_alb"
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

            fun: 0.26921801796673384
            jac: array([-3.45616986e-06,  1.68216928e-06, -3.83111556e-07,  7.56288848e-07,
        2.18962090e-05])
  kl_divergence: 0.26921801796673384
 log_likelihood: -227761.2459478381
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_hist': 83780.51569952827, 'n_alb': 292334.9561390578, 'n_cont': 142.3305860674932, 't_exp': 122410.84653078076, 't_bot': 3.612156040376977})
         status: 1
        success: True
              x: array([11.33595575, 12.58565553,  4.95815242, -2.05501812, -3.28407526])
```

## Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.052862551324543976
            jac: array([-4.77668426e-06, -5.21754995e-09,  7.47355669e-06, -1.53581565e-05,
        1.55028267e-07])
  kl_divergence: 0.052862551324543976
 log_likelihood: -160006.33669059034
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 80
            nit: 20
     parameters: ParamsDict({'n_hist': 96796.39344801658, 'n_rec': 10000000000.000004, 'n_cont': 1191.186508607882, 't_exp': 59700.433750993325, 't_rec': 111.1074208554139})
         status: 1
        success: True
              x: array([11.48036501, 23.02585093,  7.08270516, -2.94018485, -8.19430409])
```

## Pre-Albatross and historic size change, temporal only

```bash
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_temponly.set_data(sfs, length=485532)
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

            fun: 0.002952494201889239
            jac: array([ 8.22299508e-07,  2.37432450e-09, -9.70468307e-07, -4.57187288e-07,
        5.85982320e-08])
  kl_divergence: 0.002952494201889239
 log_likelihood: -47354.44920116131
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 42
            nit: 10
     parameters: ParamsDict({'n_hist': 96008.63860971128, 'n_rec': 68325436.13785093, 'n_cont': 1479.0886510236305, 't_exp': 77678.8327906517, 't_rec': 151.39101242121495})
         status: 1
        success: True
              x: array([11.47219345, 18.03979267,  7.2991814 , -2.61212003, -2.15534362])
```

## Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histchange_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histchange_temporal.add_size_param("n_hist")
model_inf_2histchange_temporal.add_size_param("n_rec")
model_inf_2histchange_temporal.add_size_param("n_cont")
model_inf_2histchange_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histchange_temporal.add_leaf("CNas",N="n_cont")
model_inf_2histchange_temporal.add_leaf("AHam", N="n_cont",t=109)
model_inf_2histchange_temporal.move_lineages("CNas","AHam",t=110)
model_inf_2histchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_2histchange_temporal.set_size("AHam", N="n_hist", t="t_exp")
#run model
model_inf_2histchange_temporal.optimize(method="TNC")

            fun: 0.24912100777893859
            jac: array([-1.84850834e-07, -2.30629830e-08, -1.05702397e-07, -6.39252516e-08,
       -1.78394589e-07])
  kl_divergence: 0.24912100777893859
 log_likelihood: -226820.30393084552
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 98
            nit: 24
     parameters: ParamsDict({'n_hist': 96581.8352802453, 'n_rec': 1666887.4817127457, 'n_cont': 6014.629678166622, 't_exp': 66711.92332913587, 't_rec': 499.72189607875646})
         status: 1
        success: True
              x: array([11.47814596, 14.32646866,  8.70195006, -2.80072944,  7.24262459])
```

## Recent, Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.052696313952353434
            jac: array([-7.23441125e-06,  1.17379435e-08,  3.53460059e-06,  1.31535815e-05,
       -3.90069390e-06, -2.72008855e-07, -1.28629935e-05])
  kl_divergence: 0.052696313952353434
 log_likelihood: -159998.81827295828
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 100
            nit: 15
     parameters: ParamsDict({'n_hist': 96758.42478939319, 'n_rec': 26819996.6949608, 'n_alb': 3600.1997020335652, 'n_cont': 1.9593220247884295, 't_exp': 59539.81153860802, 't_rec': 125.99796667046965, 't_bot': 0.11144833045270708})
         status: 1
        success: True
              x: array([11.47997269, 17.10465831,  8.1887446 ,  0.67259851, -2.9435927 ,
       -3.2163466 , -6.79824928])
```

## Recent, Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.24912083188323825
            jac: array([-1.55100696e-08, -3.93612874e-07,  6.30669530e-08,  3.19103695e-13,
       -2.43246266e-06, -1.89539429e-09,  1.47368083e-10])
  kl_divergence: 0.24912083188323825
 log_likelihood: -226820.29569540883
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 5
            nit: 1
     parameters: ParamsDict({'n_hist': 96582.65315318438, 'n_rec': 1671134.0483899915, 'n_alb': 6015.020947844817, 'n_cont': 2773475.070201748, 't_exp': 66696.88994824042, 't_rec': 499.9970220377002, 't_bot': 2.664320267868887e-05})
         status: 1
        success: True
              x: array([ 11.47815443,  14.32901302,   8.70201511,  14.83561163,
        -2.80101067,  11.78008769, -15.13814642])
```

## Recent exponential change, contemp only 

```bash
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
#specify model
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.13790008351945102
            jac: array([-3.06060106e-06,  3.57533913e-06,  1.55027091e-05])
  kl_divergence: 0.13790008351945102
 log_likelihood: -163852.3291591694
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 20
            nit: 6
     parameters: ParamsDict({'n_alb': 121406.38591723084, 'n_bot': 120.9146008702186, 't_bot': 9.019429049111494})
         status: 1
        success: True
              x: array([11.70689876,  4.79508452, -2.31126494])
```

## Recent exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.3394753064447676
            jac: array([ 1.24360859e-06, -8.53186653e-05,  1.15288499e-04])
  kl_divergence: 0.3394753064447676
 log_likelihood: -231050.69219437963
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 74
            nit: 19
     parameters: ParamsDict({'n_alb': 122347.46673919637, 'n_bot': 112.3419366209314, 't_bot': 12.725898508489117})
         status: 1
        success: True
              x: array([11.71462036,  4.72154723, -1.92541459])
```

## Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.13822238921582958
            jac: array([-9.21408293e-08, -9.89894593e-06,  5.67345621e-08])
  kl_divergence: 0.13822238921582958
 log_likelihood: -163866.9060788995
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 55
            nit: 20
     parameters: ParamsDict({'n_rec': 121768.08503020415, 'n_cont': 2367.758384328284, 't_rec': 111.01546944437655})
         status: 1
        success: True
              x: array([ 11.70987357,   7.76969896, -10.13242811])
```

## Pre-Albatross exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_recexpg_temponly.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recexpg_temponly.add_size_param("n_rec")
model_inf_recexpg_temponly.add_size_param("n_cont")
model_inf_recexpg_temponly.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recexpg_temponly.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_recexpg_temponly.set_size("AHam", g=0, t="t_rec")
#run model
model_inf_recexpg_temponly.optimize(method="TNC")

            fun: 0.06295933457864036
            jac: array([ 6.69885892e-11, -1.86417944e-05, -2.82775685e-08])
  kl_divergence: 0.06295933457864036
 log_likelihood: -48965.692872117455
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 44
            nit: 16
     parameters: ParamsDict({'n_rec': 126578.50970086116, 'n_cont': 10000000000.000004, 't_rec': 499.9934387818833})
         status: 1
        success: True
              x: array([11.74861802, 23.02585093, 10.99014148])
```

## Pre-Albatross size change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.3406784659358727
            jac: array([2.42270892e-05, 1.42259127e-04, 3.55121313e-05])
  kl_divergence: 0.3406784659358727
 log_likelihood: -231107.02412175317
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 13
     parameters: ParamsDict({'n_rec': 123617.57524677253, 'n_cont': 6834.352249907571, 't_rec': 119.09905059741142})
         status: 1
        success: True
              x: array([11.72494801,  8.82971698, -3.85079252])
```
		  
## Historic exponential change, contemp only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.1259024163754658
            jac: array([ 3.01505677e-07,  9.35448922e-08, -4.66004476e-06])
  kl_divergence: 0.1259024163754658
 log_likelihood: -163309.71066724838
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 40
            nit: 14
     parameters: ParamsDict({'n_hist': 27609.624355585274, 'n_cont': 205073.6016040889, 't_exp': 790243.4382875398})
         status: 1
        success: True
              x: array([10.2259197 , 12.23112423,  1.31365834])
 ```

## Historic exponential change, temporal only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_histexpg_temponly.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_histexpg_temponly.add_size_param("n_hist")
model_inf_histexpg_temponly.add_size_param("n_cont")
model_inf_histexpg_temponly.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_histexpg_temponly.add_leaf("AHam",N="n_cont",g=lambda params: log(params.n_cont/params.n_hist)/params.t_exp) #parameterizes exp growth rate in terms of starting and ending pop sizes
model_inf_histexpg_temponly.set_size("AHam",g=0, t="t_exp")
#run model
model_inf_histexpg_temponly.optimize(method="TNC")

            fun: 0.004111123269165564
            jac: array([-4.71156196e-06,  1.52002504e-04,  2.27518780e-05])
  kl_divergence: 0.004111123269165564
 log_likelihood: -47385.55955024675
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 59
            nit: 20
     parameters: ParamsDict({'n_hist': 76641.08098057321, 'n_cont': 381972.80538918934, 't_exp': 246061.1553413577}
         status: 1
        success: True
              x: array([11.24688852, 12.85310469, -1.16122035])
 ```

## Historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.3491548685161667
            jac: array([ 4.46150558e-06, -2.43047467e-07, -1.37010418e-06])
  kl_divergence: 0.3491548685161667
 log_likelihood: -231503.88929056254
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 60
            nit: 19
     parameters: ParamsDict({'n_hist': 26268.127902763066, 'n_cont': 204074.932999037, 't_exp': 808782.5677983493})
         status: 1
        success: True
              x: array([10.17611162, 12.22624252,  1.42967761])
```

## Recent and Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_contemp.add_size_param("n_rec")
model_inf_2recexpg_contemp.add_size_param("n_alb")
model_inf_2recexpg_contemp.add_size_param("n_cont")
model_inf_2recexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_contemp.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_contemp.set_size("CNas", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec, t= "t_bot")
model_inf_2recexpg_contemp.set_size("CNas", g=0, t="t_rec")
#run model
model_inf_2recexpg_contemp.optimize(method="TNC")

            fun: 0.13754717724927998
            jac: array([-4.26041456e-06, -9.69757615e-06, -1.28019885e-05, -1.57444884e-04,
        1.66419883e-05])
  kl_divergence: 0.13754717724927998
 log_likelihood: -163836.36826728837
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 72
            nit: 16
     parameters: ParamsDict({'n_rec': 10447.239677337022, 'n_alb': 10000000000.000004, 'n_cont': 288.04246830622253, 't_rec': 374.33419952355416, 't_bot': 66.52080118695791})
         status: 1
        success: True
              x: array([ 9.25409308, 23.02585093,  5.66310793,  0.73979794,  0.68659038])
```

## Recent and Pre-Albatross exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2recexpg_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2recexpg_temporal.add_size_param("n_rec")
model_inf_2recexpg_temporal.add_size_param("n_alb")
model_inf_2recexpg_temporal.add_size_param("n_cont")
model_inf_2recexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2recexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2recexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec)
model_inf_2recexpg_temporal.set_size("AHam", g=0, t="t_rec")
model_inf_2recexpg_temporal.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2recexpg_temporal.set_size("CNas", g=0, t="t_bot")
model_inf_2recexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2recexpg_temporal.optimize(method="TNC")

            fun: 0.33946517946458304
            jac: array([5.95361772e-06, 5.63086782e-07, 3.85334094e-06, 3.56745581e-06,
       2.48417582e-07])
  kl_divergence: 0.33946517946458304
 log_likelihood: -231050.2180491674
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 7
            nit: 2
     parameters: ParamsDict({'n_rec': 122289.82096165365, 'n_alb': 9883168247.074593, 'n_cont': 14.478106903327356, 't_rec': 291.5915826497615, 't_bot': 5.178603422151369})
         status: 1
        success: True
              x: array([11.71414909, 23.01409897,  2.67263764, -0.14326167, -2.90745968])
```

## Recent and historic size change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2changeexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, alb size, historic size (pre-alb), times of two size changes
model_inf_2changeexpg_contemp.add_size_param("n_hist")
model_inf_2changeexpg_contemp.add_size_param("n_alb")
model_inf_2changeexpg_contemp.add_size_param("n_cont")
model_inf_2changeexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_contemp.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_contemp.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_contemp.set_size("CNas",g=lambda params: log(params.n_alb/params.n_hist)/params.t_exp, t= "t_bot")
model_inf_2changeexpg_contemp.set_size("CNas",g=0,t="t_exp")
model_inf_2changeexpg_contemp.optimize(method="TNC")
            
            fun: 0.05262438215965355
            jac: array([ 7.10343128e-08, -1.33077797e-07,  2.02571643e-07,  1.97214035e-07,
        4.04132743e-06])
  kl_divergence: 0.05262438215965355
 log_likelihood: -159995.56501376984
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 83
            nit: 20
     parameters: ParamsDict({'n_hist': 94885.92329912631, 'n_alb': 5935404.20983619, 'n_cont': 32.124195669467966, 't_exp': 81005.21439920632, 't_bot': 29.877608806567558})
         status: 1
        success: True
              x: array([11.46043064, 15.59644569,  3.46960951, -2.56052713, -0.85313283])
```
			  
## Recent and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data
model_inf_2changeexpg_temporal.set_data(sfs,length=485532)
#set parameters to infer - contemp size, alb size, time of bottleneck
model_inf_2changeexpg_temporal.add_size_param("n_alb")
model_inf_2changeexpg_temporal.add_size_param("n_hist")
model_inf_2changeexpg_temporal.add_size_param("n_cont")
model_inf_2changeexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2changeexpg_temporal.add_time_param("t_bot",upper=1e2)
model_inf_2changeexpg_temporal.add_leaf("AHam",N="n_alb",g=lambda params: log(params.n_alb/params.n_hist)/params.t_exp)
model_inf_2changeexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2changeexpg_temporal.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_2changeexpg_temporal.set_size("CNas",g=0, t="t_bot")
model_inf_2changeexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2changeexpg_temporal.optimize(method="TNC")

            fun: 0.2729530068213786
            jac: array([ 5.26198901e-06, -3.70802384e-06, -3.58280526e-06, -6.68263238e-06,
        4.04778026e-05])
  kl_divergence: 0.2729530068213786
 log_likelihood: -227936.11812601256
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 55
            nit: 15
     parameters: ParamsDict({'n_alb': 311899.92426598864, 'n_hist': 17863.087844854643, 'n_cont': 377.72434791950275, 't_exp': 629882.1641845881, 't_bot': 65.75490804774162})
         status: 1
        success: True
              x: array([12.65043766,  9.79049173,  5.93416469,  0.51570797,  0.65239106])
```
			  
## Pre-Albatross and historic exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_contemp.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_contemp.add_size_param("n_hist")
model_inf_2histexpg_contemp.add_size_param("n_rec")
model_inf_2histexpg_contemp.add_size_param("n_cont")
model_inf_2histexpg_contemp.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_contemp.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_contemp.add_leaf("CNas",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_contemp.set_size("CNas",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_2histexpg_contemp.set_size("CNas", g=0, t="t_exp")
#run model
model_inf_2histexpg_contemp.optimize(method="TNC")

            fun: 0.052671522831912225
            jac: array([-1.93717557e-06, -3.50133808e-06, -9.72880940e-07,  1.58774268e-05,
        2.11774084e-05])
  kl_divergence: 0.052671522831912225
 log_likelihood: -159997.69704495408
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 65
            nit: 15
     parameters: ParamsDict({'n_hist': 93971.21738728328, 'n_rec': 6586080.121323865, 'n_cont': 305.8839126536867, 't_exp': 80295.8247566267, 't_rec': 236.89579737840785})
         status: 1
        success: True
              x: array([11.45074382, 15.70046891,  5.72320566, -2.57133966, -0.7370956 ])
```

## Pre-Albatross and historic exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_temponly.set_data(sfs, length=485532)
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

            fun: 0.0029536962398617407
            jac: array([ 5.55177118e-06,  9.25592779e-07,  5.88473205e-06,  4.55963941e-05,
       -1.35224871e-06])
  kl_divergence: 0.0029536962398617407
 log_likelihood: -47354.48147708291
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 99
            nit: 23
     parameters: ParamsDict({'n_hist': 93454.3787594531, 'n_rec': 23163338.76379235, 'n_cont': 397.8480806800519, 't_exp': 94396.28274751532, 't_rec': 387.1028312550731})
         status: 1
        success: True
              x: array([11.44522867, 16.95808136,  5.98607023, -2.37307845,  0.89429598])
```

## Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_2histexpg_temporal.set_data(sfs, length=485532)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_2histexpg_temporal.add_size_param("n_hist")
model_inf_2histexpg_temporal.add_size_param("n_rec")
model_inf_2histexpg_temporal.add_size_param("n_cont")
model_inf_2histexpg_temporal.add_time_param("t_exp",lower=1e4,upper=1e6)
model_inf_2histexpg_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_2histexpg_temporal.add_leaf("AHam",N="n_cont", g=lambda params: log(params.n_cont/params.n_rec)/params.t_rec)
model_inf_2histexpg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_2histexpg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_2histexpg_temporal.add_leaf("CNas", N="n_cont",t=109)
model_inf_2histexpg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_2histexpg_temporal.optimize(method="TNC")

            fun: 0.37136378948915577
            jac: array([-1.09193225e-05, -5.94595442e-07,  5.02916654e-07,  8.80219606e-08,
        8.14884029e-08])
  kl_divergence: 0.37136378948915577
 log_likelihood: -232543.7109705179
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 60
            nit: 15
     parameters: ParamsDict({'n_hist': 123655.32766326301, 'n_rec': 116695.9745042532, 'n_cont': 1244.0519702602262, 't_exp': 10000.569618362695, 't_rec': 111.00202835816303})
         status: 1
        success: True
              x: array([ 11.72525336,  11.66732732,   7.12612905, -14.36824833,
       -12.16410273])
```

## Recent, Pre-Albatross and historical exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
model_inf_3expg_contemp.set_size("CNas", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec, t= "t_bot")
model_inf_3expg_contemp.set_size("CNas",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_3expg_contemp.set_size("CNas", g=0, t="t_exp")
#run model
model_inf_3expg_contemp.optimize(method="TNC")

            fun: 0.05265461569347702
            jac: array([ 2.82443192e-04, -1.21511769e-04, -1.27157269e-04, -4.50866990e-04,
       -1.74320806e-04,  2.11096257e-04, -2.82571895e-05])
  kl_divergence: 0.05265461569347702
 log_likelihood: -159996.93238580407
        message: 'Max. number of function evaluations reached'
           nfev: 101
            nit: 16
     parameters: ParamsDict({'n_hist': 2023028.7411477377, 'n_rec': 184389611.1957335, 'n_alb': 5241.412647472242, 'n_cont': 184.66824724901718, 't_exp': 78363.00217560692, 't_rec': 166.8629597543654, 't_bot': 48.83347997528906})
         status: 3
        success: False
              x: array([14.52010632, 19.03256153,  8.56434633,  5.21856096, -2.60131966,
       -1.78565239, -0.04666927])
```

## Recent, Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
model_inf_3expg_temporal.add_leaf("AHam",N="n_alb", g=lambda params: log(params.n_alb/params.n_rec)/params.t_rec)
model_inf_3expg_temporal.set_size("AHam",g=lambda params: log(params.n_rec/params.n_hist)/params.t_exp, t= "t_rec")
model_inf_3expg_temporal.set_size("AHam",g=0, t="t_exp")
model_inf_3expg_temporal.add_leaf("CNas", N="n_cont", g=lambda params: log(params.n_cont/params.n_alb)/params.t_bot)
model_inf_3expg_temporal.set_size("CNas", g=0, t="t_bot")
model_inf_3expg_temporal.move_lineages("CNas","AHam",t=110)
#run model
model_inf_3expg_temporal.optimize(method="TNC")

            fun: 0.26664362364496597
            jac: array([ 3.55219329e-06, -3.80490873e-06, -2.21927336e-06,  4.92730326e-07,
        4.21761429e-06, -2.36565767e-06,  9.16277033e-06])
  kl_divergence: 0.26664362364496597
 log_likelihood: -227640.71280569292
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 4
            nit: 1
     parameters: ParamsDict({'n_hist': 13793.777288216152, 'n_rec': 338980.6199575828, 'n_alb': 5484.051265069044, 'n_cont': 4398.2840496717145, 't_exp': 631796.1571001245, 't_rec': 499.9052505728709, 't_bot': 79.95732700221512})
         status: 1
        success: True
              x: array([ 9.53197285, 12.73369822,  8.60959939,  8.38896976,  0.52397561,
        8.31985522,  1.38362943])
```
