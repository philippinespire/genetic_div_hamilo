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
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CNas",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.1547525618425751
            jac: array([3.0360747e-14])
  kl_divergence: 0.1547525618425751
 log_likelihood: -164614.51619628933
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 14
            nit: 4
     parameters: ParamsDict({'n_constant': 120741.66906623017})
         status: 1
        success: True
              x: array([11.70140858])
```

## Constant pop size, temporal only

```bash
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=485532) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_tempnly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.06514603518730182
            jac: array([6.15955195e-14])
  kl_divergence: 0.06514603518730182
 log_likelihood: -49024.407970160624
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 12
            nit: 4
     parameters: ParamsDict({'n_constant': 128060.49187888313})
         status: 1
        success: True
              x: array([11.76025802])
```

## Constant pop size, temp & contemp

```bash
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=485532)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CNas",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)#adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CNas","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.37910435048787217
            jac: array([3.65922151e-14])
  kl_divergence: 0.37910435048787217
 log_likelihood: -232906.1240364778
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 10
            nit: 4
     parameters: ParamsDict({'n_constant': 120576.95587570632})
         status: 1
        success: True
              x: array([11.70004347])
```

## Recent size change, contemp only

```bash
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.13789445892566868
            jac: array([-4.01507011e-05,  2.80658074e-04, -2.64120285e-04])
  kl_divergence: 0.13789445892566868
 log_likelihood: -163852.0747756664
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 41
            nit: 15
     parameters: ParamsDict({'n_alb': 121354.61496795954, 'n_bot': 130.1861360000737, 't_bot': 1.9199311189861523})
         status: 1
        success: True
              x: array([11.70647224,  4.86896524, -3.93349486])
```

## Recent size change, temp and contemp

```bash
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.3396215969237677
            jac: array([-5.52154128e-06, -3.19607455e-06,  1.28478593e-04])
  kl_divergence: 0.3396215969237677
 log_likelihood: -231057.54151460642
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 18
            nit: 5
     parameters: ParamsDict({'n_alb': 122460.57351754051, 'n_bot': 1100.77306629389, 't_bot': 24.94527607404725})
         status: 1
        success: True
              x: array([11.71554441,  7.003768  , -1.10153303])
```

## Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.13851078189692717
            jac: array([ 2.36501155e-07, -2.53012818e-05,  1.05139518e-07])
  kl_divergence: 0.13851078189692717
 log_likelihood: -163879.9492146875
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 61
            nit: 20
     parameters: ParamsDict({'n_rec': 121973.33730330408, 'n_cont': 7009.696468658677, 't_rec': 111.01818001038886})
         status: 1
        success: True
              x: array([11.71155775,  8.85504968, -9.97096523])
```
		  
## Pre-Albatross size change, temporal only

```bash
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.06415746441724292
            jac: array([-6.04395102e-05, -1.52049597e-05, -3.28209072e-04])
  kl_divergence: 0.06415746441724292
 log_likelihood: -48997.86385641377
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 5
     parameters: ParamsDict({'n_rec': 127347.72178550586, 'n_cont': 8405787.12315958, 't_rec': 289.05532341559933})
         status: 1
        success: True
              x: array([11.75467659, 15.94443097, -0.1695016 ])
```

## Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.3401975499302251
            jac: array([-2.81828817e-05, -2.67615943e-04,  1.69679961e-05])
  kl_divergence: 0.3401975499302251
 log_likelihood: -231084.50763436875
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 13
     parameters: ParamsDict({'n_rec': 122999.92693123457, 'n_cont': 4737.79444968098, 't_rec': 112.33761800796496})
         status: 1
        success: True
              x: array([11.71993904,  8.463327  , -5.66924439])
```

## Historic size change, contemp only

```bash
#specify model
model_inf_histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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
			    
            fun: 0.11917196629467183
            jac: array([ 4.85369513e-05, -6.83071733e-04,  1.26655604e-04])
  kl_divergence: 0.11917196629467183
 log_likelihood: -163005.3126014443
        message: 'Max. number of function evaluations reached'
           nfev: 101
            nit: 28
     parameters: ParamsDict({'n_hist': 18853.07816239883, 'n_cont': 202209.72716178323, 't_exp': 454416.0620205225})
         status: 3
        success: False
              x: array([ 9.84443148, 12.21706069, -0.20509547])
```
			  
## Historic size change, temporal only

```bash
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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
			    
            fun: 0.0035947557542230363
            jac: array([ 5.42455164e-09, -7.03216942e-09, -1.85735612e-08])
  kl_divergence: 0.0035947557542230363
 log_likelihood: -47371.694566103026
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 57
            nit: 17
     parameters: ParamsDict({'n_hist': 90761.32372061227, 'n_cont': 321857.21995896683, 't_exp': 157838.68837137494})
         status: 0
        success: True
              x: array([11.41598852, 12.68186331, -1.73984984])
```

## Historic size change, temp and contemp		  

```bash
#specify model
model_inf_histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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
			    
            fun: 0.3427121834883651
            jac: array([ 4.14251036e-06, -4.66318717e-06, -9.86911692e-05])
  kl_divergence: 0.3427121834883651
 log_likelihood: -231202.24277756087
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 59
            nit: 19
     parameters: ParamsDict({'n_hist': 27017.65761394764, 'n_cont': 205667.5128634407, 't_exp': 426822.6880108057})
         status: 1
        success: True
              x: array([10.20424592, 12.23401613, -0.31853419])
```

## Recent and Pre-Albatross size change, contemp only

```bash
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.13851070659745993
            jac: array([ 4.31611127e-12,  0.00000000e+00, -5.13068758e-11,  3.86364767e-08,
        0.00000000e+00])
  kl_divergence: 0.13851070659745993
 log_likelihood: -163879.9458091185
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 63
            nit: 16
     parameters: ParamsDict({'n_rec': 121973.72091334919, 'n_alb': 31079.97001474941, 'n_cont': 7016.069247755969, 't_rec': 111.00693840039123, 't_bot': 2.718447830549104})
         status: 1
        success: True
              x: array([ 11.7115609 ,  10.34431884,   8.8559584 , -10.93424553,
        -3.57754831])
```

## Recent and Pre-Albatross size change, temp and contemp

```bash
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

fun: 0.33934000499806277
            jac: array([-1.59555553e-06, -2.30566015e-07,  1.62212298e-07, -1.10498761e-08,
        1.73033744e-09])
  kl_divergence: 0.33934000499806277
 log_likelihood: -231044.3573806449
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 69
            nit: 16
     parameters: ParamsDict({'n_rec': 122074.77574591296, 'n_alb': 293787.41995804, 'n_cont': 1.0, 't_rec': 499.97256351242515, 't_bot': 0.02402041532077487})
         status: 1
        success: True
              x: array([11.71238905, 12.59061172,  0.        ,  9.5593903 , -8.33378112])
```

## Recent and historic size change, contemp only

```bash
#specify model
model_inf_2change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.05278546886717413
            jac: array([ 5.81097612e-06, -1.73001085e-06,  6.52464179e-06, -2.06524153e-05,
        1.50713924e-05])
  kl_divergence: 0.05278546886717413
 log_likelihood: -160002.85048229087
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 83
            nit: 21
     parameters: ParamsDict({'n_hist': 96782.88824576537, 'n_alb': 17736102.226764422, 'n_cont': 644.5205753046723, 't_exp': 83194.10013364904, 't_bot': 82.0816103558371})
         status: 1
        success: True
              x: array([11.48022548, 16.69111279,  6.46850675, -2.52778096,  1.52188646])
```

## Recent and historic size change, temp and contemp

```bash
#specify model
model_inf_2change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

			fun: 0.26911677491626107
            jac: array([-2.01934077e-06,  1.05289240e-06, -4.19036927e-07, -3.61293796e-07,
        1.23332577e-05])
  kl_divergence: 0.26911677491626107
 log_likelihood: -227756.50574821496
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_hist': 83741.87272269317, 'n_alb': 292163.3814096943, 'n_cont': 78.06315130280231, 't_exp': 170239.58015337546, 't_bot': 2.761466947060552})
         status: 1
        success: True
              x: array([11.3354944 , 12.58506845,  4.35751813, -1.64446694, -3.56140502])
```

## Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.052815129630512485
            jac: array([-1.56221658e-05, -1.14551316e-08, -2.64226945e-06,  8.57795051e-06,
        9.64853141e-09])
  kl_divergence: 0.052815129630512485
 log_likelihood: -160004.19194963438
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 86
            nit: 18
     parameters: ParamsDict({'n_hist': 96783.48722353263, 'n_rec': 3386011885.464452, 'n_cont': 855.1699905794678, 't_exp': 83107.75607286434, 't_rec': 111.00871560777078})
         status: 1
        success: True
              x: array([ 11.48023167,  21.94291863,   6.75130027,  -2.52905549,  -10.7061968 ])
```

## Pre-Albatross and historic size change, temporal only

```bash
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.003604115363886547
            jac: array([ 9.81757208e-06, -6.96921818e-07,  8.91950673e-07,  1.90309305e-07,
        7.49477471e-07])
  kl_divergence: 0.003604115363886547
 log_likelihood: -47371.9458809821
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 64
            nit: 15
     parameters: ParamsDict({'n_hist': 90738.17733601239, 'n_rec': 321056.5673283817, 'n_cont': 3720585.1336476207, 't_exp': 158114.02143142815, 't_rec': 120.88094881613168})
         status: 1
        success: True
              x: array([11.41573347, 12.67937261, 15.12939151, -1.7376622 , -3.64724173])
```

## Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.24877688401259593
            jac: array([-1.11249725e-06,  4.14608964e-07, -8.42285806e-07,  5.92393681e-06,
       -4.24617228e-07])
  kl_divergence: 0.24877688401259593
 log_likelihood: -226804.19205610536
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 13
            nit: 4
     parameters: ParamsDict({'n_hist': 96392.78850669235, 'n_rec': 1611082.737968213, 'n_cont': 4367.848543495222, 't_exp': 92520.0364055123, 't_rec': 499.34823652567314})
         status: 1
        success: True
              x: array([11.47618667, 14.29241702,  8.38202584, -2.39763036,  6.38997601])
```

## Recent, Pre-Albatross and historic size change, contemp only

```bash
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.05272923222038136
            jac: array([ 2.17849627e-06, -5.30945154e-08, -1.97942028e-06, -7.08605845e-06,
       -1.31349563e-05,  8.18900538e-07,  2.87071004e-05])
  kl_divergence: 0.05272923222038136
 log_likelihood: -160000.30706746638
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 4
            nit: 1
     parameters: ParamsDict({'n_hist': 96770.73363183733, 'n_rec': 278710750.57179576, 'n_alb': 33102.91178208107, 'n_cont': 245.01528435414053, 't_exp': 82859.37194334481, 't_rec': 225.1396175633434, 't_bot': 30.25438461879669})
         status: 1
        success: True
              x: array([11.48009989, 19.44568507, 10.40737653,  5.50132059, -2.53272965,
       -0.87884085, -0.83521344])
```

## Recent, Pre-Albatross and historic size change, temp and contemp

```bash
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.24872023084119066
            jac: array([-2.72742581e-05, -8.73114009e-06, -1.73707855e-05,  1.07479962e-05,
       -1.38578969e-05,  9.21418316e-06, -4.46450821e-06])
  kl_divergence: 0.24872023084119066
 log_likelihood: -226801.53955462016
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 9
            nit: 2
     parameters: ParamsDict({'n_hist': 96692.00717412314, 'n_rec': 6068986.414345329, 'n_alb': 3292.239287595851, 'n_cont': 4645.212381701367, 't_exp': 88971.23797343168, 't_rec': 463.06602588469616, 't_bot': 96.93916473378927})
         status: 1
        success: True
              x: array([11.47928602, 15.61870217,  8.09932325,  8.44359237, -2.44549076,
        2.2546869 ,  3.45539577])
```

## Recent exponential change, contemp only 

```bash
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
#specify model
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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
            jac: array([-8.18881287e-08, -3.69757133e-08,  3.91855817e-06])
  kl_divergence: 0.13788519624466078
 log_likelihood: -163851.65585239246
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 27
            nit: 7
     parameters: ParamsDict({'n_alb': 121393.50489082932, 'n_bot': 23.05143311412777, 't_bot': 2.951523752548452})
         status: 1
        success: True
              x: array([11.70679265,  3.13772794, -3.49288905])
```

## Recent exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.3394610411293201
            jac: array([ 3.93219019e-05,  2.74512552e-04, -2.87381148e-04])
  kl_divergence: 0.3394610411293201
 log_likelihood: -231050.02429231038
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 51
            nit: 15
     parameters: ParamsDict({'n_alb': 122301.87381768516, 'n_bot': 27.297879608801146, 't_bot': 5.128294568394272})
         status: 1
        success: True
              x: array([11.71424764,  3.30680903, -2.91775235])
```

## Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.1381143356086217
            jac: array([-4.92839974e-06, -9.78236077e-06,  1.03107939e-06])
  kl_divergence: 0.1381143356086217
 log_likelihood: -163862.01913840632
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 51
            nit: 16
     parameters: ParamsDict({'n_rec': 121645.40156151234, 'n_cont': 1597.9864211816316, 't_rec': 111.40643143890436})
         status: 1
        success: True
              x: array([11.70886555,  7.37649963, -6.86287401])
```

## Pre-Albatross exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.06356983469496447
            jac: array([ 5.75073133e-09, -1.34919850e-05, -2.28799344e-08])
  kl_divergence: 0.06356983469496447
 log_likelihood: -48982.085410740874
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 26
            nit: 9
     parameters: ParamsDict({'n_rec': 126989.09764567083, 'n_cont': 10000000000.000004, 't_rec': 499.9926661231869})
         status: 1
        success: True
              x: array([11.75185652, 23.02585093, 10.8788115 ])
```

## Pre-Albatross size change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.340282599539089
            jac: array([-5.83875833e-06, -7.88779516e-06,  3.22415067e-06])
  kl_divergence: 0.340282599539089
 log_likelihood: -231088.48965705576
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 39
            nit: 15
     parameters: ParamsDict({'n_rec': 123204.29442291227, 'n_cont': 4910.364129359067, 't_rec': 111.89499209235495})
         status: 1
        success: True
              x: array([11.72159919,  8.49910338, -6.07221634])
```
		  
## Historic exponential change, contemp only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.12593973175180864
            jac: array([ 7.60009782e-06,  1.96210147e-05, -5.33533234e-05])
  kl_divergence: 0.12593973175180864
 log_likelihood: -163311.39832977424
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 17
     parameters: ParamsDict({'n_hist': 43034.94804421936, 'n_cont': 204803.21184262287, 't_exp': 848619.6207866713})
         status: 1
        success: True
              x: array([10.66976781, 12.22980485,  1.71196149])
 ```

## Historic exponential change, temporal only 

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.004110957064094276
            jac: array([ 1.40929822e-05,  3.92958648e-06, -4.15684747e-06])
  kl_divergence: 0.004110957064094276
 log_likelihood: -47385.55508747438
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 55
            nit: 18
     parameters: ParamsDict({'n_hist': 76770.67700322662, 'n_cont': 381347.3127510388, 't_exp': 341482.0468215655})
         status: 1
        success: True
              x: array([11.24857804, 12.85146582, -0.68641813])
 ```

## Historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.3755258926679366
            jac: array([4.22703485e-06, 5.06530135e-07, 5.92425831e-06])
  kl_divergence: 0.3755258926679366
 log_likelihood: -232738.5806413484
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 34
            nit: 12
     parameters: ParamsDict({'n_hist': 126851.66636996307, 'n_cont': 93995.08426995439, 't_exp': 10008.189125897496})
         status: 1
        success: True
              x: array([ 11.7507737 ,  11.45099776, -11.70264479])
```

## Recent and Pre-Albatross exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.13752988051213424
            jac: array([ 2.66265415e-08, -4.94722763e-06,  1.81833233e-08, -9.84903524e-16,
        1.35910524e-07])
  kl_divergence: 0.13752988051213424
 log_likelihood: -163835.58598775748
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 24
            nit: 7
     parameters: ParamsDict({'n_rec': 119610.22914613984, 'n_alb': 10000000000.000004, 'n_cont': 1.0000116841391484, 't_rec': 499.99999999923705, 't_bot': 0.4218158915397544})
         status: 1
        success: True
              x: array([ 1.16919936e+01,  2.30258509e+01,  1.16840709e-05,  2.69573504e+01,
       -5.46412944e+00])
```

## Recent and Pre-Albatross exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.3394686284630602
            jac: array([2.09499379e-06, 3.21705136e-07, 1.18081206e-06, 7.17743665e-07,
       1.54370551e-06])
  kl_divergence: 0.3394686284630602
 log_likelihood: -231050.3795312761
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 7
            nit: 2
     parameters: ParamsDict({'n_rec': 122304.40989775065, 'n_alb': 7226162513.042549, 'n_cont': 10.11312956469077, 't_rec': 113.9772383233002, 't_bot': 4.786565994458967})
         status: 1
        success: True
              x: array([11.71426838, 22.70097396,  2.31383454, -4.8649002 , -2.9903078 ])
```

## Recent and historic size change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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
#run model
model_inf_2changeexpg_contemp.optimize(method="TNC")
            
             fun: 0.05262929599553802
            jac: array([ 1.28357671e-06,  7.54952630e-07,  2.04752420e-06, -1.77997814e-07,
        2.29296880e-06])
  kl_divergence: 0.05262929599553802
 log_likelihood: -159995.7872518254
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 72
            nit: 19
     parameters: ParamsDict({'n_hist': 94800.69219240596, 'n_alb': 6097876.539936607, 'n_cont': 60.61506159415043, 't_exp': 112547.97145512915, 't_bot': 74.79778737402366})
         status: 1
        success: True
              x: array([11.45953199, 15.62345116,  4.1045434 , -2.15802376,  1.08785651])
```
			  
## Recent and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.2726859087877937
            jac: array([-1.71292808e-06, -3.76311007e-07,  2.65293417e-07, -3.54170483e-14,
       -2.15491722e-10])
  kl_divergence: 0.2726859087877937
 log_likelihood: -227923.61259608011
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 80
            nit: 19
     parameters: ParamsDict({'n_alb': 311074.87473361153, 'n_hist': 11889.347000374011, 'n_cont': 1.0, 't_exp': 999999.9762611754, 't_bot': 0.4525115130898972})
         status: 1
        success: True
              x: array([12.64778892,  9.38339807,  0.        , 17.5461036 , -5.39357687])```
```
			  
## Pre-Albatross and historic exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.05958069060600997
            jac: array([-4.60836559e-07,  5.16752946e-08, -3.50304091e-07, -1.18257122e-07,
        1.02644056e-05])
  kl_divergence: 0.05958069060600997
 log_likelihood: -160310.1779758732
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 42
            nit: 10
     parameters: ParamsDict({'n_hist': 5684.428167975763, 'n_rec': 442396.68480894505, 'n_cont': 338.34615916268586, 't_exp': 940842.0406154696, 't_rec': 128.90698853933412})
         status: 1
        success: True
              x: array([ 8.64548582, 12.99996224,  5.82406951,  2.75587845, -3.03126168])
```

## Pre-Albatross and historic exponential change, temporal only

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.0029522344482584358
            jac: array([ 4.65340272e-07, -1.40498601e-06, -1.20754013e-05, -7.91999734e-06,
        4.47450431e-06])
  kl_divergence: 0.0029522344482584358
 log_likelihood: -47354.44222651657
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 87
            nit: 20
     parameters: ParamsDict({'n_hist': 94246.7033562499, 'n_rec': 201226455.2173973, 'n_cont': 144.5980535850605, 't_exp': 123801.99055464036, 't_rec': 272.2299818157171})
         status: 1
        success: True
              x: array([11.45367113, 19.11994147,  4.97395785, -2.04113209, -0.34550462])
```

## Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.29741017730626634
            jac: array([-3.13899454e-06, -6.84291359e-07,  3.67461769e-07, -3.71097306e-07,
       -3.84268969e-07])
  kl_divergence: 0.29741017730626634
 log_likelihood: -229081.202848115
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 71
            nit: 18
     parameters: ParamsDict({'n_hist': 13361.95918342115, 'n_cont': 1398.4663588312144, 'n_rec': 268654.9076525111, 't_exp': 963626.3323445306, 't_rec': 499.9956076571407})
         status: 1
        success: True
              x: array([ 9.50016708,  7.24313146, 12.50118296,  3.26642681, 11.39146057])
```

## Recent, Pre-Albatross and historical exponential change, contemp only

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.05262990230905679
            jac: array([ 3.74881654e-07, -6.42678490e-08, -6.55885688e-07, -9.44483102e-10,
        8.97555453e-06, -1.33900930e-07,  2.05862564e-06])
  kl_divergence: 0.05262990230905679
 log_likelihood: -159995.8146735669
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 80
            nit: 15
     parameters: ParamsDict({'n_hist': 43773.5547570807, 'n_rec': 2723503.907659654, 'n_alb': 12041336.694217999, 'n_cont': 64.41376192691459, 't_exp': 112953.72618410317, 't_rec': 161.38998173939441, 't_bot': 83.52084293056308})
         status: 1
        success: True
              x: array([10.68678514, 14.81742981, 16.30385601,  4.1653273 , -2.15361752,
       -1.90505668,  1.62299984])
```

## Recent, Pre-Albatross and historic exponential change, temp and contemp

```bash
from autograd.numpy import log
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1.394, muts_per_gen=2.5e-8)
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

            fun: 0.26680906207074934
            jac: array([ 1.13526187e-06, -3.84284720e-06,  4.69733332e-07, -3.05343450e-07,
        7.22824531e-08, -5.69065640e-22,  4.14308966e-07])
  kl_divergence: 0.26680906207074934
 log_likelihood: -227648.4586327881
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 53
            nit: 14
     parameters: ParamsDict({'n_hist': 9724.185285083677, 'n_rec': 333077.1549282568, 'n_alb': 4179.385321860573, 'n_cont': 2.964572944203051, 't_exp': 981951.1421432283, 't_rec': 500.0, 't_bot': 0.09943209677151356})
         status: 1
        success: True
              x: array([ 9.18237139, 12.71612944,  8.33791946,  1.08673299,  3.98622313,
       44.24387856, -6.91245568])
```
