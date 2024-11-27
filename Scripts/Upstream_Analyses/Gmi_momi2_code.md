#### RUNNING MOMI2 ON WAHAB ####
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

##Constant pop size, contemp only
#specify model
model_inf_constant_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_contemp.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_contemp.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_contemp.add_leaf("CBat",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_contemp.optimize(method="TNC")

            fun: 0.2506224208032932
            jac: array([1.75456524e-09])
  kl_divergence: 0.2506224208032932
 log_likelihood: -18372.83344348638
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 7
            nit: 3
     parameters: ParamsDict({'n_constant': 16165.187024312596})
         status: 1
        success: True
              x: array([9.69061526])

##Constant pop size, temporal only
#specify model
model_inf_constant_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8) #this sets the model input
#add data
model_inf_constant_temponly.set_data(sfs, length=467359) #gives the data to simulate and # of SNPs that go into it
#set parameter to infer - contemp size
model_inf_constant_temponly.add_size_param("n_constant") #says, want to estimate n
model_inf_constant_tempnly.add_leaf("AHam",N="n_constant") #wants to estimate n in this leaf (population)
#run model
model_inf_constant_temponly.optimize(method="TNC")

            fun: 0.16220847010228948
            jac: array([1.81722797e-08])
  kl_divergence: 0.16220847010228948
 log_likelihood: -9636.23215911425
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 10
            nit: 3
     parameters: ParamsDict({'n_constant': 16023.938035203233})
         status: 0
        success: True
              x: array([9.68183901])

##Constant pop size, temp & contemp
#specify model 
model_inf_constant_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_constant_temporal.set_data(sfs, length=467359)
#set parameter to infer - contemp size
model_inf_constant_temporal.add_size_param("n_constant")
model_inf_constant_temporal.add_leaf("CBat",N="n_constant")
model_inf_constant_temporal.add_leaf("AHam",N="n_constant",t=109)#adds another population (leaf) at a specific time
model_inf_constant_temporal.move_lineages("CBat","AHam",t=110) #says move ALL indv from one pop to another at this time
#run model
model_inf_constant_temporal.optimize(method="TNC")

            fun: 0.7696632421752071
            jac: array([1.6781456e-08])
  kl_divergence: 0.7696632421752071
 log_likelihood: -31097.921819451207
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 9
            nit: 4
     parameters: ParamsDict({'n_constant': 15840.993766767238})
         status: 0
        success: True
              x: array([9.6703564])

##Recent size change, contemp only
#specify model
model_inf_change_contemp =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.2425199739037029
            jac: array([-3.08523481e-07, -8.47790489e-09, -2.77818790e-07])
  kl_divergence: 0.2425199739037029
 log_likelihood: -18316.391798383833
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 32
            nit: 11
     parameters: ParamsDict({'n_alb': 15097.22655230098, 'n_bot': 10000000000.000004, 't_bot': 99.99578469002915})
         status: 1
        success: True
              x: array([ 9.62226633, 23.02585093, 10.07416018])

##Recent size change, temp and contemp
#specify model
model_inf_change_temporal =  momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.7684898103203458
            jac: array([-6.83794877e-06,  1.67297770e-05, -5.20990089e-06])
  kl_divergence: 0.7684898103203458
 log_likelihood: -31089.158630359103
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 31
            nit: 9
     parameters: ParamsDict({'n_alb': 15444.172798208174, 'n_bot': 32670.92731603795, 't_bot': 99.15907444075927})
         status: 1
        success: True
              x: array([ 9.64498705, 10.39424089,  4.76997751])

##Pre-Albatross size change, contemp only
#specify model
model_inf_recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.22654886185976097
            jac: array([ 2.48961138e-05, -1.75049734e-05, -2.76545656e-06])
  kl_divergence: 0.22654886185976097
 log_likelihood: -18205.137031885733
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 35
            nit: 14
     parameters: ParamsDict({'n_rec': 12987.889061408505, 'n_cont': 43104.34749502873, 't_rec': 499.9289129897196})
         status: 1
        success: True
              x: array([ 9.47177259, 10.67137914,  8.60724724])

##Pre-Albatross size change, temporal only
#specify model
model_inf_recchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.12494447712486616
            jac: array([-5.12761069e-08, -5.96396091e-07, -1.97279655e-06])
  kl_divergence: 0.12494447712486616
 log_likelihood: -9463.55081565687
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 24
            nit: 9
     parameters: ParamsDict({'n_rec': 12628.279029959891, 'n_cont': 657386203.2316313, 't_rec': 499.9661009936787})
         status: 1
        success: True
              x: array([ 9.44369395, 20.30378223,  9.34786177])

##Pre-Albatross size change, temp and contemp
#specify model
model_inf_recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
#add data to model
model_inf_recchange_temporal.set_data(sfs, length=467359)
#set parameters to infer - contemp size, historic size (pre-alb), time of size changes
model_inf_recchange_temporal.add_size_param("n_rec")
model_inf_recchange_temporal.add_size_param("n_cont")
model_inf_recchange_temporal.add_time_param("t_rec",lower=111,upper=5e2)
model_inf_recchange_temporal.add_leaf("CBat",N="n_cont")
model_inf_recchange_temporal.add_leaf("AHam", N=n_cont", t=109)
model_inf_recchange_temporal.set_size("AHam", N="n_rec", t="t_rec")
model_inf_recchange_temporal.move_lineages("CBat","AHam",t=110)
#run model
model_inf_recchange_temporal.optimize(method="TNC")

            fun: 0.767800084989135
            jac: array([1.45082272e-07, 8.54720726e-09, 1.81240509e-07])
  kl_divergence: 0.767800084989135
 log_likelihood: -31084.00776158562
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 38
            nit: 13
     parameters: ParamsDict({'n_rec': 15397.262350731062, 'n_cont': 1192739352.5015318, 't_rec': 177.4406641774952})
         status: 1
        success: True
              x: array([ 9.641945  , 20.89951848, -1.57997783])

##Historic size change, contemp only
#specify model
model_inf_histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
			    
            fun: 0.2506224208032948
            jac: array([0.0000000e+00, 1.5014996e-08, 0.0000000e+00])
  kl_divergence: 0.2506224208032948
 log_likelihood: -18372.83344348639
        message: 'Local minimum reached (|pg| ~= 0)'
           nfev: 14
            nit: 4
     parameters: ParamsDict({'n_hist': 1001.0870563875804, 'n_cont': 16165.189072342851, 't_exp': 14696.249163052715})
         status: 0
        success: True
              x: array([ 6.90884174,  9.69061539, -5.34618583])

##Historic size change, temporal only
#specify model
model_inf_histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
			    
            fun: 0.015552432185710185
            jac: array([-6.56902725e-07, -9.64603101e-07, -6.83699969e-06])
  kl_divergence: 0.015552432185710185
 log_likelihood: -8956.628079408822
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 26
            nit: 9
     parameters: ParamsDict({'n_hist': 8351.923956326065, 'n_cont': 53612.64977742627, 't_exp': 19910.802595319445})
         status: 1
        success: True
              x: array([ 9.03024721, 10.88954032, -4.59401825])

##Historic size change, temp and contemp		  
#specify model
model_inf_histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
			    
            fun: 0.616050475240657
            jac: array([-3.57517511e-05,  1.65277726e-05,  7.31727008e-05])
  kl_divergence: 0.616050475240657
 log_likelihood: -29950.741675983987
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 24
            nit: 8
     parameters: ParamsDict({'n_hist': 7027.546799876036, 'n_cont': 53199.71664471141, 't_exp': 23531.55736417384})
         status: 1
        success: True
              x: array([ 8.85759296, 10.88180835, -4.27891789])
			  
##Recent and Pre-Albatross size change, contemp only
#specify model
model_inf_2recchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.21835952900316113
            jac: array([ 1.36615244e-06, -1.89179476e-08,  1.96080794e-06, -1.37952856e-16,
       -7.67112434e-07])
  kl_divergence: 0.21835952900316113
 log_likelihood: -18148.09013920666
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 86
            nit: 19
     parameters: ParamsDict({'n_rec': 10321.311278486823, 'n_alb': 10000000000.000004, 'n_cont': 5.8121765201238915, 't_rec': 499.9999999999981, 't_bot': 0.07732679401832793})
         status: 1
        success: True
              x: array([ 9.24196609, 23.02585093,  1.75995512, 32.96293489, -7.16411138])

##Recent and Pre-Albatross size change, temp and contemp
#specify model
model_inf_2recchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.7572102743615556
            jac: array([ 5.97514259e-08,  9.82667725e-08,  7.21539829e-08, -4.98858980e-33,
        1.50638417e-07])
  kl_divergence: 0.7572102743615556
 log_likelihood: -31004.923055818857
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 99
            nit: 24
     parameters: ParamsDict({'n_rec': 13256.873433905437, 'n_alb': 34866.31799275578, 'n_cont': 3.7265487518979814, 't_rec': 500.0, 't_bot': 0.02184808441553456})
         status: 1
        success: True
              x: array([ 9.49227145, 10.45927654,  1.31548254, 70.02923745, -8.42859371])
			  
##Recent and historic size change, contemp only
#specify model
model_inf_2change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.04903734535483152
            jac: array([2.13010416e-05, 5.28465144e-06, 9.01765628e-06, 1.28374611e-05,
       1.41348942e-05])
  kl_divergence: 0.04903734535483152
 log_likelihood: -16968.591807912395
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 29
            nit: 7
     parameters: ParamsDict({'n_hist': 9490.91434000129, 'n_alb': 771759.8181135856, 'n_cont': 957.0627868212428, 't_exp': 14941.405806292314, 't_bot': 75.38228319133373})
         status: 1
        success: True
              x: array([ 9.15809023, 13.55642866,  6.863869  , -5.29505126,  1.1191059 ])

##Recent and historic size change, temp and contemp
model_inf_2change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.6102557335010643
            jac: array([-3.53774093e-06, -1.22418575e-06, -1.41995675e-06,  8.98593405e-06,
        3.05384641e-05])
  kl_divergence: 0.6102557335010643
 log_likelihood: -29907.46654467271
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 92
            nit: 21
     parameters: ParamsDict({'n_hist': 7487.517117494892, 'n_alb': 57765.66012735797, 'n_cont': 943.8819573352519, 't_exp': 22546.714173565226, 't_bot': 7.0259438744215235})
         status: 1
        success: True
              x: array([ 8.92099253, 10.96414976,  6.85000111, -4.35549169, -2.58271092])

##Pre-Albatross and historic size change, contemp only
#specify model
model_inf_2histchange_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.049078560024754615
            jac: array([ 3.31864455e-07, -8.46546855e-10,  2.81523259e-07,  4.41101284e-06,
        1.34899236e-05])
  kl_divergence: 0.049078560024754615
 log_likelihood: -16968.87890930308
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 66
            nit: 16
     parameters: ParamsDict({'n_hist': 9496.455552862495, 'n_rec': 10000000000.000004, 'n_cont': 1458.3658507256755, 't_exp': 14831.512561679683, 't_rec': 128.76343938974517})
         status: 1
        success: True
              x: array([ 9.15867391, 23.02585093,  7.28507181, -5.3176531 , -3.03969711])

##Pre-Albatross and historic size change, temporal only
#specify model
model_inf_2histchange_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.014850762683237269
            jac: array([-6.47945431e-07,  8.10223017e-07, -1.50499383e-06, -4.00304460e-06,
        2.73576862e-06])
  kl_divergence: 0.014850762683237269
 log_likelihood: -8953.376542934362
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 66
            nit: 15
     parameters: ParamsDict({'n_hist': 8912.735465400281, 'n_rec': 66253.41178404902, 'n_cont': 9530.521868285698, 't_exp': 17895.513003662207, 't_rec': 184.43906682679813})
         status: 1
        success: True
              x: array([ 9.09523648, 11.10124224,  9.16225476, -4.82340308, -1.45789576])

##Pre-Albatross and historic size change, temp and contemp
#specify model
model_inf_2histchange_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.5887142038624853
            jac: array([-1.20814613e-07, -1.75552276e-06,  4.19191433e-06,  2.05196235e-07,
       -1.27904379e-07])
  kl_divergence: 0.5887142038624853
 log_likelihood: -29746.5944013318
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 93
            nit: 23
     parameters: ParamsDict({'n_hist': 9006.351042039356, 'n_rec': 100558.33180281635, 'n_cont': 11866.582311396021, 't_exp': 18295.06367421879, 't_rec': 499.99069561436403})
         status: 1
        success: True
              x: array([ 9.10568528, 11.51849325,  9.38148152, -4.77363019, 10.64082484])

##Recent, Pre-Albatross and historic size change, contemp only
#specify model
model_inf_3change_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.04900483811208187
            jac: array([-3.99820794e-05,  3.27448695e-06,  2.24950412e-11,  8.28687278e-05,
       -1.97520308e-07, -1.21004674e-08, -1.38108023e-05])
  kl_divergence: 0.04900483811208187
 log_likelihood: -16968.3653624594
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 81
            nit: 15
     parameters: ParamsDict({'n_hist': 9392.74276820225, 'n_rec': 294325.0226077957, 'n_alb': 2582598395.3152575, 'n_cont': 730.887761602192, 't_exp': 15294.828691802155, 't_rec': 368.87847554106077, 't_bot': 48.328503366735085})
         status: 1
        success: True
              x: array([ 9.14769262, 12.59243995, 21.67206186,  6.59425991, -5.22561165,
        0.67636389, -0.06688479])

##Recent, Pre-Albatross and historic size change, temp and contemp
#specify model
model_inf_3change_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

			fun: 0.5809993684152883
            jac: array([-3.50336185e-08, -1.25668741e-06, -1.04551347e-09, -4.44582607e-09,
        2.78917334e-06,  8.12552982e-07, -7.89041400e-09])
  kl_divergence: 0.5809993684152883
 log_likelihood: -29688.980010212133
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 21
            nit: 5
     parameters: ParamsDict({'n_hist': 9739.561236742396, 'n_rec': 33843854.860583164, 'n_alb': 4229.099029186769, 'n_cont': 4426231959.381977, 't_exp': 14880.739632266565, 't_rec': 499.4237876392127, 't_bot': 99.99986304431384})
         status: 1
        success: True
              x: array([ 9.18395135, 17.337268  ,  8.34974425, 22.21081449, -5.30746594,
        6.51337598, 13.50102196])

##Recent exponential change, contemp only 
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
model_inf_expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.24301935510238604
            jac: array([-3.26432646e-09, -3.72113555e-05, -4.51973485e-06])
  kl_divergence: 0.24301935510238604
 log_likelihood: -18319.87048781386
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 24
            nit: 8
     parameters: ParamsDict({'n_alb': 15136.956876782267, 'n_bot': 10000000000.000004, 't_bot': 99.9284957147097})
         status: 1
        success: True
              x: array([ 9.62489451, 23.02585093,  7.24245278])

###Recent exponential change, temp and contemp
from autograd.numpy import log
model_inf_expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.7688111009232282
            jac: array([-1.28216753e-07,  7.02812826e-07, -9.21253197e-07])
  kl_divergence: 0.7688111009232282
 log_likelihood: -31091.55802858143
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 27
            nit: 10
     parameters: ParamsDict({'n_alb': 15557.926286891963, 'n_bot': 894494097.6003059, 't_bot': 55.21621377608118})
         status: 1
        success: True
              x: array([ 9.65232552, 20.61176886,  0.20941048])
			  
##Pre-Albatross exponential change, contemp only
from autograd.numpy import log
model_inf_recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.23106195824597867
            jac: array([ 4.35130216e-06, -1.11648354e-06, -6.44028809e-07])
  kl_divergence: 0.23106195824597867
 log_likelihood: -18236.575261312126
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 34
            nit: 12
     parameters: ParamsDict({'n_rec': 14183.029903121771, 'n_cont': 173902.95337906881, 't_rec': 499.97274022677465})
         status: 1
        success: True
              x: array([ 9.55980145, 12.06625268,  9.56585244])

##Pre-Albatross exponential change, temporal only
from autograd.numpy import log
model_inf_recexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.12714954326387662
            jac: array([ 3.51200971e-09, -1.67538607e-04, -9.55834371e-07])
  kl_divergence: 0.12714954326387662
 log_likelihood: -9473.769092145045
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 37
            nit: 11
     parameters: ParamsDict({'n_rec': 12723.580764438773, 'n_cont': 10000000000.000004, 't_rec': 499.98306478063813})
         status: 1
        success: True
              x: array([ 9.4512123 , 23.02585093, 10.04189565])

##Pre-Albatross size change, temp and contemp
from autograd.numpy import log
model_inf_recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.7642337435428218
            jac: array([ 1.36905706e-07,  3.45129155e-06, -1.33289123e-06])
  kl_divergence: 0.7642337435428218
 log_likelihood: -31057.374323664553
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 40
            nit: 15
     parameters: ParamsDict({'n_rec': 14926.343796826617, 'n_cont': 31327.609988315973, 't_rec': 499.90342698801373})
         status: 1
        success: True
              x: array([ 9.61088297, 10.3522551 ,  8.30078701])
			  
##Historic exponential change, contemp only 
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
model_inf_histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.08419721015653103
            jac: array([-6.02938695e-06,  1.57594618e-06, -2.55480808e-05])
  kl_divergence: 0.08419721015653103
 log_likelihood: -17213.515426121034
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 30
            nit: 14
     parameters: ParamsDict({'n_hist': 1216.9062442341651, 'n_cont': 67742.30539432236, 't_exp': 72532.36255075423})
         status: 1
        success: True
              x: array([ 7.10406705, 11.12346616, -2.69677368])

##Historic exponential change, temporal only 
from autograd.numpy import log #otherwise won't recognize log function in model (can say np.log in growth function but that doesn't run right either??)
model_inf_histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.017447309934639943
            jac: array([-1.47567045e-06,  3.38609087e-06, -1.48392274e-05])
  kl_divergence: 0.017447309934639943
 log_likelihood: -8965.408942897362
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 56
            nit: 21
     parameters: ParamsDict({'n_hist': 1526.7929179733308, 'n_cont': 69239.86701225713, 't_exp': 66741.68842100815})
         status: 1
        success: True
              x: array([ 7.33092468, 11.14533209, -2.80017284])
 
###Historic exponential change, temp and contemp
from autograd.numpy import log
model_inf_histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.6326068232396757
            jac: array([ 4.29603488e-07, -9.54174317e-07,  1.35456896e-06])
  kl_divergence: 0.6326068232396757
 log_likelihood: -30074.384482840658
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 57
            nit: 18
     parameters: ParamsDict({'n_hist': 910.1810111082997, 'n_cont': 57693.07469830108, 't_exp': 86984.34703590035})
         status: 1
        success: True
              x: array([ 6.81364349, 10.96289242, -2.47315091])

##Recent and Pre-Albatross exponential change, contemp only
#specify model
model_inf_2recexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.24112521208547433
            jac: array([-0.00403567, -0.01233989,  0.01028191,  0.00513463, -0.02302775])
  kl_divergence: 0.24112521208547433
 log_likelihood: -18306.675887558053
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 77
            nit: 8
     parameters: ParamsDict({'n_rec': 15.911938205751282, 'n_alb': 10000000000.000004, 'n_cont': 9999999809.07156, 't_rec': 132.2570204315168, 't_bot': 44.26062198357669})
         status: 1
        success: True
              x: array([ 2.76706966, 23.02585093, 23.02585091, -2.85069706, -0.23059148])

##Recent and Pre-Albatross exponential change, temp and contemp
#specify model
model_inf_2recexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.7632236379509274
            jac: array([ 4.55052668e-09,  2.98617851e-08,  1.90472780e-08, -4.32340070e-22,
        4.86232982e-07])
  kl_divergence: 0.7632236379509274
 log_likelihood: -31049.830855104286
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 97
            nit: 22
     parameters: ParamsDict({'n_rec': 14704.208216420157, 'n_alb': 61220.677402732, 'n_cont': 6.127423800517535, 't_rec': 500.0, 't_bot': 0.2627321881212179})
         status: 1
        success: True
              x: array([ 9.59588901, 11.02224028,  1.8127744 , 44.23440042, -5.93915947])

##Recent and historic exponential size change, contemp only
#specify model
model_inf_2changeexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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
            
            fun: 0.04869863747238678
            jac: array([ 2.04965750e-07, -9.18906010e-07,  1.44008478e-07, -1.61801755e-06,
        3.48207713e-06])
  kl_divergence: 0.04869863747238678
 log_likelihood: -16966.232368803285
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 68
            nit: 16
     parameters: ParamsDict({'n_hist': 8872.338093070972, 'n_alb': 2827887.103759936, 'n_cont': 40.147695165238936, 't_exp': 19175.530252114644, 't_bot': 31.41940931901135})
         status: 1
        success: True
              x: array([ 9.09069364, 14.85504038,  3.69256503, -4.67185333, -0.78058373])

###Recent and historic exponential change, temp and contemp
from autograd.numpy import log
model_inf_2changeexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.6246235648226108
            jac: array([-1.03384228e-07, -5.28548008e-08,  6.15486516e-07, -3.40628491e-07,
       -1.12571203e-07])
  kl_divergence: 0.6246235648226108
 log_likelihood: -30014.76550898202
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 10
            nit: 3
     parameters: ParamsDict({'n_alb': 67389.94700510419, 'n_hist': 20.973133896042523, 'n_cont': 5.645018962046234, 't_exp': 155634.9858730654, 't_bot': 0.4566909819828753})
         status: 1
        success: True
              x: array([11.11825113,  3.04324228,  1.73077356, -1.75748149, -5.38434112])

##Pre-Albatross and historic exponential change, contemp only
#specify model
model_inf_2histexpg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.04873334285223303
            jac: array([-2.41416984e-07,  3.05645652e-07, -1.31995224e-07,  5.84964942e-08,
        1.55092106e-05])
  kl_divergence: 0.04873334285223303
 log_likelihood: -16966.474126479294
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 59
            nit: 14
     parameters: ParamsDict({'n_hist': 8410.42522840673, 'n_rec': 3048294.347072656, 'n_cont': 327.2730616372988, 't_exp': 19042.738627432576, 't_rec': 212.83025965869498})
         status: 1
        success: True
              x: array([ 9.03722731, 14.93009276,  5.79079487, -4.68656682, -1.03676616])

##Pre-Albatross and historic exponential change, temporal only
#specify model
model_inf_2histexpg_temponly = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.014855925896606705
            jac: array([-6.72245919e-07, -2.03270431e-06, -2.27019129e-07, -9.97924398e-06,
       -1.68376574e-07])
  kl_divergence: 0.014855925896606705
 log_likelihood: -8953.400469265116
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 53
            nit: 14
     parameters: ParamsDict({'n_hist': 7758.707099655648, 'n_rec': 146159.3861099045, 'n_cont': 1075.8155345954685, 't_exp': 27220.812970932537, 't_rec': 189.2615360839331})
         status: 1
        success: True
              x: array([ 8.95657099, 11.89245299,  6.98083429, -4.03403841, -1.37889536])

##Pre-Albatross and historic exponential change, temp and contemp
#specify model
model_inf_2histexpg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.5924238885810141
            jac: array([-5.82475277e-07, -2.68504110e-07,  6.28782572e-07, -5.16131485e-06,
       -7.73787465e-07])
  kl_divergence: 0.5924238885810141
 log_likelihood: -29774.298326809774
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 99
            nit: 23
     parameters: ParamsDict({'n_hist': 1.2840416566575636, 'n_rec': 90561.53208781501, 'n_cont': 2926.970031062013, 't_exp': 162774.10222560682, 't_rec': 499.9857731432342})
         status: 1
        success: True
              x: array([ 0.25001265, 11.41378481,  7.98172305, -1.70113355, 10.21616655])

##Recent, Pre-Albatross and historical exponential change, contemp only
#specify model
model_inf_3expg_contemp = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.04869768407579971
            jac: array([ 6.85549746e-07,  1.77776374e-06,  6.47255256e-08, -9.27646744e-08,
        5.20840822e-06,  1.43464917e-07,  2.82945822e-06])
  kl_divergence: 0.04869768407579971
 log_likelihood: -16966.22572744266
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 9
            nit: 2
     parameters: ParamsDict({'n_hist': 10565.734524957936, 'n_rec': 3414569.203428931, 'n_alb': 973212.1845126101, 'n_cont': 31.954787310107495, 't_exp': 19162.215528422676, 't_rec': 140.28375563332077, 't_bot': 23.14353595528605})
         status: 1
        success: True
              x: array([ 9.26537145, 15.04356189, 13.78835741,  3.46432201, -4.67331907,
       -2.50828256, -1.20022406])

##Recent, Pre-Albatross and historic exponential change, temp and contemp
#specify model
model_inf_3expg_temporal = momi.DemographicModel(N_e=NeConstant, gen_time=1, muts_per_gen=2.5e-8)
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

            fun: 0.5967828303899166
            jac: array([ 1.77892796e-06,  4.42980033e-07, -9.70554885e-08, -5.41348884e-05,
        1.45430997e-05, -2.05752830e-09, -4.41565619e-05])
  kl_divergence: 0.5967828303899166
 log_likelihood: -29806.850904238658
        message: 'Converged (|f_n-f_(n-1)| ~= 0)'
           nfev: 69
            nit: 15
     parameters: ParamsDict({'n_hist': 3.3230222016211797, 'n_rec': 92961.20326399336, 'n_alb': 3449.589114381881, 'n_bot': 10000000000.000004, 't_exp': 155060.49285569612, 't_rec': 499.99995953987764, 't_bot': 99.61507933173557})
         status: 1
        success: True
              x: array([ 1.20087467, 11.43993752,  8.14601041, 23.02585093, -1.76211419,
       16.07877294,  5.55603157])