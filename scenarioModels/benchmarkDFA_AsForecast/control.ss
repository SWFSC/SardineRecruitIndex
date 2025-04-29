#C file created using an r4ss function
#C file write time: 2025-04-24  17:20:57
#
1 # 0 means do not read wtatage.ss; 1 means read and usewtatage.ss and also read and use growth parameters
1 #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern
4 # recr_dist_method for parameters
1 # not yet implemented; Future usage:Spawner-Recruitment; 1=global; 2=by area
1 # number of recruitment settlement assignments 
0 # unused option
# for each settlement assignment:
#_GPattern	month	area	age
1	1	1	0	#_recr_dist_pattern1
#
#_Cond 0 # N_movement_definitions goes here if N_areas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
4 #_Nblock_Patterns
1 16 9 5 #_blocks_per_pattern
#_begin and end years of blocks
2004 2004
2007 2007 2008 2008 2009 2009 2010 2010 2011 2011 2012 2012 2013 2013 2014 2014 2015 2015 2016 2016 2017 2017 2018 2018 2019 2019 2021 2021 2022 2022 2023 2023
2006 2006 2007 2007 2008 2008 2009 2009 2010 2010 2011 2011 2012 2012 2013 2013 2014 2019
2015 2019 2020 2020 2021 2021 2022 2022 2023 2023
#
# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement
#
6 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate;_5=Maunder_M;_6=Age-range_Lorenzen
0 #_minimum age for Age-range Lorenzen M;
8 #_maximum age for Age-range Lorenzen M;
#_ #_later read 1P per Sex x G Morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr;5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0 #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
5 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
0 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env_var&link	dev_link	dev_minyr	dev_maxyr	dev_PH	Block	Block_Fxn
0.199	   0.936	       0.7	-0.393	0.31	3	  2	0	0	0	0	0	0	0	#_NatM_p_1_Fem_GP_1  
    3	      15	        10	     0	  99	0	 -3	0	0	0	0	0	0	0	#_L_at_Amin_Fem_GP_1 
   20	      30	        25	     0	  99	0	 -3	0	0	0	0	0	0	0	#_L_at_Amax_Fem_GP_1 
 0.05	    0.99	       0.4	     0	  99	0	 -3	0	0	0	0	0	0	0	#_VonBert_K_Fem_GP_1 
 0.05	     0.5	      0.14	     0	  99	0	 -3	0	0	0	0	0	0	0	#_CV_young_Fem_GP_1  
 0.01	     0.1	      0.05	     0	  99	0	 -3	0	0	0	0	0	0	0	#_CV_old_Fem_GP_1    
   -3	       3	7.5242e-06	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Wtlen_1_Fem_GP_1   
   -3	       5	    3.2332	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Wtlen_2_Fem_GP_1   
    9	      19	     15.44	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Mat50%_Fem_GP_1    
  -20	       3	  -0.89252	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Mat_slope_Fem_GP_1 
    0	      10	         1	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Eggs_alpha_Fem_GP_1
   -1	       5	         0	     0	  99	0	 -3	0	0	0	0	0	0	0	#_Eggs_beta_Fem_GP_1 
  0.1	      10	         1	     1	   1	6	 -1	0	0	0	0	0	0	0	#_CohortGrowDev      
1e-06	0.999999	       0.5	   0.5	 0.5	0	-99	0	0	0	0	0	0	0	#_FracFemale_GP_1    
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; 2=Ricker; 3=std_B-H; 4=SCAA;5=Hockey; 6=B-H_flattop; 7=survival_3Parm;8=Shepard_3Parm
0 # 0/1 to use steepness in initial equ recruitment calculation
0 # future feature: 0/1 to make realized sigmaR a function of SR curvature
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn # parm_name
  3	25	14.2	0	99	0	 1	0	0	0	0	0	0	0	#_SR_LN(R0)  
0.2	 1	 0.6	0	99	0	-5	0	0	0	0	0	0	0	#_SR_BH_steep
  0	 2	 1.2	0	99	0	-3	0	0	0	0	0	0	0	#_SR_sigmaR  
-15	15	   0	0	99	0	-1	0	0	0	0	0	1	2	#_SR_regime  
  0	 0	   0	0	99	0	-3	0	0	0	0	0	0	0	#_SR_autocorr
# timevary SR parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-15	15	2.01545	0	99	0	4	#_SR_regime_BLK1repl_2004
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
2005 # first year of main recr_devs; early devs can preceed this era
2022 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase
1 # (0/1) to read 13 advanced options
-6 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
2 #_recdev_early_phase
-1 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
1 #_lambda for Fcast_recr_like occurring before endyr+1
1992.5 #_last_yr_nobias_adj_in_MPD; begin of ramp
2002.9 #_first_yr_fullbias_adj_in_MPD; begin of plateau
2020.8 #_last_yr_fullbias_adj_in_MPD
2022.9 #_end_yr_for_ramp_in_MPD (can be in forecast to shape ramp, but SS sets bias_adj to 0.0 for fcast yrs)
0.9249 #_max_bias_adj_in_MPD (-1 to override ramp and set biasadj=1.0 for all estimated recdevs)
0 #_period of cycles in recruitment (N parms read below)
-5 #min rec_dev
5 #max rec_dev
2 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
#_Year	recdev
2023	-0.447411	#_1
2024	-0.447411	#_2
#
#Fishing Mortality info
0.1 # F ballpark
-2006 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
10 # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0	3	2.25258	0	99	0	1	#_InitF_seas_1_flt_1MexCal_S1
#
#_Q_setup for fleets with cpue or survey data
#_fleet	link	link_info	extra_se	biasadj	float  #  fleetname
    4	1	0	0	0	0	#_AT_Survey 
-9999	0	0	0	0	0	#_terminator
#_Q_parms(if_any);Qunits_are_ln(q)
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE	env-var	use_dev	dev_mnyr	dev_mxyr	dev_PH	Block	Blk_Fxn  #  parm_name
-3	3	0	0	99	0	-2	0	0	0	0	0	4	2	#_LnQ_base_AT_Survey(4)
# timevary Q parameters
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
-3	3	-0.311	0	99	0	-1	#_LnQ_base_AT_Survey(4)_BLK4repl_2015
-3	3	 -0.53	0	99	0	-1	#_LnQ_base_AT_Survey(4)_BLK4repl_2020
-3	3	-0.311	0	99	0	-1	#_LnQ_base_AT_Survey(4)_BLK4repl_2021
-3	3	     0	0	99	0	-1	#_LnQ_base_AT_Survey(4)_BLK4repl_2022
-3	3	     0	0	99	0	-1	#_LnQ_base_AT_Survey(4)_BLK4repl_2023
# info on dev vectors created for Q parms are reported with other devs after tag parameter section
#
#_size_selex_patterns
#_Pattern	Discard	Male	Special
0	0	0	0	#_1 MexCal_S1
0	0	0	0	#_2 MexCal_S2
0	0	0	0	#_3 PNW      
0	0	0	0	#_4 AT_Survey
#
#_age_selex_patterns
#_Pattern	Discard	Male	Special
17	0	0	8	#_1 MexCal_S1
17	0	0	8	#_2 MexCal_S2
12	0	0	0	#_3 PNW      
17	0	0	1	#_4 AT_Survey
#
#_SizeSelex
#_No size_selex_parm
#_AgeSelex
-7	 9	     1.23	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_1_MexCal_S1(1)
-7	 9	   3.7717	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_2_MexCal_S1(1)
-7	 9	 0.809132	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_3_MexCal_S1(1)
-7	 9	 -1.29578	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_4_MexCal_S1(1)
-7	 9	-0.223181	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_5_MexCal_S1(1)
-7	 9	  -1.1357	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_6_MexCal_S1(1)
-7	 9	 -0.16435	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_7_MexCal_S1(1)
-7	 9	-0.694264	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_8_MexCal_S1(1)
-7	 9	0.0568794	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_9_MexCal_S1(1)
-7	 9	  1.99999	-1	99	0	-3	0	0	0	0	0	0	0	#_AgeSel_P_1_MexCal_S2(2)
-7	 9	 0.670819	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_2_MexCal_S2(2)
-7	 9	-0.967526	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_3_MexCal_S2(2)
-7	 9	 -0.42972	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_4_MexCal_S2(2)
-7	 9	-0.629359	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_5_MexCal_S2(2)
-7	 9	 0.472316	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_6_MexCal_S2(2)
-7	 9	-0.224747	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_7_MexCal_S2(2)
-7	 9	 0.299453	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_8_MexCal_S2(2)
-7	 9	 -0.65828	-1	99	0	 3	0	0	0	0	0	0	0	#_AgeSel_P_9_MexCal_S2(2)
 0	10	  3.40597	 0	99	0	 4	0	0	0	0	0	3	2	#_AgeSel_P_1_PNW(3)      
-5	15	  1.36441	 0	99	0	 4	0	0	0	0	0	0	0	#_AgeSel_P_2_PNW(3)      
 0	 9	        0	-1	99	0	-3	0	0	0	0	0	0	0	#_AgeSel_P_1_AT_Survey(4)
 0	 9	      0.1	-1	99	0	 4	0	0	0	0	0	2	2	#_AgeSel_P_2_AT_Survey(4)
# timevary selex parameters 
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2006      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2007      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2008      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2009      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2010      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2011      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2012      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2013      
0	10	3.40597	 0	99	0	4	#_AgeSel_P_1_PNW(3)_BLK3repl_2014      
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2007
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2008
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2009
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2010
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2011
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2012
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2013
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2014
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2015
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2016
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2017
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2018
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2019
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2021
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2022
0	 9	    0.1	-1	99	0	4	#_AgeSel_P_2_AT_Survey(4)_BLK2repl_2023
# info on dev vectors created for selex parms are reported with other devs after tag parameter section
#
1 #  use 2D_AR1 selectivity(0/1):  experimental feature
#_specifications for 2D_AR1 and associated parameters
#_fleet	ymin	ymax	amin	amax	sig_amax	use_rho	l1/a2	devphase	before_range	after_range
1	2006	2014	0	3	0	1	2	3	0	0	#_1
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0.001	  2	1	0	99	0	-3	#_sigma_sel_fleet1_age0
 -0.8	0.8	0	0	99	0	-3	#_rho_year_fleet1      
 -0.8	0.8	0	0	99	0	-3	#_rho_age_fleet1       
#_fleet	ymin	ymax	amin	amax	sig_amax	use_rho	l1/a2	devphase	before_range	after_range
2	2006	2014	0	4	0	1	2	3	0	0	#_2
#_LO	HI	INIT	PRIOR	PR_SD	PR_type	PHASE
0.001	  2	1	0	99	9	-3	#_sigma_sel_fleet2_age0
 -0.8	0.8	0	0	99	0	-3	#_rho_year_fleet2      
 -0.8	0.8	0	0	99	0	-3	#_rho_age_fleet2       
-9999 1 1 1 1 1 1 1 1 1 1 # Terminator 
# Tag loss and Tag reporting parameters go next
0 # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# Input variance adjustments factors: 
#_Data_type Fleet Value
-9999 1 0 # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 16 changes to default Lambdas (default value is 1.0)
#_like_comp	fleet	phase	value	sizefreq_method
    1	4	1	1	1	#_Surv_AT_Survey_Phz1                    
    4	4	1	0	1	#_length_AT_Survey_sizefreq_method_1_Phz1
    5	1	1	1	1	#_age_MexCal_S1_Phz1                     
    5	2	1	1	1	#_age_MexCal_S2_Phz1                     
    5	3	1	1	1	#_age_PNW_Phz1                           
    5	4	1	1	1	#_age_AT_Survey_Phz1                     
    8	1	1	1	1	#_catch_MexCal_S1_Phz1                   
    8	2	1	1	1	#_catch_MexCal_S2_Phz1                   
    8	3	1	1	1	#_catch_PNW_Phz1                         
    9	1	1	0	1	#_init_equ_catch_MexCal_S1_Phz1          
    9	2	1	0	1	#_init_equ_catch_MexCal_S2_Phz1          
    9	3	1	0	1	#_init_equ_catch_PNW_Phz1                
   18	1	1	0	1	#_Regime-shift_MexCal_S1_Phz1            
   18	2	1	0	1	#_Regime-shift_MexCal_S2_Phz1            
   18	3	1	0	1	#_Regime-shift_PNW_Phz1                  
   18	4	1	0	1	#_Regime-shift_AT_Survey_Phz1            
-9999	0	0	0	0	#_terminator                             
#
2 # 0/1 read specs for more stddev reporting
0 2 0 0 # selex_fleet, 1=len/2=age/3=both, year, N selex bins
0 0       # 0 or Growth pattern, N growth ages
0 0 0    # 0 or NatAge_area(-1 for sum), NatAge_yr, N Natages
0 0 0 1     # Mortality, Dyn B0 (>3.30.16), SmryBio (>3.30.16) 
#
999
