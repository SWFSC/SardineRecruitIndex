#V3.30.13.00-trans;_2019_03_09;_Stock_Synthesis_by_Richard_Methot_(NOAA)_using_ADMB_12.0
#Stock Synthesis (SS) is a work of the U.S. Government and is not subject to copyright protection in the United States.
#Foreign copyrights may apply. See copyright.txt for more information.
#_user_support_available_at:NMFS.Stock.Synthesis@noaa.gov
#_user_info_available_at:https://vlab.ncep.noaa.gov/group/stock-synthesis
#_data_and_control_files: ALT_19.dat // ALT_19.ctl
1  # 0 means do not read wtatage.ss; 1 means read and use wtatage.ss and also read and use growth parameters
1  #_N_Growth_Patterns
1 #_N_platoons_Within_GrowthPattern 
#_Cond 1 #_Morph_between/within_stdev_ratio (no read if N_morphs=1)
#_Cond  1 #vector_Morphdist_(-1_in_first_val_gives_normal_approx)
#
4 # recr_dist_method for parameters:  2=main effects for GP, Settle timing, Area; 3=each Settle entity; 4=none, only when N_GP*Nsettle*pop==1
1 # not yet implemented; Future usage: Spawner-Recruitment: 1=global; 2=by area
1 #  number of recruitment settlement assignments 
0 # unused option
#GPattern month  area  age (for each settlement assignment)
 1 1 1 0
#
#_Cond 0 # N_movement_definitions goes here if Nareas > 1
#_Cond 1.0 # first age that moves (real age at begin of season, not integer) also cond on do_migration>0
#_Cond 1 1 1 2 4 10 # example move definition for seas=1, morph=1, source=1 dest=2, age1=4, age2=10
#
4 #_Nblock_Patterns
 1 16 9 5 #_blocks_per_pattern 
# begin and end years of blocks
 2004 2004

 #Survey blocks; No survey in 2006 and 2005 is the base selex pattern
 2007 2007
 2008 2008
 2009 2009
 2010 2010
 2011 2011
 2012 2012
 2013 2013
 2014 2014
 2015 2015
 2016 2016 
 2017 2017
 2018 2018
 2019 2019
 2021 2021
 2022 2022
 2023 2023

 #Fishery blocks; 2005 is base selex pattern
 2006 2006 
 2007 2007
 2008 2008 
 2009 2009
 2010 2010
 2011 2011
 2012 2012
 2013 2013
 2014 2019
#
#Q prior block
2015 2019
2020 2020
2021 2021
2022 2022
2023 2023

# controls for all timevary parameters 
1 #_env/block/dev_adjust_method for all time-vary parms (1=warn relative to base parm bounds; 3=no bound check)
#
# AUTOGEN
 1 1 1 1 1 # autogen: 1st element for biology, 2nd for SR, 3rd for Q, 4th reserved, 5th for selex
# where: 0 = autogen all time-varying parms; 1 = read each time-varying parm line; 2 = read then autogen if parm min==-12345
#
#_Available timevary codes
#_Block types: 0: P_block=P_base*exp(TVP); 1: P_block=P_base+TVP; 2: P_block=TVP; 3: P_block=P_block(-1) + TVP
#_Block_trends: -1: trend bounded by base parm min-max and parms in transformed units (beware); -2: endtrend and infl_year direct values; -3: end and infl as fraction of base range
#_EnvLinks:  1: P(y)=P_base*exp(TVP*env(y));  2: P(y)=P_base+TVP*env(y);  3: null;  4: P(y)=2.0/(1.0+exp(-TVP1*env(y) - TVP2))
#_DevLinks:  1: P(y)*=exp(dev(y)*dev_se;  2: P(y)+=dev(y)*dev_se;  3: random walk;  4: zero-reverting random walk with rho;  21-24 keep last dev for rest of years
#
#
#
# setup for M, growth, maturity, fecundity, recruitment distibution, movement 
#
6 #_natM_type:_0=1Parm; 1=N_breakpoints;_2=Lorenzen;_3=agespecific;_4=agespec_withseasinterpolate
  #_no additional input for selected M option; read 1P per morph
#
0 #_minimum age for Lorenzen
8 #_maximum age for Lorenzen; read 1P per morph
1 # GrowthModel: 1=vonBert with L1&L2; 2=Richards with L1&L2; 3=age_specific_K_incr; 4=age_specific_K_decr; 5=age_specific_K_each; 6=NA; 7=NA; 8=growth cessation
0.5 #_Age(post-settlement)_for_L1;linear growth below this
999 #_Growth_Age_for_L2 (999 to use as Linf)
-999 #_exponential decay for growth above maxage (value should approx initial Z; -999 replicates 3.24; -998 to not allow growth above maxage)
0  #_placeholder for future growth feature
#
0 #_SD_add_to_LAA (set to 0.1 for SS2 V1.x compatibility)
0 #_CV_Growth_Pattern:  0 CV=f(LAA); 1 CV=F(A); 2 SD=F(LAA); 3 SD=F(A); 4 logSD=F(A)
#
5 #_maturity_option:  1=length logistic; 2=age logistic; 3=read age-maturity matrix by growth_pattern; 4=read age-fecundity; 5=disabled; 6=read length-maturity
#_Age_Fecundity by growth pattern from wt-at-age.ss now invoked by read bodywt flag
0 #_First_Mature_Age
1 #_fecundity option:(1)eggs=Wt*(a+b*Wt);(2)eggs=a*L^b;(3)eggs=a*Wt^b; (4)eggs=a+b*L; (5)eggs=a+b*W
0 #_hermaphroditism option:  0=none; 1=female-to-male age-specific fxn; -1=male-to-female age-specific fxn
1 #_parameter_offset_approach (1=none, 2= M, G, CV_G as offset from female-GP1, 3=like SS2 V1.x)
#
#_growth_parms
#_ LO HI INIT PRIOR PR_SD PR_type PHASE env_var&link dev_link dev_minyr dev_maxyr dev_PH Block Block_Fxn
# Sex: 1  BioPattern: 1  NatMort
0.199 0.936 0.7 -0.3930 0.31 3 2 0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
# 0.3 0.8 0.6 0 99 0 3 0 0 0 0 0 0 0 # NatM_p_1_Fem_GP_1
# Sex: 1  BioPattern: 1  Growth 
 3 15 10 0 99 0 -3 0 0 0 0 0 0 0 # L_at_Amin_Fem_GP_1
 20 30 25 0 99 0 -3 0 0 0 0 0 0 0 # L_at_Amax_Fem_GP_1
 0.05 0.99 0.4 0 99 0 -3 0 0 0 0 0 0 0 # VonBert_K_Fem_GP_1
 0.05 0.5 0.14 0 99 0 -3 0 0 0 0 0 0 0 # CV_young_Fem_GP_1
 0.01 0.1 0.05 0 99 0 -3 0 0 0 0 0 0 0 # CV_old_Fem_GP_1
# Sex: 1  BioPattern: 1  WtLen
 -3 3 7.5242e-006 0 99 0 -3 0 0 0 0 0 0 0 # Wtlen_1_Fem
 -3 5 3.2332 0 99 0 -3 0 0 0 0 0 0 0 # Wtlen_2_Fem
# Sex: 1  BioPattern: 1  Maturity&Fecundity
 9 19 15.44 0 99 0 -3 0 0 0 0 0 0 0 # Mat50%_Fem
 -20 3 -0.89252 0 99 0 -3 0 0 0 0 0 0 0 # Mat_slope_Fem
 0 10 1 0 99 0 -3 0 0 0 0 0 0 0 # Eggs/kg_inter_Fem
 -1 5 0 0 99 0 -3 0 0 0 0 0 0 0 # Eggs/kg_slope_wt_Fem
# Hermaphroditism
#  Recruitment Distribution  
 #-4 4 0 0 99 0 -3 0 0 0 0 0 0 0 # RecrDist_GP_1
 #-4 4 1 0 99 0 -3 0 0 0 0 0 0 0 # RecrDist_Area_1
 #-4 4 1 0 99 0 -3 0 0 0 0 0 0 0 # RecrDist_timing_1
#  Cohort growth dev base
 0.1 10 1 1 1 6 -1 0 0 0 0 0 0 0 # CohortGrowDev
#  Movement
#  Age Error from parameters
#  catch multiplier
#  fraction female, by GP
 0.000001 0.999999 0.5 0.5  0.5 0 -99 0 0 0 0 0 0 0 # FracFemale_GP_1
#
#_no timevary MG parameters
#
#_seasonal_effects_on_biology_parms
 0 0 0 0 0 0 0 0 0 0 #_femwtlen1,femwtlen2,mat1,mat2,fec1,fec2,Malewtlen1,malewtlen2,L1,K
#_ LO HI INIT PRIOR PR_SD PR_type PHASE
#_Cond -2 2 0 0 -1 99 -2 #_placeholder when no seasonal MG parameters
#
3 #_Spawner-Recruitment; Options: 2=Ricker; 3=std_B-H; 4=SCAA; 5=Hockey; 6=B-H_flattop; 7=survival_3Parm; 8=Shepherd_3Parm; 9=RickerPower_3parm
0  # 0/1 to use steepness in initial equ recruitment calculation
0  #  future feature:  0/1 to make realized sigmaR a function of SR curvature
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn #  parm_name
             3            25          14.2             0            99             0          1          0          0          0          0          0          0          0 # SR_LN(R0)
           0.2             1           0.6               0            99             0         -5          0          0          0          0          0          0          0 # SR_BH_steep
             0             2           1.2             0            99             0         -3          0          0          0          0          0          0          0 # SR_sigmaR
           -15            15             0             0            99             0         -1          0          0          0          0          0          1          2 # SR_regime
             0             0             0             0            99             0         -3          0          0          0          0          0          0          0 # SR_autocorr
#Next are short parm lines for timevary 
#_LO HI INIT      PRIOR PR_SD PR_type PHASE
-15 15  2.01545       0   99       0    4 # SR_regime_BLK1add_2004_2004
1 #do_recdev:  0=none; 1=devvector (R=F(SSB)+dev); 2=deviations (R=F(SSB)+dev); 3=deviations (R=R0*dev; dev2=R-f(SSB)); 4=like 3 with sum(dev2) adding penalty
2005 # first year of main recr_devs; early devs can preceed this era
2023 # last year of main recr_devs; forecast devs start in following year
1 #_recdev phase 
1 # (0/1) to read 13 advanced options
 -6 #_recdev_early_start (0=none; neg value makes relative to recdev_start)
 2 #_recdev_early_phase
 0 #_forecast_recruitment phase (incl. late recr) (0 value resets to maxphase+1)
 1 #_lambda for Fcast_recr_like occurring before endyr+1
1992.5   #_last_early_yr_nobias_adj_in_MPD 
2002.9   #_first_yr_fullbias_adj_in_MPD 
2020.8   #_last_yr_fullbias_adj_in_MPD 
2023.1   #_first_recent_yr_nobias_adj_in_MPD 
0.9249  #_max_bias_adj_in_MPD (1.0 to mimic pre-2009 models)
 0 #_period of cycles in recruitment (N parms read below)
 -5 #min rec_dev
 5 #max rec_dev
 0 #_read_recdevs
#_end of advanced SR options
#
#_placeholder for full parameter lines for recruitment cycles
# read specified recr devs
#_Yr Input_value
#
# all recruitment deviations
#  1999E 2000E 2001E 2002E 2003E 2004E 2005R 2006R 2007R 2008R 2009R 2010R 2011R 2012R 2013R 2014R 2015R 2016R 2017R 2018F 2019F
#  -0.268996 -0.317947 -0.407416 0.0690791 0.892756 0.15518 -0.192824 -0.475354 -0.955199 -0.180308 0.206281 -1.63497 -2.4461 -1.71318 0.281394 0.547216 0.534335 1.04716 4.98155 0 0
# implementation error by year in forecast:  0
#
#Fishing Mortality info 
0.1 # F ballpark
-2006 # F ballpark year (neg value to disable)
3 # F_Method:  1=Pope; 2=instan. F; 3=hybrid (hybrid is recommended)
4 # max F or harvest rate, depends on F_Method
# no additional F input needed for Fmethod 1
# if Fmethod=2; read overall start F value; overall phase; N detailed inputs to read
# if Fmethod=3; read N iterations for tuning for Fmethod 3
10  # N iterations for tuning F in hybrid method (recommend 3 to 7)
#
#_initial_F_parms; count = 1
#_ LO HI INIT PRIOR PR_SD  PR_type  PHASE
 0 3 2.25258 0 99 0 1 # InitF_seas_1_flt_1MexCal_S1
#2019 2074
# F rates by fleet
# Yr:  2005 2005 2006 2006 2007 2007 2008 2008 2009 2009 2010 2010 2011 2011 2012 2012 2013 2013 2014 2014 2015 2015 2016 2016 2017 2017 2018 2018 2019 2019
# seas:  1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2
# MexCal_S1 0.0438255 0 0.0660767 0 0.176544 0 0.187782 0 0.135454 0 0.13695 0 0.216923 0 0.0187663 0 0.0356888 0 0.175152 0 0.00049481 0 0.0142116 0 0.00991665 0 0.000951018 0 5.05195e-05 0
# MexCal_S2 0 0.0821925 0 0.119609 0 0.224956 0 0.258034 0 0.275264 0 0.168161 0 0.199346 0 0.356185 0 0.335375 0 0.053009 0 0.00880949 0 0.488743 0 0.487494 0 0.0176186 0 0.00958493
# PNW 1.38581 0.00373427 0.44885 0 0.367948 0 0.256983 0 0.339426 0.0145672 0.582787 1.3703e-06 0.571456 0.12506 2.17521 0.070177 1.7188 0.0480624 0.59509 0.135034 0.00409881 0.000115873 0.0074476 1.15541e-05 0.000104843 0.000262809 0.000606653 0.000212787 0.000314888 0.000108745
#
#_Q_setup for fleets with cpue or survey data
#_1:  fleet number
#_2:  link type: (1=simple q, 1 parm; 2=mirror simple q, 1 mirrored parm; 3=q and power, 2 parm; 4=mirror with offset, 2 parm)
#_3:  extra input for link, i.e. mirror fleet# or dev index number
#_4:  0/1 to select extra sd parameter
#_5:  0/1 for biasadj or not
#_6:  0/1 to float
#_   fleet      link link_info  extra_se   biasadj     float  #  fleetname
         4         1         0         0         0         0  #  AT_Survey
         5         1         0         1         0         0  #  BEUTI
         #5         1         0         0         0         0  #  Lisa_Marie
-9999 0 0 0 0 0
#
#_Q_parms(if_any);Qunits_are_ln(q)
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn   #parm_name
            -3             3             0             0            99             0          -2          0          0          0          0          0          4          2  #LnQ_base_ATM_Summer(4)            
            -3             3             -0.05             0            99             0          -2          0          0          0          0          0          0          0  #LnQ_dfaT1(5) 
            0.001 1.3  0.949309   0   99  0    2  0   0   0   0   0   0   0   #_Q_extraSD_dfaT1(5) # extra_se estimate, borrowed from sablefish assessment           
            #-3             3            -0.4815         0            99             0          -2          0          0          0          0          0          0          0  #LnQ_base_Lisa_Marie            
            #-3             3             0             0            .1             6          -2          0          0          0          0          0          4          2  #  LnQ_base_ATM_Summer(4)
            #-3             3      0.157183             0            99             0          4          0          0          0          0          0          0          0  #  LnQ_base_AT_Survey(4)
#_no timevary Q parameters
#
-3 3 -0.311 0 99  0 -1 # Q block
-3 3 -0.530 0 99  0 -1 # Q block spring 2020
-3 3 -0.311 0 99  0 -1 # Q block 2021
-3 3 0      0 99 0 -1 # Q block 2022
-3 3 0       0 99 0 -1 # Q block 2023
#_size_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for all sizes
#Pattern:_1; parm=2; logistic; with 95% width specification
#Pattern:_5; parm=2; mirror another size selex; PARMS pick the min-max bin to mirror
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_6; parm=2+special; non-parm len selex
#Pattern:_43; parm=2+special+2;  like 6, with 2 additional param for scaling (average over bin range)
#Pattern:_8; parm=8; New doublelogistic with smooth transitions and constant above Linf option
#Pattern:_9; parm=6; simple 4-parm double logistic with starting length; parm 5 is first length; parm 6=1 does desc as offset
#Pattern:_21; parm=2+special; non-parm len selex, read as pairs of size, then selex
#Pattern:_22; parm=4; double_normal as in CASAL
#Pattern:_23; parm=6; double_normal where final value is directly equal to sp(6) so can be >1.0
#Pattern:_24; parm=6; double_normal with sel(minL) and sel(maxL), using joiners
#Pattern:_25; parm=3; exponential-logistic in size
#Pattern:_27; parm=3+special; cubic spline 
#Pattern:_42; parm=2+special+3; // like 27, with 2 additional param for scaling (average over bin range)
#_discard_options:_0=none;_1=define_retention;_2=retention&mortality;_3=all_discarded_dead;_4=define_dome-shaped_retention
#_Pattern Discard Male Special
 0 0 0 0 # 1 MexCal_S1
 0 0 0 0 # 2 MexCal_S2
 0 0 0 0 # 3 PNW
 0 0 0 0 # 4 AT_Survey
 0 0 0 0 # 5 BEUTI
# 0 0 0 0 # 5 AT_Survey
#
#_age_selex_patterns
#Pattern:_0; parm=0; selex=1.0 for ages 0 to maxage
#Pattern:_10; parm=0; selex=1.0 for ages 1 to maxage
#Pattern:_11; parm=2; selex=1.0  for specified min-max age
#Pattern:_12; parm=2; age logistic
#Pattern:_13; parm=8; age double logistic
#Pattern:_14; parm=nages+1; age empirical
#Pattern:_15; parm=0; mirror another age or length selex
#Pattern:_16; parm=2; Coleraine - Gaussian
#Pattern:_17; parm=nages+1; empirical as random walk  N parameters to read can be overridden by setting special to non-zero
#Pattern: _41; parm=2+nages+1; // like 17, with 2 additional param for scaling (average over bin range)
#Pattern:_18; parm=8; double logistic - smooth transition
#Pattern:_19; parm=6; simple 4-parm double logistic with starting age
#Pattern:_20; parm=6; double_normal,using joiners
#Pattern:_26; parm=3; exponential-logistic in age
#Pattern:_27; parm=3+special; cubic spline in age
#Pattern:_42; parm=2+special+3; // cubic spline; with 2 additional param for scaling (average over bin range)
#_Pattern Discard Male Special
 17 0 0 8 # 1 MexCal_S1
 17 0 0 8 # 2 MexCal_S2
 12 0 0 0 # 3 PNW
 17 0 0 1 # 4 AT_Survey
 11  0 0 0 # 5 BEUTI index on recruits only
# 17 0 0 1 # 5 Lisa Marie
#
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type      PHASE    env-var    use_dev   dev_mnyr   dev_mxyr     dev_PH      Block    Blk_Fxn  #  parm_name
# 1   MexCal_S1 LenSelex
# 2   MexCal_S2 LenSelex
# 3   PNW LenSelex
# 4   AT_Survey LenSelex
# 1   MexCal_S1 AgeSelex
            -7             9          1.23            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P1_MexCal_S1(1)
            -7             9        3.7717            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P2_MexCal_S1(1)
            -7             9      0.809132            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P3_MexCal_S1(1)
            -7             9      -1.29578            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P4_MexCal_S1(1)
            -7             9     -0.223181            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P5_MexCal_S1(1)
            -7             9       -1.1357            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P6_MexCal_S1(1)
            -7             9      -0.16435            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P7_MexCal_S1(1)
            -7             9     -0.694264            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P8_MexCal_S1(1)
            -7             9     0.0568794            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P9_MexCal_S1(1)

# 2   MexCal_S2 AgeSelex
            -7             9       1.99999            -1            99             0         -3          0          0          0          0          0          0          0  #  AgeSel_P1_MexCal_S2(2)
            -7             9      0.670819            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P2_MexCal_S2(2)
            -7             9     -0.967526            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P3_MexCal_S2(2)
            -7             9      -0.42972            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P4_MexCal_S2(2)
            -7             9     -0.629359            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P5_MexCal_S2(2)
            -7             9      0.472316            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P6_MexCal_S2(2)
            -7             9     -0.224747            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P7_MexCal_S2(2)
            -7             9      0.299453            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P8_MexCal_S2(2)
            -7             9      -0.65828            -1            99             0          3          0          0          0          0          0          0          0  #  AgeSel_P9_MexCal_S2(2)

# 3   PNW AgeSelex
             0            10       3.40597            0            99             0          4          0          0          0          0          0          3          2  #  AgeSel_P1_PNW(3)
            -5            15       1.36441            0            99             0          4          0          0          0          0          0          0          0  #  AgeSel_P2_PNW(3)
# 4   AT_Survey AgeSelex
            0             9             0            -1            99             0         -3          0          0          0          0          0          0          0  #AgeSel_P1_AT_Survey(4)
            0             9             .1           -1            99             0          4          0          0          0          0          0          2          2  #AgeSel_P2_AT_Survey(4)

# 5 dfaT1
            0             9             0            -1            99             0         -3          0          0          0          0          0          0          0  #AgeSel_P1_dfaT1_Survey(5)
            0             9             0            -1            99             0         -3          0          0          0          0          0          0          0  #AgeSel_P2_dfaT1_Survey(5)

# timevary selex parameters 
#_          LO            HI          INIT         PRIOR         PR_SD       PR_type    PHASE  #  parm_name
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2006
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2007
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2008
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2009
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2010
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2011
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2012
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2013
            #-7             9        3.7717            -1            99             0      3  # AgeSel_P2_MexCal_S1(1)_BLK3repl_2014
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2006
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2007
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2008
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2009
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2010
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2011
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2012
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2013
            #-7             9      0.809132            -1            99             0      3  # AgeSel_P3_MexCal_S1(1)_BLK3repl_2014
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2006
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2007
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2008
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2009
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2010
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2011
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2012
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2013
            #-7             9      -1.29578            -1            99             0      3  # AgeSel_P4_MexCal_S1(1)_BLK3repl_2014
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2006
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2007
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2008
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2009
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2010
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2011
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2012
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2013
            #-7             9     -0.223181            -1            99             0      3  # AgeSel_P5_MexCal_S1(1)_BLK3repl_2014
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2006
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2007
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2008
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2009
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2010
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2012
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2013
            #-7             9       -1.1357            -1            99             0      3  # AgeSel_P6_MexCal_S1(1)_BLK3repl_2014
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2006
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2007
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2008
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2009
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2010
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2011
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2012
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2013
            #-7             9      0.670819            -1            99             0      3  # AgeSel_P2_MexCal_S2(2)_BLK3repl_2014
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2006
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2007
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2008
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2009
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2010
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2011
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2012
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2013
            #-7             9     -0.967526            -1            99             0      3  # AgeSel_P3_MexCal_S2(2)_BLK3repl_2014
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2006
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2007
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2008
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2009
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2010
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2011
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2012
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2013
            #-7             9      -0.42972            -1            99             0      3  # AgeSel_P4_MexCal_S2(2)_BLK3repl_2014
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2006
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2007
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2008
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2009
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2010
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2011
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2012
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2013
             0            10       3.40597             0            99             0      4  # Age_inflection_PNW(3)_BLK3repl_2014
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2007
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2008
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2009
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2010
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2011
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2012
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2013
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2014
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2015
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2016
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2017
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2018
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2019
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2021
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2022
             0             9           0.1            -1            99             0      4  # AgeSel_P2_AT_Survey(4)_BLK2repl_2023

1   #  use 2D_AR1 selectivity(0/1):  experimental feature
#_specifications for 2D_AR1 and associated parameters
#_specs:  fleet, ymin, ymax, amin, amax, sigma_amax, use_rho, len1/age2, devphase, before_range, after_range
#_sigma_amax>amin means create sigma parm for each bin from min to sigma_amax; sigma_amax<0 means just one sigma parm is read and used for all bins
 1 2006 2014 0 3 0 1 2 3 0 0  #  2d_AR specs for fleet: FISHERY AGE
 0.001 2  1 0 99 0 -3  # sigma_sel for fleet:_1; AGE_0
 -0.8 0.8 0 0 99 0 -3  # rho_year for fleet:_1
 -0.8 0.8 0 0 99 0 -3  # rho_AGE for fleet:_1

 2 2006 2014 0 4 0 1 2 3 0 0  #  2d_AR specs for fleet: FISHERY AGE
 0.001 2 1 0 99 9 -3  # sigma_sel for fleet:_1; AGE_0
 -0.8 0.8 0 0 99 0 -3  # rho_year for fleet:_1
 -0.8 0.8 0 0 99 0 -3  # rho_AGE for fleet:_1
-9999  0 0 0 0 0 0 0 0 0 0 # terminator
# Tag loss and Tag reporting parameters go next
0  # TG_custom:  0=no read; 1=read if tags exist
#_Cond -6 6 1 1 2 0.01 -4 0 0 0 0 0 0 0  #_placeholder if no parameters
#
# no timevary parameters
#
#
# Input variance adjustments factors: 
 #_1=add_to_survey_CV
 #_2=add_to_discard_stddev
 #_3=add_to_bodywt_CV
 #_4=mult_by_lencomp_N
 #_5=mult_by_agecomp_N
 #_6=mult_by_size-at-age_N
 #_7=mult_by_generalized_sizecomp
#_Factor  Fleet  Value
 -9999   1    0  # terminator
#
1 #_maxlambdaphase
1 #_sd_offset; must be 1 if any growthCV, sigmaR, or survey extraSD is an estimated parameter
# read 9 changes to default Lambdas (default value is 1.0)
# Like_comp codes:  1=surv; 2=disc; 3=mnwt; 4=length; 5=age; 6=SizeFreq; 7=sizeage; 8=catch; 9=init_equ_catch; 
# 10=recrdev; 11=parm_prior; 12=parm_dev; 13=CrashPen; 14=Morphcomp; 15=Tag-comp; 16=Tag-negbin; 17=F_ballpark; 18=initEQregime
#like_comp fleet  phase  value  sizefreq_method
 1 4 1 1 1
 1 5 1 1 1
 4 4 1 0 1
 5 1 1 1 1
 5 2 1 1 1
 5 3 1 1 1
 5 4 1 1 1
 #5 5 1 1 1
 8 1 1 1 1
 8 2 1 1 1
 8 3 1 1 1
 9 1 1 0 1
 9 2 1 0 1
 9 3 1 0 1
 #12 1 1 1 1
 #12 2 1 1 1
 18 1 1 0 1
 18 2 1 0 1
 18 3 1 0 1
 18 4 1 0 1
-9999  1  1  1  1  #  terminator
#
# lambdas (for info only; columns are phases)
#  0 #_CPUE/survey:_1
#  0 #_CPUE/survey:_2
#  0 #_CPUE/survey:_3
#  1 #_CPUE/survey:_4
#  0 #_lencomp:_1
#  0 #_lencomp:_2
#  0 #_lencomp:_3
#  0 #_lencomp:_4
#  1 #_agecomp:_1
#  1 #_agecomp:_2
#  1 #_agecomp:_3
#  1 #_agecomp:_4
#  0 #_init_equ_catch
#  1 #_recruitments
#  1 #_parameter-priors
#  1 #_parameter-dev-vectors
#  1 #_crashPenLambda
#  0 # F_ballpark_lambda
2 # (0/1) read specs for more stddev reporting 
0 2 0 0 # Selectivity: (1) fleet, (2) 1=len/2=age/3=both, (3) year, (4) N selex bins
0 0 # Growth: (1) growth pattern, (2) growth ages
0 0 0 # Numbers-at-age: (1) area(-1 for all), (2) year, (3) N ages
0 0 # Mortality: (1) 0 to skip or growth pattern, (2) N ages for mortality; NOTE: does each sex
0 # Dyn Bzero: 0 to skip, 1 to include, or 2 to add recr
1 # SmryBio: 0 to skip, 1 to include
 # 0 0 0 0 0 0 0 0 0 # placeholder for # selex_fleet, 1=len/2=age/3=both, year, N selex bins, 0 or Growth pattern, N growth ages, 0 or NatAge_area(-1 for all), NatAge_yr, N Natages
 # placeholder for vector of selex bins to be reported
 # placeholder for vector of growth ages to be reported
 # placeholder for vector of NatAges ages to be reported
999

