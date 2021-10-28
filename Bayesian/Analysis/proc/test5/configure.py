#!/usr/bin/env python

#-> Configuration file of Bayeisan Inversion System. <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

###################################
#         Primary setting         #
###################################

#--- Case name. ---#
CASENAME = "test5"

#--- Experiment name. ---#
#EXPNAME = "independent_priorFnet"
EXPNAME = "const_emiss"

#--- Experiment type ("OSSE" or "REAL" currently) ---#
EXPTYPE = "OSSE"

#--- Optimization core. ---#
CORE = "3DVar"

#--- Start and end time definition. ---#
START = "2018-11-29_00:00:00" # Start (YYYY-mm-dd_HH:MM:SS)
END = "2018-12-03_00:00:00" # End (YYYY-mm-dd_HH:MM:SS)
DTHRS = 1 # Time interval (unit: hour)

#--- Optimization parameters definition. ---#
WINDOW = 24 # Window of assimilation.
NBACK = 24 # Backward time of STILT model. (unit: hour)

LS = 6 # Spatial length. (unit: km)
LT = 48 # Temporal length. (unit: hour)

LSBCK = 1000 # Spatial length of background. (unit: km)
LTBCK = 24 * 10 # Temporal length of background. (unit: hour)

POSTEMISSCOL = "copy" # BIS output files type ("copy" or "link")

OUTPUTTYPE = "mean" # One of "sum" or "mean"

#--- Receptors definition. ---#

from .stations.picarro import RECEPTORS_DIC
#from .stations.pis import RECEPTORS_DIC

RECEPTORS = list(RECEPTORS_DIC.keys())
RECEPTORS_ERROR = {staName: RECEPTORS_DIC[staName][3] for staName in RECEPTORS}

DRAWLIST = [
    #"movie", 
    "spatial_plots",
    "site_timeSeries",
    "pointRegion_timeSeries",
]

SPATIAL_TIMESCALELIST = [
    "daily",
    #"weekly",
]

PRTS_TIMESCALELIST = [
    #"hourly",
    "daily", 
    #"weekly",
] 


#--- Directories definition. ---#
ROOTDIR = "/home/tangwh/modeling/BIS_v6.0" # by init.sh
CASESDIR = "/home/tangwh/modeldata/BIS_cases"
#FOOTDIR = "/home/tangwh/datasets/STILT_out/PIS_DJF2019"
FOOTDIR = "/home/tangwh/datasets/STILT_out/DJF2019/stilt_out"

PICARRODIR = "/home/tangwh/datasets/SENSE/data/min/L2.0"

BENCHDIR = "/home/tangwh/public_html/BIS_show"

SRCDIR = ROOTDIR + "/Bayesian"

CASEDIR = CASESDIR + "/" + CASENAME
CASELNK = ROOTDIR + "/cases/" + CASENAME

EMISSDIR = CASEDIR + "/emiss"
SIGMADIR = CASEDIR + "/sigma"
DDIR = CASEDIR + "/D_save"
EDIR = CASEDIR + "/E_save"
HDIR = CASEDIR + "/H_save"
RDIR = CASEDIR + "/R_save"
OBSDIR = CASEDIR + "/obs"
WORKDIR = CASEDIR + "/work"
BCKDIR = CASEDIR + "/bck"
SBCKDIR = CASEDIR + "/sigmaBck"
INITCONCDIR = CASEDIR + "/initConc"
FINALCONCDIR = CASEDIR + "/finalConc"
OUTDIR = CASEDIR + "/output"

#--- Option of whether the background will be optimized. ---#
OPTBCK = True

#--- Multiprocessing parameters. ---#
NPROC = 18 # CPU cores

#--- Screen printing parameters. ---#
isShowProg = True

#--- Draw BenchPlots. ---#
isDrawBench = True


###################################
#        Advanced setting         #
###################################

#--- Core optimization settings. ---#
CoreConfig = {

    #--- Settings for core_optimize. (the base class core) ---#
    "optimize": {
        "nBack": NBACK, # Backward time of STILT model, defaulted as same as above.
        "workDir": WORKDIR, # Directory of work folder, defaulted as same as above.

        "HQ_Prefix": "HQ", # HQ files prefix.
        "HQHt_Prefix": "HQHt",
        "HQHtR_Prefix": "HQHtR",
        "INV_Prefix": "INV",
        "HxPrior_Prefix": "Hx_prior",
        "HxProc_Prefix": "Hx_proc",
        "HxBckPrior_Prefix": "HxBck_prior",
        "HxBckProc_Prefix": "HxBck_proc",
        "d_Prefix": "d",
        "dsigma2_Prefix": "dSigma2",

        #--- Option of whether background should be optimized. ---#
        "optBck": OPTBCK, 

        
    },

    #--- Settings for core_3DVar. ---#
    "3DVar": {
        "window": WINDOW, # Assimilation window (unit: hour), defaulted as same as above.
        #"iteration": "3DVarWindow", # Assimilation iteration method.
        "iteration": "3DVar_specific_hour", # Assimilation iteration method.
        "iteration_kwargs":
        {
            "hourList": [14, 15, 16, 17],
            "UTC": +8,
        }
    },

    #--- Settings for core_4DVar. ---#
    "4DVar": {
    },

}

#--- Experiments settings. ---#
ExpConfig = {

    #--- Setting for ExpBIS ---#
    "BIS": 
    {
        #--- Start & end time ---#
        "start": START, # Start time, defaulted as same as above.
        "end": END, # End time, defaulted as same as above.
        "dtHrs": DTHRS, # Time interval (unit: hour), defaulted as same as above.

        #--- CPU numbers ---#
        "nProc": NPROC, # CPU cores, defaulted as same as above.

        #--- Optimization core ---#
        "core": CORE, # Core optimization, defaulted as same as above.

        #--- BIS output file type (copy or link) ---#
        "postEmissCoL": POSTEMISSCOL,

        #--- output value type for different time scale ---#
        "outputType": OUTPUTTYPE,

        #--- Default sigle sector setting. ---#
        "sectors": ["total"], # Emission sectors definition. (using apython list format)
        "typeToSector_H": {"site" : ["total"]}, # H matrix type -> sector correspondence.
        ### format {"<H type>": ["<sector1>", "<sector2>", ...]}
        "typeToSector_E": {"site" : ["total"]}, # E matrix type -> sector correspondence.
        ### format {"<E type>": ["<sector1>", "<sector2>", ...]}
        "typeToSector_D": {"site" : ["total"]}, # D matrix type -> sector correspondence.
        ### format {"<D type>": ["<sector1>", "<sector2>", ...]}

        #"Ls_typeE": {"site": LS},
        #"Lt_typeD": {"site": LT},

        #--- Define footprint type for each observation operator type. ---#
        "HtypeToFootprint": {"site": "stilt"},

        #--- Prefix definition. ---#
        "emissPrior_Prefix": "Prior", # Prior emission files prefix.
        "emissProc_Prefix": "Proc", # Processing emission files prefix.
        "emissPost_Prefix": "BISOUT", # Posterior emission files prefix.
        "sigmaPrior_Prefix": "Sigma", # Uncertainty (stardard deviation / sigma) of prior emission files prefix.
        "sigmaProc_Prefix": "SigmaProc", # Uncertainty processing files prefix.
        "H_Prefix": "H",  # H matrix files prefix.
        "obs_Prefix": "OBS", # Observation files prefix.
        "D_Prefix": "D",
        "E_Prefix": "E",
        "Dbck_Prefix": "Dbck",
        "Ebck_Prefix": "Ebck",
        "R_Prefix": "R",
        "Qbck_Prefix": "Qbck",
        "initHx_Prefix": "Hx_init",
        "initHxBck_Prefix": "HxBck_init",
        "finalHx_Prefix": "Hx_final",
        "finalHxBck_Prefix": "HxBck_final",

        "postHourly_Prefix": "BIS-Hourly",
        "postDaily_Prefix": "BIS-Daily",
        "postWeekly_Prefix": "BIS-Weekly",
        "postMonthly_Prefix": "BIS-Monthly",
        "postAll_Prefix": "BIS-All",
        "spostHourly_Prefix": "BIS-Sigma-Hourly",
        "spostDaily_Prefix": "BIS-Sigma-Daily",
        "spostWeekly_Prefix": "BIS-Sigma-Weekly",
        "spostMonthly_Prefix": "BIS-Sigma-Monthly",
        "spostAll_Prefix": "BIS-Sigma-All",

        "priorHourly_Prefix": "Prior-Hourly",
        "priorDaily_Prefix": "Prior-Daily",
        "priorWeekly_Prefix": "Prior-Weekly",
        "priorMonthly_Prefix": "Prior-Monthly",
        "priorAll_Prefix": "Prior-All",
        "spriorHourly_Prefix": "Prior-Sigma-Hourly",
        "spriorDaily_Prefix": "Prior-Sigma-Daily",
        "spriorWeekly_Prefix": "Prior-Sigma-Weekly",
        "spriorMonthly_Prefix": "Prior-Sigma-Monthly",
        "spriorAll_Prefix": "Prior-Sigma-All",

        #--- Case subdirectories definition. ---#
        "emissDir": EMISSDIR, # Directory of emission, defaulted as same as above.
        "emissPriorDir": EMISSDIR + "/prior", # Directory of prior emission.
        "emissProcDir": EMISSDIR + "/processing", # Directory of processing emission.
        "emissPostDir": EMISSDIR + "/posterior", # Direcotory of posterior emission.
        "sigmaDir": SIGMADIR, # Directory of the uncertainty files.
        "sigmaPriorDir": SIGMADIR + "/prior", # Directory of the prior uncertainty files.
        "sigmaProcDir": SIGMADIR + "/processing", # Directory of hte processing uncertainty files.
        "HDir": HDIR, # Directory of H files, defaulted as same as above.
        "DDir": DDIR, # Directory of D files, defaulted as same as above.
        "DbckDir": DDIR,
        "EDir": EDIR, # Directory of E files, defaulted as same as above.
        "EbckDir": EDIR,
        "RDir": RDIR, # Directory of R files, defaulted as same as above.
        "obsDir": OBSDIR, # Directory of observation files, defaulted as same as above.
        "initConcDir": INITCONCDIR,
        "finalConcDir": FINALCONCDIR,
        "outputDir": OUTDIR,
        "outHourlyDir": OUTDIR + "/hourly",
        "outDailyDir": OUTDIR + "/daily",
        "outWeeklyDir": OUTDIR + "/weekly",
        "outMonthlyDir": OUTDIR + "/monthly",
        "outAllDir": OUTDIR + "/all",

        "output_timeScales": ["hourly", "daily", "weekly", "monthly", "all"]

    },

    #--- Setting for ExpOSSE ---#
    "OSSE":
    {
        "emissTruth_Prefix": "Truth", # Truth emission files prefix.
        "emissTruthDir": EMISSDIR + "/truth", # Directory of truth emission files.
        "truthHourly_Prefix": "Truth-Hourly",
        "truthDaily_Prefix": "Truth-Daily",
        "truthWeekly_Prefix": "Truth-Weekly",
        "truthMonthly_Prefix": "Truth-Monthly",
        "truthAll_Prefix": "Truth-All",
    },

    #--- Setting for ExpREAL ---#
    "REAL":
    {
        
        "bckPrior_Prefix": "BCK",
        "sbckPrior_Prefix": "BCK_Error",
        "bckProc_Prefix": "BCK_Proc",
        "sbckProc_Prefix": "BCK_Prioc_Error",
        "bckPriorDir": BCKDIR + "/prior",
        "bckProcDir": BCKDIR + "/processing",
        "sbckPriorDir": SBCKDIR + "/prior",
        "sbckProcDir": SBCKDIR + "/processing",

    },

    ###########################

    #--- Setting for "one_another" OSSE experiment. ---#
    "one_another":
    {
        #--- Prior emission setting. ---#
        "PriorName": "constant",
        #"PriorName": "MEIC",
        "PriorKwargs": {"const": 0},

        #--- Truth emission setting. ---#
        "TruthName": "constant",
        "TruthKwargs": {"const": 1},

        #--- Prior sigma setting. ---#
        "SigmaMethod": "constant",
        "SigmaKwargs": {"const": 10},

        #--- E matrix setting. ---#
        "EMethod": "expdis2D",
        "EKwargs": {"Ls": LS},

        #--- D matrix setting. ---#
        "DMethod": "expLag",
        "DKwargs": {"Lt": LT},

        #--- R matrix setting. ---#
        "RMethod": "diagError",
        "RKwargs": {"error_inReceptors": RECEPTORS_ERROR}
    },

    "const_emiss":
    {

        #--- Truth emission setting. ---#
        "TruthName": "MEIC",
        "TruthKwargs": {},

        #--- Prior sigma setting. ---#
        "SigmaMethod": "difference",
        "SigmaKwargs": {},

        #--- E matrix setting. ---#
        "EMethod": "expdis2D",
        "EKwargs": {"Ls": LS},

        #--- D matrix setting. ---#
        "DMethod": "expLag",
        "DKwargs": {"Lt": LT},

        #--- R matrix setting. ---#
        "RMethod": "diagError",
        "RKwargs": {"error_inReceptors": RECEPTORS_ERROR},
    },

    "independent":
    {
        #--- Observation option. ---#
        "OBS": "picarro",

        #--- Background option. ---#
        #"BCK": "wrfbck",
        "BCK": "constant",
        
        #--- Prior emission setting. ---#
        "PriorName": "constant",
        #"PriorName": "MEIC",
        "PriorKwargs": {"const": 0},

        #--- Prior sigma setting. ---#
        "SigmaMethod": "constant",
        "SigmaKwargs": {"const": 10},

        #--- E matrix setting. ---#
        "EMethod": "expdis2D",
        "EKwargs": {"Ls": LS},

        #--- D matrix setting. ---#
        "DMethod": "expLag",
        "DKwargs": {"Lt": LT},

        #--- Ebck matrix setting. ---#
        "EbckMethod": "expdis2D",
        "EbckKwargs": {"LsBck": LSBCK},

        #--- Dbck matrix setting. ---#
        "DbckMethod": "expLag",
        "DbckKwargs": {"LtBck": LTBCK},

        #--- R matrix setting. ---#
        "RMethod": "diagError",
        "RKwargs": {"error_inReceptors": RECEPTORS_ERROR},
    },

    "independent_priorFnet": 
    {
        #--- Prior fossil fuel emission setting. ---#

        # constant
        #"ffeName": "constant",
        #"ffeKwargs": {"const": 0},

        # MEIC
        "ffeName": "MEIC",
        "ffeKwargs": {},

        #--- Prior fossil fuel sigma setting. ---#

        # constant
        #"ffeSigmaMethod": "constant",
        #"ffeSigmaKwargs": {"const": 1},

        # ratio
        "ffeSigmaMethod": "emiss_ratio",
        "ffeSigmaKwargs": {"ratio": 1.0},

        #--- Prior cfta setting. ---#

        #cfta_1x1
        "cftaName": "cfta_1x1",
        "cftaKwargs": {},

        #--- Prior cfta sigma setting. ---#
        # constant
        #"cftaSigmaMethod": "constant",
        #"cftaSigmaKwargs": {"const": 2},

        # ratio
        "cftaSigmaMethod": "emiss_ratio",
        "cftaSigmaKwargs": {"ratio": 1.0},

        "ffeDir": EMISSDIR + "/fossil_fuel",
        "cftaDir": EMISSDIR + "/cfta",
        "ffeSigmaDir": SIGMADIR + "/fossil_fuel",
        "cftaSigmaDir": SIGMADIR + "/cfta",

        "ffe_Prefix": "FFE",
        "cfta_Prefix": "CFTA",
        "ffeSigma_Prefix": "FFE_Sigma",
        "cftaSigma_Prefix": "CFTA_Sigma",

        "ffeHourly_Prefix": "FFE_Hourly",
        "sffeHourly_Prefix": "SFFE_Hourly",
        "cftaHourly_Prefix": "CFTA_Hourly",
        "scftaHourly_Prefix": "SCFTA_Hourly",

        "ffeDaily_Prefix": "FFE_Daily",
        "sffeDaily_Prefix": "SFFE_Daily",
        "cftaDaily_Prefix": "CFTA_Daily",
        "scftaDaily_Prefix": "SCFTA_Daily",

        "ffeWeekly_Prefix": "FFE_Weekly",
        "sffeWeekly_Prefix": "SFFE_Weekly",
        "cftaWeekly_Prefix": "CFTA_Weekly",
        "scftaWeekly_Prefix": "SCFTA_Weekly",

        "ffeMonthly_Prefix": "FFE_Monthly",
        "sffeMonthly_Prefix": "SFFE_Monthly",
        "cftaMonthly_Prefix": "CFTA_Monthly",
        "scftaMonthly_Prefix": "SCFTA_Monthly",

        "ffeAll_Prefix": "FFE_All",
        "sffeAll_Prefix": "SFFE_All",
        "cftaAll_Prefix": "CFTA_All",
        "scftaAll_Prefix": "SCFTA_All",



    }
}


#--- Footprint settings. ---#
FootConfig = {
    #--- STILT footprint settings. ---#
    "stilt": 
    {
        "receptors": RECEPTORS, # Receptors definitions, defaulted as same as above.
        "footDir": FOOTDIR, # STILT output directory.
    }
}

#--- Observation setting. ---#
obsConfig = {
    "picarro":
    {
        "directory": PICARRODIR,
        "sites": RECEPTORS,
        "location": RECEPTORS_DIC,
        "errors": RECEPTORS_ERROR,
    }
}

#--- Background setting. ---#
bckConfig = {

    "wrfbck":
    {
        "Receptors": RECEPTORS_DIC,
        "Wrfco2Dir": "/home/tangwh/modeling/WRF-CO2-v3.0/cases/DJF2019/output/wrfco2",
        "Wrfco2Prefix": "wrfco2",
        "Domid": 1,
        "errors": 10,
        "offset": 0,
    },

    "constant":
    {
        "const": 420,
        "Receptors": RECEPTORS_DIC,
        "errors": 10,
    },

}
