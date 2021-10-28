#!/usr/bin/env python

#-> CASE of Bayesian Inversion System <-#

# Authors:
#   Wenhan TANG - 07/2021
#   ...

#from . import configure as conf

from ..experiment.get_expClass import get_expClass
from ..utils.module2dic import module2dic

class BisCase(object):

    """
        The class represents a case.
        Its object is the conductor of the whole program.
    """
    def __init__(self, conf, directoryReserve = False):

        """
            Main steps:
            1. Loading some global setting (mainly about the case) from configure module (conf).
            2. Create directory(ies) for the case.
            3. Create the object of the experiment (objExp).
            4. ...
        """

        # Step 1:
        self.CaseConfig = module2dic(conf)
        self.CASENAME = self.CaseConfig["CASENAME"]
        self.CASEDIR = self.CaseConfig["CASESDIR"] + "/" + self.CASENAME
        self.EXPTYPE = self.CaseConfig["EXPTYPE"]
        self.EXPNAME = self.CaseConfig["EXPNAME"]
        self.ISFREE = self.CaseConfig["ISFREE"]
        self.BASEEXP = self.CaseConfig["BASEEXP"]

        if self.ISFREE:
            loc = locals()
            exec("from ..experiment.exp_free." + self.EXPNAME + ".main_class import ExpFree" "_" + self.EXPNAME + " as exp")
            self.expClass = loc["exp"]

        else:
            loc = locals()
            exec("from ..experiment." + self.EXPTYPE + "." + self.EXPNAME + ".main_class import Exp" + self.EXPTYPE + "_" + self.EXPNAME + " as exp")
            self.expClass = loc["exp"]
        #self.expClass = get_expClass(EXPNAME = self.EXPNAME, EXPTYPE = self.EXPTYPE, ISFREE = self.ISFREE, BASEEXP = self.BASEEXP)
        # Step 2:
        if directoryReserve == False:
            self.mkcase()

        # Step 3:
        self.objExp = self.expClass(configure = self.CaseConfig)


    def mkcase(self):

        import os

        # Create the root directory for the case.
        if not os.path.exists(self.CaseConfig["CASEDIR"]):
            os.makedirs(self.CaseConfig["CASEDIR"])
            print("New directory " + self.CaseConfig["CASEDIR"] + " has been created.")

        os.system("cp -p " + self.CaseConfig["SRCDIR"] + "/main/configure.py " + self.CaseConfig["CASEDIR"])
        os.system("cp -rp " + self.CaseConfig["SRCDIR"] + "/main/stations " + self.CaseConfig["CASEDIR"])

        # Link the case directory to $ROOTDIR/cases.
        if os.path.exists(self.CaseConfig["CASELNK"]):
            os.system("unlink " + self.CaseConfig["CASELNK"])
        os.system("ln -sf " + self.CaseConfig["CASEDIR"] + " " + self.CaseConfig["CASELNK"])



    def app_run(self, *args, **kwargs):
        
        # Start running the experiment.
        self.objExp.exp_run(*args, **kwargs)

