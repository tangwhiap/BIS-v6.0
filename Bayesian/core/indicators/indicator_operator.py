#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...


from .indicator import Indicator
from .types_indicator import TypesIndicator
from .sectors_indicator import SectorsIndicator

from pdb import set_trace

import numpy as np

class IndicatorOperator(object):

    def __init__(self, *args, **kwargs):

        if "data" in kwargs:
            self.data = kwargs["data"]
            self.typeNames = list(self.data.keys())
            self.indicatorClass = "sectors" # Default class of indicator included in.

            if "indicator_class" in kwargs:
                self.indicatorClass = kwargs["indicator_class"]

            if "sectorToType" in kwargs:
                self.sectorToType = kwargs["sectorToType"]

            if "typeToSector" in kwargs:
                self.typeToSector = kwargs["typeToSector"]

        elif "indicator" in kwargs:
            objIndicator = kwargs["indicator"]
            assert isinstance(objIndicator, Indicator)
            self.data = objIndicator.data

            if isinstance(objIndicator, TypesIndicator):
                self.indicatorClass = "types"

            if isinstance(objIndicator, SectorsIndicator):
                self.indicatorClass = "sectors"

            if self.indicatorClass == "types":
                self.typeToSector = objIndicator.typeToSector
                self.sectorToType = objIndicator.sectorToType

            self.typeNames = objIndicator.typeNames

        else:
            assert False, "Error! Class \"IndicatorOperator\" constructor: Arguments error."
            

        self.nType = len(self.typeNames)
        assert self.indicatorClass  in ["types", "sectors"],  "Error! IndicatorOperator.indicatorClass must be one of \"types\" and \"sectors\"."

    def get_data(self):
        return self.data

    def get_matrix(self, sect):
        return self.data[sect]

    def get_sum(self, sectList = None):
        if sectList is None:
            sectList = self.typeNames
        res = 0
        for sect in sectList:
            res += self.data[sect]
        return res

    def __caculate_operator(self, cac_left, cac_right, strOperator):

        if isinstance(cac_left, np.ndarray) and isinstance(cac_right, np.ndarray) and strOperator == "*":
            return np.matmul(cac_left, cac_right)
        loc = locals()
        exec("result = cac_left " + strOperator + " cac_right") 
        result = loc["result"]
        return result

    def __caculator_with_ndarray(self, other, strOperator, self_in = "left"):

        self_in = self_in.lower()
        assert self_in in ["left", "right"]
        isLeft = self_in == "left"
        resData = {}

        for sector in self.typeNames:
            resData[sector] = self.__caculate_operator(self.data[sector], other, strOperator) if isLeft else self.__caculate_operator(other, self.data[sector], strOperator)

        return IndicatorOperator(data = resData, indicator_class = "sectors")

    def __caculator_self_in_left(self, other, strOperator):

        assert isinstance(other, IndicatorOperator)
        assert self.indicatorClass == "sectors" or other.indicatorClass == "sectors", "Error! two objects of \"IndicatorOperator\" couldn't be both \"types\" class"

        if self.indicatorClass == "types":
            # Case: types * sectors
            resData = {}
            for sector in other.typeNames:
                #resData[sector] = self.data[self.sectorToType[sector]] * other.data[sector]
                #exec("resData[sector] = self.data[self.sectorToType[sector]] " + strOperator + " other.data[sector]")
                resData[sector] = self.__caculate_operator(self.get_matrix(self.sectorToType[sector]), other.get_matrix(sector), strOperator)


        elif other.indicatorClass == "sectors":
            # Case: sectors * sectors
            resData = {}
            for sector in self.typeNames:
                assert sector in other.typeNames, "Error! Sector: " + sector + " in left object couldn't be found in right object."
                #resData[sector] = self.data[sector] * other.data[sector]
                #exec("resData[sector] = self.data[sector] " + strOperator + " other.data[sector]")
                resData[sector] = self.__caculate_operator(self.get_matrix(sector), other.get_matrix(sector), strOperator)

        elif other.indicatorClass == "types":
            # Case: sectors * types
            resData = {}
            for sector in self.typeNames:
                #resData[sector] = self.data[sector] * other.data[other.sectorToType[sector]]
                #exec("resData[sector] = self.data[sector] " + strOperator + " other.data[other.sectorToType[sector]]")
                resData[sector] = self.__caculate_operator(self.get_matrix(sector), other.get_matrix(other.sectorToType[sector]), strOperator)

        else:
            assert False, "Error! Logical contradiction."

        return IndicatorOperator(data = resData, indicator_class = "sectors")

    def __getitem__(self, key):
        resData = {}
        for sector in self.typeNames:
            resData[sector] = self.data[sector][key]
        return IndicatorOperator(data = resData, indicator_class = "sectors")

    def __mul__(self, other):

        if isinstance(other, IndicatorOperator_Zero):
            # a * 0 = 0
            return other
        if isinstance(other, IndicatorOperator_One):
            # a * 1 = a
            return self

        if isinstance(other, np.ndarray):
            return self.__caculator_with_ndarray(other, "*")

        return self.__caculator_self_in_left(other, "*")

    def __rmul__(self, other):
        #print("__rmul__")
        #print(other)
        #print(type(other))
        if isinstance(other, np.ndarray):
            return self.__caculator_with_ndarray(other, "*", self_in = "right")
        else:
            assert False

    def __add__(self, other):
        if isinstance(other, IndicatorOperator_Zero):
            # a + 0 = a
            return self
        return self.__caculator_self_in_left(other, "+")

    def __sub__(self, other):
        if isinstance(other, IndicatorOperator_Zero):
            # a + 0 = a
            return self
        return self.__caculator_self_in_left(other, "-")

    def __radd__(self, other):
        return self.__add__(other)

    
        
class IndicatorOperator_Zero(object):
    """
        An abstract class representing the "Zero".
        The object of this class can participate in calculation with objects of "IndicatorOperator".
    """
    pass

class IndicatorOperator_One(object):
    """
        An abstract class representing the "One".
        The object of this class can participate in calculation with objects of "IndicatorOperator".
    """
    pass
    #def __rmul__(self, other):
#        assert isinstance(other, IndicatorOperator)
#        assert self.indicatorClass == "sectors" or other.indicatorClass == "sectors",
#            "Error! two objects of \"IndicatorOperator\" couln'd be both \"types\" class"
#        if self.indicatorClass == "sectors":
#            return self.__mul__(other)


                
    
