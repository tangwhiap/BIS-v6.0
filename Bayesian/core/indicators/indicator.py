#!/usr/bin/env python

# Authors:
#   Wenhan TANG - 07/2021
#   ...

class Indicator(object):
    """
         Base class of the indicator for a type of matries group combined with different sectors/types. 
    """
    def __init__(self, *args, **kwargs):
        assert False, "Defined in the child class."
        pass

    def dim_val(self):
        for sect in self.typeNames:
            assert self.data[sect].shape == (self.dim1Dict[sect], self.dim2Dict[sect])

    def check(self):
        self.dim_val()

    def get_data(self):
        return self.data

    def get_matrix(self, sect):
        return self.data[sect]

    def get_dim(self, sect):
        return self.dim1Dict[sect], self.dim2Dict[sect]



