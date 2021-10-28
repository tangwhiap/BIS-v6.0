#!/usr/bin/env python

# Authors:
#    Wenhan TANG - 08/2021
#    ...

class LocPoint(object):

    def __init__(self, staName, staLon, staLat):

        self.staName = staName
        self.staLon = staLon
        self.staLat = staLat
        self.parentLON = None
        self.parentLAT = None
        self.isAnchor = False
        #self.nLat = LAT.shape[0]        
        #self.nLon = LON.shape[1]
        #self.parentLonlist = self.parentLON[0]
        #self.parentLatlist = self.parentLAT[:, 0]
        #assert self.isIn(), "This point: " + self.staName + " is out of range."
        #self.anchor()

    def isIn(self):
        if self.staLon >= self.parentLonlist[0] and self.staLon <= self.parentLonlist[-1]:
            if self.staLat >= self.parentLatlist[0] and self.staLat <= self.parentLatlist[-1]:
                return True
        return False

    def parentLONLAT_check(self, LON, LAT):
        lonList = LON[0]
        latList = LAT[:, 0]
        if np.sum(np.abs(lonList - self.parentLonlist)) <= 1e-6:
            if np.sum(np.abs(latList - self.parentLatlist)) <= 1e-6:
                return True
        return False

    def anchor(self, LON, LAT):
        assert LON.shape == LAT.shape
        lon = self.staLon; lat = self.staLat
        
        self.nLat = LAT.shape[0]
        self.nLon = LON.shape[1]
        self.parentLON = LON
        self.parentLAT = LAT
        self.parentLonlist = self.parentLON[0]
        self.parentLatlist = self.parentLAT[:, 0]
        assert self.isIn(), "This point: " + self.staName + " is out of range."
        lonList = self.parentLonlist
        latList = self.parentLatlist
        #dLonList = lon - lonList; dLatList = lat - latList
        #sign_dLonList = dLonList[1:] * dLonList[:-1]
        #sign_dLatList = dLatList[1:] * dLatList[:-1]
        #assert np.sum(sign_dLonList < 0) == 1
        for iLon, leftLon, rightLon in zip(range(len(lonList) - 1), lonList[:-1], lonList[1:]):
            if leftLon <= lon and lon <= rightLon:
                self.wel_i = iLon; self.wer_i = iLon + 1
                self.wel_lon = leftLon; self.wer_lon = rightLon

        for iLat, southLat, northLat in zip(range(len(latList) - 1), latList[:-1], latList[1:]):
            if southLat <= lat and lat <= northLat:
                self.sns_i = iLat; self.snn_i = iLat + 1
                self.sns_lat = southLat; self.snn_lat = northLat
        
        self.wel = lon - self.wel_lon; self.wer = self.wer_lon - lon
        self.sns = lat - self.sns_lat; self.snn = self.snn_lat - lat

        self.we_nearest_i = self.wel_i if self.wel < self.wer else self.wer_i
        self.we_nearest_lon = self.wel_lon if self.wel < self.wer else self.wer_lon
        self.sn_nearest_i = self.sns_i if self.sns < self.snn else self.snn_i
        self.sn_nearest_lat = self.sns_lat if self.sns < self.snn else self.snn_lat

        self.isAnchor = True

        """
         (wel)  (wer) 
        +--------------+ snn_lat
        |A   |        B|
        |    |         | (snn)
        |----o---------|
        |    |         | (sns)
        |D   |        C|
        +--------------+ sns_lat
      wel_lon        wer_lon
        
        """

    def nearest(self, arr):
        assert arr.shape == (self.nLat, self.nLon)
        return arr[self.sn_nearest_i, self.we_nearest_i]

    def belinear(self, arr):

        wel_i = self.wel_i
        wer_i = self.wer_i
        sns_i = self.sns_i
        snn_i = self.snn_i

        wel = self.wel
        wer = self.wer
        sns = self.sns
        snn = self.snn

        wer_lon = self.wer_lon
        wel_lon = self.wel_lon
        sns_lat = self.sns_lat
        snn_lat = self.snn_lat

        assert arr.shape == (self.nLat, self.nLon)
        valueA = arr[snn_i, wel_i]
        valueB = arr[snn_i, wer_i]
        valueC = arr[sns_i, wer_i]
        valueD = arr[sns_i, wel_i]
        S = (wer_lon - wel_lon) * (snn_lat - sns_lat)

        weightA = sns * wer / S
        weightB = sns * wel / S
        weightC = snn * wel / S
        weightD = snn * wer / S

        return valueA * weightA + valueB * weightB + valueC * weightC + valueD * weightD
        
         
        
