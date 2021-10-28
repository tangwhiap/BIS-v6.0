#!/usr/bin/env python
# [===============> Showing the progress... ] 40.0%
# Author: Wenhan TANG - 05/2021

import numpy as np
last = "None"
def show_progress(i, N):
    global last
    Nshow = 80
    #if np.mod(i, 100) != 0:
    #    return
    progress = (i + 1) / N
    str_progress = "%4.1f" % (progress * 100)
    if str_progress == last:
        if i >= N-1:
            print("")
        return
    Nok = round(Nshow * progress)
    if Nok < Nshow:
        string = "[" + ("=") * Nok + ">" + (" ") * (Nshow - Nok - 1) + "] " + str_progress + "%"
    elif i >= N-1:
        string = "[" + ("=") * Nok + (" ") * (Nshow - Nok - 1) + "] " + str_progress + "% Complete!"
    else:
        string = "[" + ("=") * Nok + (" ") * (Nshow - Nok - 1) + "] " + str_progress + "%"
    print(string, end = "\r")
    if i >= N-1:
        print("")
    last = str_progress
    return

# for testing ...
if __name__ == "__main__":
    print("Computing ...")
    for i in range(1000000):
        show_progress(i, 1000000)
