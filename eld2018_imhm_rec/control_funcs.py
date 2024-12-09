import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model\eld2018_imhm_rec")

from params_funcs import t_GS32,t_GS39

def noControl(time,cl_arr,ch_arr, maxCl, maxCh):
    controlName = 'No Control'
    Cl = 0; Ch = 0
    return Cl, Ch, controlName

def altLowHigh(time, cl_arr, ch_arr, maxCl, maxCh):
    controlName = 'Alt Low-High'
    if time < t_GS32:
        Cl = 0; Ch = 0
    elif time == t_GS32:
        Cl = maxCl
        Ch = 0
        # print('Added Cl at t_GS32, time:', t_GS32)
    elif time == t_GS39:
        Cl = cl_arr[-1]
        Ch = maxCh
        # print('Added Ch at t_GS39, time:', t_GS39)
    else:
        Cl = cl_arr[-1]
        Ch = ch_arr[-1]
    return Cl, Ch, controlName

def altHighLow(time, cl_arr, ch_arr, maxCl, maxCh):
    controlName = 'Alt High-Low'
    if time < t_GS32:
        Cl = 0; Ch = 0
    elif time == t_GS32:
        Cl = 0
        Ch = maxCh
        # print('Added Ch at t_GS32, time:', t_GS32)
    elif time == t_GS39:
        Cl = maxCl
        Ch = ch_arr[-1]
        # print('Added Cl at t_GS39, time:', t_GS39)
    else:
        Cl = cl_arr[-1]
        Ch = ch_arr[-1]
    return Cl, Ch, controlName

def mixture(time, cl_arr, ch_arr, maxCl, maxCh):
    controlName = 'Mixture'
    if time < t_GS32:
        Cl = 0; Ch = 0
    elif time == t_GS32:
        Cl = maxCl/2
        Ch = maxCh/2
        # print('Added Cl and Ch at t_GS32, time:', t_GS32)
    elif time == t_GS39:
        Cl = maxCl/2 + cl_arr[-1]
        Ch = maxCh/2 + ch_arr[-1]
        # print('Added Cl and Ch at t_GS39, time:', t_GS39)
    else:
        Cl = cl_arr[-1]
        Ch = ch_arr[-1]
    return Cl, Ch, controlName