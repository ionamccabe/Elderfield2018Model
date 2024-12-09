""" Comparing my recreation to the elderfield 2018 model on GitHub"""
import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model")

# My recreation
from eld2018_imhm_rec.model_funcs import *
from eld2018_imhm_rec.control_funcs import *
from eld2018_imhm_rec.plot_funcs import *
# running the simulation
numS = 1
maxCl = 1; maxCh = 1
simLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh)
simHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh)
simMix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh)

simLH_rec[2]

# GitHub Model
from fungicide_model import Parameters, SprayStrategies, Model, Simulation
# Build a parameters object, specifying which model we intend to use
params = Parameters(Model.WheatSeptoria)
# Choose which doses to apply (NB: should be halved under mixture relative to alternation)
params.highRiskDose = 1
params.lowRiskDose = 1
# Provide the chosen parameters to the simulation
params.strategy = SprayStrategies.AltLoHi
simLH_github = Simulation(params)
params.strategy = SprayStrategies.AltHiLo
simHL_github = Simulation(params)
params.strategy = SprayStrategies.Mixture
simMix_github = Simulation(params)
# Do the calculation
simLH_github.run()
simHL_github.run()
simMix_github.run()


len(simLH_github.output["S"])
len(simLH_rec[0][0][-1690:])
simLH_rec_short = tuple(arr[1211:1211+1689] for arr in simLH_rec[0])
len(simLH_rec_short[0])

simLH_github.yields

type(simLH_github.output["S"])

