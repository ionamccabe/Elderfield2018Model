""" Comparing my recreation to the elderfield 2018 model on GitHub"""
import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model")

import pandas as pd

numS = 30

### My recreation
from eld2018_imhm_rec.model_funcs import *
from eld2018_imhm_rec.control_funcs import *
from eld2018_imhm_rec.plot_funcs import *
# running the simulation
maxCl = 1; maxCh = 1
simLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh)
simHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh)
simMix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh)

### GitHub Model
from fungicide_model import Parameters, SprayStrategies, Model, Simulation
# Build a parameters object, specifying which model we intend to use
params = Parameters(Model.WheatSeptoria)
# Choose which doses to apply (NB: should be halved under mixture relative to alternation)
params.highRiskDose = 1
params.lowRiskDose = 1
params.maxSeasons = numS
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

### Making DFs to compare state classes (raw nums and perc diff)
state_classes = ['S', 'ER', 'ES', 'IR', 'IS', 'R', 'PR', 'PS']
LH_state_classes_df = pd.DataFrame()
LH_perc_diff_df = pd.DataFrame()
HL_state_classes_df = pd.DataFrame()
HL_perc_diff_df = pd.DataFrame()
Mix_state_classes_df = pd.DataFrame()
Mix_perc_diff_df = pd.DataFrame()
for i, state in enumerate(state_classes):
    LH_state_classes_df[f'GitHub {state}'] = simLH_github.output[state]
    LH_state_classes_df[f'Rec {state}'] = simLH_rec[0][i]
    HL_state_classes_df[f'GitHub {state}'] = simHL_github.output[state]
    HL_state_classes_df[f'Rec {state}'] = simHL_rec[0][i]
    Mix_state_classes_df[f'GitHub {state}'] = simMix_github.output[state]
    Mix_state_classes_df[f'Rec {state}'] = simMix_rec[0][i]

    LH_avg = (simLH_github.output[state] + simLH_rec[0][i])/2
    LH_abs_diff = abs(simLH_github.output[state] - simLH_rec[0][i])
    LH_perc_diff_df[f'Diff {state}'] = LH_abs_diff / LH_avg * 100

    HL_avg = (simHL_github.output[state] + simHL_rec[0][i])/2
    HL_abs_diff = abs(simHL_github.output[state] - simHL_rec[0][i])
    HL_perc_diff_df[f'Diff {state}'] = HL_abs_diff / HL_avg * 100

    Mix_avg = (simMix_github.output[state] + simMix_rec[0][i])/2
    Mix_abs_diff = abs(simMix_github.output[state] - simMix_rec[0][i])
    Mix_perc_diff_df[f'Diff {state}'] = Mix_abs_diff / Mix_avg * 100

### Making Percent Difference Yield DF
yield_df = pd.DataFrame()

# calculations for perc diff
LH_avg_yield = (simLH_github.yields + simLH_rec[2] / 2)/2 # LH
LH_abs_yield = abs(simLH_github.yields - simLH_rec[2] / 2)
HL_avg_yield = (simHL_github.yields + simHL_rec[2] / 2)/2 # HL
HL_abs_yield = abs(simHL_github.yields - simHL_rec[2] / 2)
Mix_avg_yield = (simMix_github.yields + simMix_rec[2] / 2)/2 # Mix
Mix_abs_yield = abs(simMix_github.yields - simMix_rec[2] / 2)

# Adding perc diff for all
yield_df['Perc Diff (Mod Rec) LH'] = LH_abs_yield / LH_avg_yield * 100
yield_df['Perc Diff (Mod Rec) HL'] = HL_abs_yield / HL_avg_yield * 100
yield_df['Perc Diff (Mod Rec) Mix'] = Mix_abs_yield / Mix_avg_yield * 100

# Adding LH
yield_df['GitHub Yield LH'] = simLH_github.yields
yield_df['Modified Rec Yield LH'] = simLH_rec[2] / 2
yield_df['Rec Yield LH'] = simLH_rec[2]

# Adding HL
yield_df['GitHub Yield HL'] = simHL_github.yields
yield_df['Modified Rec Yield HL'] = simHL_rec[2] / 2
yield_df['Rec Yield HL'] = simHL_rec[2]

# Adding Mix
yield_df['GitHub Yield Mix'] = simMix_github.yields
yield_df['Modified Rec Yield Mix'] = simMix_rec[2] / 2
yield_df['Rec Yield Mix'] = simMix_rec[2]

### Making DF to Compare LY Yield
# finding lty for my model
X_altLH, percYield_altLH, totalYield_altLH, dfYield_altLH, Psi_altLH, SR_altLH_temp, controlName_altLH, t_altLH = simLH_rec
X_altHL, percYield_altHL, totalYield_altHL, dfYield_altHL, Psi_altHL, SR_altHL_temp, controlName_altHL, t_altHL = simHL_rec
X_mix, percYield_mix, totalYield_mix, dfYield_mix, Psi_mix, SR_mix_temp, controlName_mix, t_mix = simMix_rec

lty_LH_rec = calcLifetimeYield(percYield_altLH,totalYield_altLH,dfYield_altLH)
lty_HL_rec = calcLifetimeYield(percYield_altHL,totalYield_altHL,dfYield_altHL)
lty_mix_rec = calcLifetimeYield(percYield_mix,totalYield_mix,dfYield_mix)

# making a df
# LH_avg_lty = (simLH_github.yieldTilCriticalLoss + lty_LH_rec)/2 # LH
# LH_abs_lty = abs(simLH_github.yields - lty_LH_rec)
print('LH Yield to Critical Loss GitHub: ', simLH_github.yieldTilCriticalLoss)
print('LH LTY Calc (YTCL / DFY) GitHub: ', simLH_github.yieldTilCriticalLoss/simLH_github.diseaseFreeYield)
print('LH Lifetime Yield Rec: ', lty_LH_rec)
# lty_df = pd.DataFrame()
# lty_df['Perc Diff LH'] = 



