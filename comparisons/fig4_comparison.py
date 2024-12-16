#%% Fig 4 plots C, D and E

import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model")

import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt


from eld2018_imhm_rec.model_funcs import systemWithControl
from eld2018_imhm_rec.control_funcs import *
from eld2018_imhm_rec.yield_funcs import *


doseHR = np.arange(0,1.025,0.025)
fsy_altLH_rec = np.array([]); fsy_altHL_rec = np.array([]); fsy_mix_rec = np.array([])
SR_altLH_rec = np.array([]); SR_altHL_rec = np.array([]); SR_mix_rec = np.array([])
maxCl = 1
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 1
    
    altLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)
    altHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)
    mix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)

    X_altLH_rec, percYield_altLH_rec, totalYield_altLH_rec, dfYield_altLH_rec, Psi_altLH_rec, SR_altLH_rec_temp, controlName_altLH_rec, t_altLH_rec = altLH_rec
    X_altHL_rec, percYield_altHL_rec, totalYield_altHL_rec, dfYield_altHL_rec, Psi_altHL_rec, SR_altHL_rec_temp, controlName_altHL_rec, t_altHL_rec = altHL_rec
    X_mix_rec, percYield_mix_rec, totalYield_mix_rec, dfYield_mix_rec, Psi_mix_rec, SR_mix_rec_temp, controlName_mix_rec, t_mix_rec = mix_rec
    
    ### Plot C
    fsy_altLH_rec_temp = totalYield_altLH_rec[0] / dfYield_altLH_rec
    fsy_altHL_rec_temp = totalYield_altHL_rec[0] / dfYield_altHL_rec
    fsy_mix_rec_temp = totalYield_mix_rec[0] / dfYield_mix_rec
    
    fsy_altLH_rec = np.append(fsy_altLH_rec, fsy_altLH_rec_temp)
    fsy_altHL_rec = np.append(fsy_altHL_rec, fsy_altHL_rec_temp)
    fsy_mix_rec = np.append(fsy_mix_rec, fsy_mix_rec_temp)
    ###
    
    ### Plot D
    SR_altLH_rec = np.append(SR_altLH_rec, SR_altLH_rec_temp)
    SR_altHL_rec = np.append(SR_altHL_rec, SR_altHL_rec_temp)
    SR_mix_rec = np.append(SR_mix_rec, SR_mix_rec_temp)
    ###

# Plot E
ly_altLH_rec = np.array([]); ly_altHL_rec = np.array([]); ly_mix_rec = np.array([])
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 30
    
    altLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    altHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    mix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    
    X_altLH_rec, percYield_altLH_rec, totalYield_altLH_rec, dfYield_altLH_rec, Psi_altLH_rec, SR_altLH_rec_temp, controlName_altLH_rec, t_altLH_rec = altLH_rec
    X_altHL_rec, percYield_altHL_rec, totalYield_altHL_rec, dfYield_altHL_rec, Psi_altHL_rec, SR_altHL_rec_temp, controlName_altHL_rec, t_altHL_rec = altHL_rec
    X_mix_rec, percYield_mix_rec, totalYield_mix_rec, dfYield_mix_rec, Psi_mix_rec, SR_mix_rec_temp, controlName_mix_rec, t_mix_rec = mix_rec
    
    ly_altLH_rec_temp = calcLifetimeYield(percYield_altLH_rec, totalYield_altLH_rec, dfYield_altLH_rec)
    ly_altHL_rec_temp = calcLifetimeYield(percYield_altHL_rec, totalYield_altHL_rec, dfYield_altHL_rec)
    ly_mix_rec_temp = calcLifetimeYield(percYield_mix_rec, totalYield_mix_rec, dfYield_mix_rec)
    
    ly_altLH_rec = np.append(ly_altLH_rec, ly_altLH_rec_temp)
    ly_altHL_rec = np.append(ly_altHL_rec, ly_altHL_rec_temp)
    ly_mix_rec = np.append(ly_mix_rec, ly_mix_rec_temp)
###

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: First Season Yield
axes[0].grid()
axes[0].set_title("First Season Yield")
axes[0].plot(doseHR, fsy_mix_rec, 'r', label='mixture')
axes[0].plot(doseHR, fsy_altHL_rec, 'g', label='Alt High Low')
axes[0].plot(doseHR, fsy_altLH_rec, 'b', label='Alt Low High')
axes[0].set_xlabel('High-Risk Dose')
axes[0].set_ylabel('First Season Yield (%)')
axes[0].legend()

# Plot 2: Selection Ratio (swapped with Lifetime Yield)
axes[1].grid()
axes[1].set_title("Selection Ratio")
axes[1].plot(doseHR, SR_mix_rec, 'r', label='mixture')
axes[1].plot(doseHR, SR_altHL_rec, 'g', label='Alt High Low')
axes[1].plot(doseHR, SR_altLH_rec, 'b', label='Alt Low High')
axes[1].set_xlabel('High-Risk Dose')
axes[1].set_ylabel('Selection Ratio')
axes[1].legend()

# Plot 3: Lifetime Yield (swapped with Selection Ratio)
axes[2].grid()
axes[2].set_title("Lifetime Yield")
axes[2].plot(doseHR, ly_mix_rec, 'r', label='mixture')
axes[2].plot(doseHR, ly_altHL_rec, 'g', label='Alt High Low')
axes[2].plot(doseHR, ly_altLH_rec, 'b', label='Alt Low High')
axes[2].set_xlabel('High-Risk Dose')
axes[2].set_ylabel('Lifetime Yield')
axes[2].legend()

# Add overall title for the entire figure
fig.suptitle('Rec', fontsize=16)

# Adjust layout
plt.tight_layout()
plt.show()

#%%
import numpy as np
from eld2018_imhm_rec.model_funcs import systemWithControl
from eld2018_imhm_rec.control_funcs import *
from tqdm import tqdm
from eld2018_imhm_rec.yield_funcs import *
import matplotlib.pyplot as plt
from fungicide_model import Parameters, SprayStrategies, Model, Simulation

doseHR = np.arange(0,1.025,0.025)
fsy_altLH_github = np.array([]); fsy_altHL_github = np.array([]); fsy_mix_github = np.array([])
SR_altLH_github = np.array([]); SR_altHL_github = np.array([]); SR_mix_github = np.array([])
# Build a parameters object, specifying which model we intend to use
params = Parameters(Model.WheatSeptoria)
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 1

    # Choose which doses to apply (NB: should be halved under mixture relative to alternation)
    params.highRiskDose = doseHighRisk
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
    
    ### Plot C
    fsy_altLH_github_temp = simLH_github.yields[0] / simLH_github.diseaseFreeYield
    fsy_altHL_github_temp = simHL_github.yields[0] / simHL_github.diseaseFreeYield
    fsy_mix_github_temp = simMix_github.yields[0] / simMix_github.diseaseFreeYield
    
    fsy_altLH_github = np.append(fsy_altLH_github, fsy_altLH_github_temp)
    fsy_altHL_github = np.append(fsy_altHL_github, fsy_altHL_github_temp)
    fsy_mix_github = np.append(fsy_mix_github, fsy_mix_github_temp)
    ###
    
    ### Plot D
    SR_altLH_github = np.append(SR_altLH_github, simLH_github.selectionRatio)
    SR_altHL_github = np.append(SR_altHL_github, simHL_github.selectionRatio)
    SR_mix_github = np.append(SR_mix_github, simMix_github.selectionRatio)
    ###

# Plot E
ly_altLH_github = np.array([]); ly_altHL_github = np.array([]); ly_mix_github = np.array([])
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 30
    
    # Choose which doses to apply (NB: should be halved under mixture relative to alternation)
    params.highRiskDose = doseHighRisk
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

    ly_altLH_github_temp = simLH_github.yieldTilCriticalLoss/simLH_github.diseaseFreeYield
    ly_altHL_github_temp = simHL_github.yieldTilCriticalLoss/simHL_github.diseaseFreeYield
    ly_mix_github_temp = simMix_github.yieldTilCriticalLoss/simMix_github.diseaseFreeYield
    
    ly_altLH_github = np.append(ly_altLH_github, ly_altLH_github_temp)
    ly_altHL_github = np.append(ly_altHL_github, ly_altHL_github_temp)
    ly_mix_github = np.append(ly_mix_github, ly_mix_github_temp)
###

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: First Season Yield
axes[0].grid()
axes[0].set_title("First Season Yield")
axes[0].plot(doseHR, fsy_mix_github, 'r', label='mixture')
axes[0].plot(doseHR, fsy_altHL_github, 'g', label='Alt High Low')
axes[0].plot(doseHR, fsy_altLH_github, 'b', label='Alt Low High')
axes[0].set_xlabel('High-Risk Dose')
axes[0].set_ylabel('First Season Yield (%)')
axes[0].legend()

# Plot 2: Selection Ratio (swapped with Lifetime Yield)
axes[1].grid()
axes[1].set_title("Selection Ratio")
axes[1].plot(doseHR, SR_mix_github, 'r', label='mixture')
axes[1].plot(doseHR, SR_altHL_github, 'g', label='Alt High Low')
axes[1].plot(doseHR, SR_altLH_github, 'b', label='Alt Low High')
axes[1].set_xlabel('High-Risk Dose')
axes[1].set_ylabel('Selection Ratio')
axes[1].legend()

# Plot 3: Lifetime Yield (swapped with Selection Ratio)
axes[2].grid()
axes[2].set_title("Lifetime Yield")
axes[2].plot(doseHR, ly_mix_github, 'r', label='mixture')
axes[2].plot(doseHR, ly_altHL_github, 'g', label='Alt High Low')
axes[2].plot(doseHR, ly_altLH_github, 'b', label='Alt Low High')
axes[2].set_xlabel('High-Risk Dose')
axes[2].set_ylabel('Lifetime Yield')
axes[2].legend()

# Add overall title for the entire figure
fig.suptitle('GitHub', fontsize=16)

# Adjust layout
plt.tight_layout()
plt.show()

# %%
import pandas as pd

fig4_df = pd.DataFrame()

# Percent diff of fsy
LH_fsy_avg = (fsy_altLH_github + fsy_altLH_rec)/2
LH_fsy_abs = abs(fsy_altLH_github - fsy_altLH_rec)
LH_fsy_percDiff = LH_fsy_abs / LH_fsy_avg * 100
fig4_df['FSY LH'] = LH_fsy_percDiff

HL_fsy_avg = (fsy_altHL_github + fsy_altHL_rec)/2
HL_fsy_abs = abs(fsy_altHL_github - fsy_altHL_rec)
HL_fsy_percDiff = HL_fsy_abs / HL_fsy_avg * 100
fig4_df['FSY HL'] = HL_fsy_percDiff

mix_fsy_avg = (fsy_mix_github + fsy_mix_rec)/2
mix_fsy_abs = abs(fsy_mix_github - fsy_mix_rec)
mix_fsy_percDiff = mix_fsy_abs / mix_fsy_avg * 100
fig4_df['FSY Mix'] = mix_fsy_percDiff

# Percent diff of SR
LH_SR_avg = (SR_altLH_github + SR_altLH_rec)/2
LH_SR_abs = abs(SR_altLH_github - SR_altLH_rec)
LH_SR_percDiff = LH_SR_abs / LH_SR_avg * 100
fig4_df['SR LH'] = LH_SR_percDiff

HL_SR_avg = (SR_altHL_github + SR_altHL_rec)/2
HL_SR_abs = abs(SR_altHL_github - SR_altHL_rec)
HL_SR_percDiff = HL_SR_abs / HL_SR_avg * 100
fig4_df['SR HL'] = HL_SR_percDiff

mix_SR_avg = (SR_mix_github + SR_mix_rec)/2
mix_SR_abs = abs(SR_mix_github - SR_mix_rec)
mix_SR_percDiff = mix_SR_abs / mix_SR_avg * 100
fig4_df['SR Mix'] = mix_SR_percDiff

# Percent diff of ly
LH_ly_avg = (ly_altLH_github + ly_altLH_rec)/2
LH_ly_abs = abs(ly_altLH_github - ly_altLH_rec)
LH_ly_percDiff = LH_ly_abs / LH_ly_avg * 100
fig4_df['LTY LH'] = LH_ly_percDiff

HL_ly_avg = (ly_altHL_github + ly_altHL_rec)/2
HL_ly_abs = abs(ly_altHL_github - ly_altHL_rec)
HL_ly_percDiff = HL_ly_abs / HL_ly_avg * 100
fig4_df['LTY HL'] = HL_ly_percDiff

mix_ly_avg = (ly_mix_github + ly_mix_rec)/2
mix_ly_abs = abs(ly_mix_github - ly_mix_rec)
mix_ly_percDiff = mix_ly_abs / mix_ly_avg * 100
fig4_df['LTY Mix'] = mix_ly_percDiff

# DF yield
LH_df_abs = abs(simLH_github.diseaseFreeYield - dfYield_altLH_rec/2)
HL_df_abs = abs(simHL_github.diseaseFreeYield - dfYield_altHL_rec/2)
mix_df_abs = abs(simMix_github.diseaseFreeYield - dfYield_mix_rec/2)

LH_df_avg = (simLH_github.diseaseFreeYield + dfYield_altLH_rec/2)/2
HL_df_avg = (simHL_github.diseaseFreeYield + dfYield_altHL_rec/2)/2
mix_df_avg = (simMix_github.diseaseFreeYield + dfYield_mix_rec/2)/2

LH_df_percDiff = LH_df_abs / LH_df_avg * 100
HL_df_percDiff = HL_df_abs / HL_df_avg * 100
mix_df_percDiff =mix_df_abs / mix_df_avg * 100

fig4_df['DFY LH'] = LH_df_percDiff
fig4_df['DFY HL'] = HL_df_percDiff
fig4_df['DFY Mix'] = mix_df_percDiff



# Calculate the max value for each column
max_values = fig4_df.max()

# Calculate the count of values >= 5 for each column
count_values_ge_5 = (fig4_df >= 5).sum()

# Create a new DataFrame
fig4_max_df = pd.DataFrame({
    'Max Value': max_values,
    'Count >= 5': count_values_ge_5
}).transpose()

# Rename the rows to indicate what they represent
fig4_max_df.index = ['Max Value', 'Count >= 5']

print('DFY GitHub: ',simLH_github.diseaseFreeYield,'    DFY Rec: ',dfYield_altLH_rec/2)

# %%
