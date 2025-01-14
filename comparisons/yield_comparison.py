""" 2024-12-16, compare lifetime yields"""

# %% Import Functions and run models
### Setup
import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model")
# Recreation
from eld2018_imhm_rec.model_funcs import *
from eld2018_imhm_rec.control_funcs import *
from eld2018_imhm_rec.plot_funcs import *
# GitHub
from fungicide_model import Parameters, SprayStrategies, Model, Simulation

### Running Models
numS = 30
# Recreation
maxCl = 1; maxCh = 1
simLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh)
simHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh)
simMix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh)
# GitHub
params = Parameters(Model.WheatSeptoria)
params.highRiskDose = 1
params.lowRiskDose = 1
params.maxSeasons = numS
params.strategy = SprayStrategies.AltLoHi
simLH_github = Simulation(params)
simLH_github.run()
params.strategy = SprayStrategies.AltHiLo
simHL_github = Simulation(params)
simHL_github.run()
params.strategy = SprayStrategies.Mixture
simMix_github = Simulation(params)
simMix_github.run()

### Function setup
def percDiff(val1, val2):
    return abs(val1 - val2) / ((val1+val2)/2) * 100

#%% Calculate and Compare Yield over 30 Seasons
""" Calc Percent Diff in Yield for Each Season """
import pandas as pd

# Calculate the Diff and Percent Diff for All Treatments
LH_diff = abs(simLH_rec[2]/2 - simLH_github.yields)
LH_perc_diff = percDiff(simLH_rec[2]/2, simLH_github.yields)

HL_diff = abs(simHL_rec[2]/2 - simHL_github.yields)
HL_perc_diff = percDiff(simHL_rec[2]/2, simHL_github.yields)

Mix_diff = abs(simMix_rec[2]/2 - simMix_github.yields)
Mix_perc_diff = percDiff(simMix_rec[2]/2, simMix_github.yields)

# Create a DataFrame
A_yield_diff_df = pd.DataFrame({ # labelled A for easy access (alphabetically first)
    'LH_perc_diff': LH_perc_diff,
    'HL_perc_diff': HL_perc_diff,
    'Mix_perc_diff': Mix_perc_diff,
    'LH_diff': LH_diff,
    'HL_diff': HL_diff,
    'Mix_diff': Mix_diff
})

# Compare DFY
dfyRec = simLH_rec[3]/2
dfyGithub = simLH_github.diseaseFreeYield
print('DF Yield Percent Diff: ', percDiff(dfyRec,dfyGithub))
print('Rec DFY: ', dfyRec, '    GitHub DFY: ', dfyGithub)

#%% Calculate and Compare Percent Yield over 30 Seasons
""" Calc Percent Yield for Each Season """
LH_perc_yield_rec = (simLH_rec[2]/2) / dfyRec * 100
LH_perc_yield_github = simLH_github.yields / dfyGithub * 100

HL_perc_yield_rec = (simHL_rec[2]/2) / dfyRec * 100
HL_perc_yield_github = simHL_github.yields / dfyGithub * 100

Mix_perc_yield_rec = (simMix_rec[2]/2) / dfyRec * 100
Mix_perc_yield_github = simMix_github.yields / dfyGithub * 100

# %% Calculate and Compare Lifetime Yield
""" Calc Lifetime Yield """
LH_lty_rec = calcLifetimeYield((simLH_rec[2]/2),dfyRec)
HL_lty_rec = calcLifetimeYield((simHL_rec[2]/2),dfyRec)
Mix_lty_rec = calcLifetimeYield((simMix_rec[2]/2),dfyRec)

LH_lty_github = calcLifetimeYield(simLH_github.yields,dfyGithub)
HL_lty_github = calcLifetimeYield(simHL_github.yields,dfyGithub)
Mix_lty_github = calcLifetimeYield(simMix_github.yields,dfyGithub)

print('Lifetime Yield')
print('Low-High -  Rec:',LH_lty_rec,' GitHub:',LH_lty_github)
print('High-Low -  Rec:',HL_lty_rec,' GitHub:',HL_lty_github)
print('Mixture  -  Rec:',Mix_lty_rec,' GitHub:',Mix_lty_github)

#%% Calculate and Compare LTY as maxCH changes
""" Calc LTY as maxCh Changes """

from tqdm import tqdm

doseHR = np.arange(0,1.025,0.025)
maxCl = 1
numS = 30

params = Parameters(Model.WheatSeptoria)
params.lowRiskDose = maxCl
params.maxSeasons = numS

ly_altLH_rec = np.array([]); ly_altHL_rec = np.array([]); ly_mix_rec = np.array([])
ly_altLH_github = np.array([]); ly_altHL_github = np.array([]); ly_mix_github = np.array([])
for dose in tqdm(doseHR):
    maxCh = dose
    # Recreation
    simLH_rec = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh, printDone=False)
    simHL_rec = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh, printDone=False)
    simMix_rec = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh, printDone=False)
    # GitHub
    params.highRiskDose = maxCh
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

    # # Caclulating Percent Yield
    # LH_perc_yield_rec = (simLH_rec[2]/2) / dfyRec * 100
    # LH_perc_yield_github = simLH_github.yields / dfyGithub * 100
    # HL_perc_yield_rec = (simHL_rec[2]/2) / dfyRec * 100
    # HL_perc_yield_github = simHL_github.yields / dfyGithub * 100
    # Mix_perc_yield_rec = (simMix_rec[2]/2) / dfyRec * 100
    # Mix_perc_yield_github = simMix_github.yields / dfyGithub * 100

    # # Calculating Lifetime Yield
    # LH_lty_rec_temp = calcLifetimeYield(LH_perc_yield_rec,(simLH_rec[2]/2),dfyRec)
    # HL_lty_rec_temp = calcLifetimeYield(HL_perc_yield_rec,(simHL_rec[2]/2),dfyRec)
    # Mix_lty_rec_temp = calcLifetimeYield(LH_perc_yield_rec,(simLH_rec[2]/2),dfyRec)
    # LH_lty_github_temp = calcLifetimeYield(LH_perc_yield_github,simLH_github.yields,dfyGithub)
    # HL_lty_github_temp = calcLifetimeYield(HL_perc_yield_github,simHL_github.yields,dfyGithub)
    # Mix_lty_github_temp = calcLifetimeYield(LH_perc_yield_github,simMix_github.yields,dfyGithub)

    # Calculating Lifetime Yield
    LH_lty_rec_temp = calcLifetimeYield((simLH_rec[2]/2),dfyRec)
    HL_lty_rec_temp = calcLifetimeYield((simHL_rec[2]/2),dfyRec)
    Mix_lty_rec_temp = calcLifetimeYield((simMix_rec[2]/2),dfyRec)
    LH_lty_github_temp = calcLifetimeYield(simLH_github.yields,dfyGithub)
    HL_lty_github_temp = calcLifetimeYield(simHL_github.yields,dfyGithub)
    Mix_lty_github_temp = calcLifetimeYield(simMix_github.yields,dfyGithub)

    # Adding Value to List
    ly_altLH_rec = np.append(ly_altLH_rec, LH_lty_rec_temp)
    ly_altHL_rec = np.append(ly_altHL_rec, HL_lty_rec_temp)
    ly_mix_rec = np.append(ly_mix_rec, Mix_lty_rec_temp)
    ly_altLH_github = np.append(ly_altLH_github, LH_lty_github_temp)
    ly_altHL_github = np.append(ly_altHL_github, HL_lty_github_temp)
    ly_mix_github = np.append(ly_mix_github, Mix_lty_github_temp)


# %% Plotting
""" Plotting Lifetime Yield for GitHub Model and Recreation """

fig, axes = plt.subplots(1, 2, figsize=(18, 6))

axes[0].grid()
axes[0].set_title("Github")
axes[0].plot(doseHR, ly_mix_github, 'r', label='mixture')
axes[0].plot(doseHR, ly_altHL_github, 'g', label='Alt High Low')
axes[0].plot(doseHR, ly_altLH_github, 'b', label='Alt Low High')
axes[0].set_xlabel('High-Risk Dose')
axes[0].set_ylabel('Lifetime Yield')
axes[0].legend()

axes[1].grid()
axes[1].set_title("Recreation")
axes[1].plot(doseHR, ly_mix_rec, 'r', label='mixture')
axes[1].plot(doseHR, ly_altHL_rec, 'g', label='Alt High Low')
axes[1].plot(doseHR, ly_altLH_rec, 'b', label='Alt Low High')
axes[1].set_xlabel('High-Risk Dose')
axes[1].set_ylabel('Lifetime Yield')
axes[1].legend()

# Add overall title for the entire figure
fig.suptitle('Lifetime Yield', fontsize=16)

# Adjust layout
plt.tight_layout()
plt.show()

# %%
