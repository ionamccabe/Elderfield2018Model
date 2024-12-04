
# %% 
# 1 Year
from model_funcs import *
from control_funcs import *
from plot_funcs import *

# running the simulation
numS = 1
maxCl = 1; maxCh = 1
altLH = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh)
altHL = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh)
mix = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh)

# plotting dynamics
plotDynamics(numS, altLH)
plotDynamics(numS, altHL)
plotDynamics(numS, mix)

# 30 years
numS = 30
altLH = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = maxCh)
altHL = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = maxCh)
mix = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = maxCh)

# plotting comparison
plotComparison(numS,mix,altHL,altLH)

# %%
