#%% Fig 4 plots C, D and E
import numpy as np
from model_funcs import systemWithControl
from control_funcs import *
from tqdm import tqdm
from yield_funcs import *
import matplotlib.pyplot as plt

doseHR = np.arange(0,1.025,0.025)
fsy_altLH = np.array([]); fsy_altHL = np.array([]); fsy_mix = np.array([])
SR_altLH = np.array([]); SR_altHL = np.array([]); SR_mix = np.array([])
maxCl = 1
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 1
    
    altLH = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)
    altHL = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)
    mix = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = doseHighRisk, printDone = False)

    X_altLH, percYield_altLH, totalYield_altLH, dfYield_altLH, Psi_altLH, SR_altLH_temp, controlName_altLH = altLH
    X_altHL, percYield_altHL, totalYield_altHL, dfYield_altHL, Psi_altHL, SR_altHL_temp, controlName_altHL = altHL
    X_mix, percYield_mix, totalYield_mix, dfYield_mix, Psi_mix, SR_mix_temp, controlName_mix = mix
    
    ### Plot C
    fsy_altLH_temp = totalYield_altLH[0] / dfYield_altLH
    fsy_altHL_temp = totalYield_altHL[0] / dfYield_altHL
    fsy_mix_temp = totalYield_mix[0] / dfYield_mix
    
    fsy_altLH = np.append(fsy_altLH, fsy_altLH_temp)
    fsy_altHL = np.append(fsy_altHL, fsy_altHL_temp)
    fsy_mix = np.append(fsy_mix, fsy_mix_temp)
    ###
    
    ### Plot D
    SR_altLH = np.append(SR_altLH, SR_altLH_temp)
    SR_altHL = np.append(SR_altHL, SR_altHL_temp)
    SR_mix = np.append(SR_mix, SR_mix_temp)
    ###

# Plot E
ly_altLH = np.array([]); ly_altHL = np.array([]); ly_mix = np.array([])
for i in tqdm(doseHR):
    doseHighRisk = i
    numS = 30
    
    altLH = systemWithControl(numSeasons = numS, controlFunction = altLowHigh, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    altHL = systemWithControl(numSeasons = numS, controlFunction = altHighLow, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    mix = systemWithControl(numSeasons = numS, controlFunction = mixture, maxCl = maxCl, maxCh = doseHighRisk, printDone = False, terminateEarly = False)
    
    X_altLH, percYield_altLH, totalYield_altLH, dfYield_altLH, Psi_altLH, SR_altLH_temp, controlName_altLH = altLH
    X_altHL, percYield_altHL, totalYield_altHL, dfYield_altHL, Psi_altHL, SR_altHL_temp, controlName_altHL = altHL
    X_mix, percYield_mix, totalYield_mix, dfYield_mix, Psi_mix, SR_mix_temp, controlName_mix = mix
    
    ly_altLH_temp = calcLifetimeYield(percYield_altLH, totalYield_altLH, dfYield_altLH)
    ly_altHL_temp = calcLifetimeYield(percYield_altHL, totalYield_altHL, dfYield_altHL)
    ly_mix_temp = calcLifetimeYield(percYield_mix, totalYield_mix, dfYield_mix)
    
    ly_altLH = np.append(ly_altLH, ly_altLH_temp)
    ly_altHL = np.append(ly_altHL, ly_altHL_temp)
    ly_mix = np.append(ly_mix, ly_mix_temp)
###

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot 1: First Season Yield
axes[0].grid()
axes[0].set_title("First Season Yield")
axes[0].plot(doseHR, fsy_mix, 'r', label='Mixture')
axes[0].plot(doseHR, fsy_altHL, 'g', label='Alt High Low')
axes[0].plot(doseHR, fsy_altLH, 'b', label='Alt Low High')
axes[0].set_xlabel('High-Risk Dose')
axes[0].set_ylabel('First Season Yield (%)')
axes[0].legend()

# Plot 2: Selection Ratio (swapped with Lifetime Yield)
axes[1].grid()
axes[1].set_title("Selection Ratio")
axes[1].plot(doseHR, SR_mix, 'r', label='Mixture')
axes[1].plot(doseHR, SR_altHL, 'g', label='Alt High Low')
axes[1].plot(doseHR, SR_altLH, 'b', label='Alt Low High')
axes[1].set_xlabel('High-Risk Dose')
axes[1].set_ylabel('Selection Ratio')
axes[1].legend()

# Plot 3: Lifetime Yield (swapped with Selection Ratio)
axes[2].grid()
axes[2].set_title("Lifetime Yield")
axes[2].plot(doseHR, ly_mix, 'r', label='Mixture')
axes[2].plot(doseHR, ly_altHL, 'g', label='Alt High Low')
axes[2].plot(doseHR, ly_altLH, 'b', label='Alt Low High')
axes[2].set_xlabel('High-Risk Dose')
axes[2].set_ylabel('Lifetime Yield')
axes[2].legend()

# Adjust layout
plt.tight_layout()
plt.show()

#%% Fig 4 Plots A and B
doseHR = np.arange(0.6,1.05,0.05)
doseLR = np.arange(0.25,1.05,0.05)

# Prepare a grid to store Z_SR and Z_LY values for plotting
Z_SR_grid = np.zeros((len(doseLR), len(doseHR)))  # Initialize an empty grid for Z_SR
# Loop through doseLR and doseHR to populate the grids, plot A
for j, doseLowRisk in enumerate(tqdm(doseLR)):
    doseLowRisk = doseLR[j]
    for i, doseHighRisk in enumerate(doseHR):
        doseHighRisk = doseHR[i]
        numS = 1
        
        # Run system simulations with different control functions
        altLH = systemWithControl(numSeasons=numS, controlFunction=altLowHigh, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False)
        altHL = systemWithControl(numSeasons=numS, controlFunction=altHighLow, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False)
        mix = systemWithControl(numSeasons=numS, controlFunction=mixture, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False)
    
        # Extract data from systemWithControl function
        X_altLH, percYield_altLH, totalYield_altLH, dfYield_altLH, Psi_altLH, SR_altLH_temp, controlName_altLH = altLH
        X_altHL, percYield_altHL, totalYield_altHL, dfYield_altHL, Psi_altHL, SR_altHL_temp, controlName_altHL = altHL
        X_mix, percYield_mix, totalYield_mix, dfYield_mix, Psi_mix, SR_mix_temp, controlName_mix = mix
        
        # Plot A (for SR values)
        mixSR = SR_mix_temp
        altSR = min(SR_altHL_temp, SR_altLH_temp)
        
        Z_SR_temp = altSR / (altSR + mixSR)
        Z_SR_grid[j, i] = Z_SR_temp  # Store in grid

Z_LY_grid6 = np.zeros((len(doseLR), len(doseHR)))  # Initialize an empty grid for Z_LY
# Loop through doseLR and doseHR to populate the grids, plot B
for j, doseLowRisk in enumerate(tqdm(doseLR)):
    doseLowRisk = doseLR[j]
    for i, doseHighRisk in enumerate(doseHR):
        doseHighRisk = doseHR[i]
        numS = 25
        
        # print(doseLowRisk, doseHighRisk)
        
        # Run system simulations with different control functions
        altLH = systemWithControl(numSeasons=numS, controlFunction=altLowHigh, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False, terminateEarly = True)
        altHL = systemWithControl(numSeasons=numS, controlFunction=altHighLow, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False, terminateEarly = True)
        mix = systemWithControl(numSeasons=numS, controlFunction=mixture, maxCl=doseLowRisk, maxCh=doseHighRisk, printDone=False, terminateEarly = True)
    
        # Extract data from systemWithControl function
        X_altLH, percYield_altLH, totalYield_altLH, dfYield_altLH, Psi_altLH, SR_altLH_temp, controlName_altLH = altLH
        X_altHL, percYield_altHL, totalYield_altHL, dfYield_altHL, Psi_altHL, SR_altHL_temp, controlName_altHL = altHL
        X_mix, percYield_mix, totalYield_mix, dfYield_mix, Psi_mix, SR_mix_temp, controlName_mix = mix

        # Plot B (for LY values)
        ly_altLH_temp = calcLifetimeYield(percYield_altLH, totalYield_altLH, dfYield_altLH)
        ly_altHL_temp = calcLifetimeYield(percYield_altHL, totalYield_altHL, dfYield_altHL)
        ly_mix_temp = calcLifetimeYield(percYield_mix, totalYield_mix, dfYield_mix)
        
        mixLY = ly_mix_temp
        altLY = max(ly_altHL_temp, ly_altLH_temp)
        
        Z_LY_temp = mixLY / (altLY + mixLY)
        Z_LY_grid6[j, i] = Z_LY_temp  # Store in grid

# Plotting the contour plots for SR and LY
# Create a meshgrid for plotting
X, Y = np.meshgrid(doseHR, doseLR)

# Define custom colors for the intervals
colors = ['#4D330F', '#9C6F2C', '#E0CD95', '#FFFFFF', '#AAD5CC', '#548780', '#233E35']  # Updated colors

# Create a custom colormap
cmap = ListedColormap(colors)

# Define updated contour levels, adding 0.45 and 0.55
levels = [0.0, 0.1, 0.3, 0.498, 0.502, 0.7, 0.9, 1.0]  # 8 boundaries for 7 intervals

# Compute midpoints of intervals for labeling
midpoints = [(levels[i] + levels[i+1]) / 2 for i in range(len(levels) - 1)]

### Plot for SR
plt.figure(figsize=(10, 6))
cp = plt.contourf(Y, X, Z_SR_grid, levels=levels, cmap=cmap, extend='both')

# Add color bar
cbar = plt.colorbar(cp, ticks=midpoints)  # Use midpoints for ticks
cbar.ax.set_yticklabels(['0.0-0.1', '0.1-0.3', '0.3-0.498', '0.498-0.502', '0.502-0.7', '0.7-0.9', '0.9-1.0'])  # Updated labels

# Add labels and title
plt.xlabel('Low Risk Dose')
plt.ylabel('High Risk Dose')
plt.title('Z-Metric for Selection Ratio')
plt.show()

### Plot for LY
plt.figure(figsize=(10, 6))
cp = plt.contourf(Y, X, Z_LY_grid6, levels=levels, cmap=cmap, extend='both')

# Add color bar
cbar = plt.colorbar(cp, ticks=midpoints)  # Use midpoints for ticks
cbar.ax.set_yticklabels(['0.0-0.1', '0.1-0.3', '0.3-0.498', '0.498-0.502', '0.502-0.7', '0.7-0.9', '0.9-1.0'])  # Updated labels

# Add labels and title
plt.xlabel('Low Risk Dose')
plt.ylabel('High Risk Dose')
plt.title('Z-Metric for Lifetime Yield')
plt.show()
# %%
