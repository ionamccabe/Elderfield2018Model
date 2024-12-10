import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model\eld2018_imhm_rec")

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from params_funcs import *
from yield_funcs import calcLifetimeYield

def plotDynamics(numSeasons,sWC_output):
    
    X, percYield, totalYield, dfYield, Psi, SR, controlName, time = sWC_output
    
    S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A = X
    
    t = np.arange(0,seasonLength*numSeasons,1) # t for plotting
    
    season_ticks = np.arange(0, seasonLength * numSeasons+1, seasonLength) # make x axis seasons
    season_labels = np.arange(0,numSeasons+1,1)
    
    # Define important time points and labels
    important_points = [1212, 1456, 1700, 2066, 2900]
    important_labels = ['Emerge', 'GS32', 'GS39', 'GS61', 'GS87']
    
    # Create a 2x2 grid of subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 10))
    
    # Function to add vertical dashed lines and labels
    def add_vertical_lines(ax,numSeasons):
        if numSeasons == 1:
            for i, point in enumerate(important_points):
                ax.axvline(x=point, color='gray', linestyle='--', linewidth=1)
                ax.text(point, ax.get_ylim()[1] * 0.9,  # Position label near the top
                        important_labels[i],
                        rotation=90, va='top', ha='center', fontsize=9, color='gray')
        else:
            None
    
    # Plot 1: Susceptible and Removed
    axs[0, 0].grid()
    axs[0, 0].set_title("Susceptible and Removed")
    axs[0, 0].plot(t, S, 'orange', label='Susceptible')
    axs[0, 0].plot(t, R, 'b', label='Removed')
    add_vertical_lines(axs[0, 0],numSeasons)
    axs[0, 0].set_xticks(season_ticks)
    axs[0, 0].set_xticklabels(season_labels)
    axs[0, 0].set_xlabel('Seasons')
    axs[0, 0].set_ylabel('LAI')
    axs[0, 0].legend()
    
    # Plot 2: Infections
    axs[0, 1].grid()
    axs[0, 1].set_title("Infections")
    axs[0, 1].plot(t, Er, 'r', label='Latent: Resistant')
    axs[0, 1].plot(t, Es, 'g', label='Latent: Sensitive')
    axs[0, 1].plot(t, Ir, 'pink', label='Infected: Resistant')
    axs[0, 1].plot(t, Is, 'b', label='Infected: Sensitive')
    add_vertical_lines(axs[0, 1],numSeasons)
    axs[0, 1].set_xticks(season_ticks)
    axs[0, 1].set_xticklabels(season_labels)
    axs[0, 1].set_xlabel('Seasons')
    axs[0, 1].set_ylabel('LAI')
    axs[0, 1].legend()
    
    # Plot 3: Fungicide Concentration
    axs[1, 0].grid()
    axs[1, 0].set_title("Fungicide Concentration")
    axs[1, 0].plot(t, Ch, 'r', label='High-Risk Fungicide')
    axs[1, 0].plot(t, Cl, 'b', label='Low-Risk Fungicide')
    add_vertical_lines(axs[1, 0],numSeasons)
    axs[1, 0].set_xticks(season_ticks)
    axs[1, 0].set_xticklabels(season_labels)
    axs[1, 0].set_xlabel('Seasons')
    axs[1, 0].set_ylabel('LAI')
    axs[1, 0].legend()
    
    # Plot 4: Initial Inoculum
    axs[1, 1].grid()
    axs[1, 1].set_title("Initial Inoculum")
    axs[1, 1].plot(t, Ps, 'g', label='Primary Sensitive')
    axs[1, 1].plot(t, Pr, 'b', label='Primary Resistant')
    add_vertical_lines(axs[1, 1],numSeasons)
    axs[1, 1].set_xticks(season_ticks)
    axs[1, 1].set_xticklabels(season_labels)
    axs[1, 1].set_xlabel('Seasons')
    axs[1, 1].set_ylabel('LAI')
    axs[1, 1].legend()
    
    # Adjust layout and display the plot
    plt.tight_layout()
    fig.suptitle(controlName)
    fig.subplots_adjust(top=0.935)
    plt.show()

def plotComparison(numSeasons, sWC_output1, sWC_output2, sWC_output3):
    X1, percYield1, totalYield1, dfYield1, Psi1, SR1, controlName1, t1 = sWC_output1
    X2, percYield2, totalYield2, dfYield2, Psi2, SR2, controlName2, t2 = sWC_output2
    X3, percYield3, totalYield3, dfYield3, Psi3, SR3, controlName3, t3 = sWC_output3
    
    tYield = np.arange(0,numSeasons,1) # for 
    
    
    ### calculating lifetime yield
    ltY1 = calcLifetimeYield(percYield1, totalYield1, dfYield1)
    ltY2 = calcLifetimeYield(percYield2, totalYield2, dfYield2)
    ltY3 = calcLifetimeYield(percYield3, totalYield3, dfYield3)
    # making stuff for bar plot
    treatments = [controlName1,controlName2,controlName3]
    values = [ltY1,ltY2,ltY3]
    

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Yield
    axs[0].grid()
    axs[0].set_title("Yield")
    axs[0].plot(tYield, percYield1, 'r', label=controlName1)
    axs[0].plot(tYield, percYield2, 'g', label=controlName2)
    axs[0].plot(tYield, percYield3, 'b', label=controlName3)
    axs[0].set_xlabel('Seasons')
    axs[0].set_ylabel('Yield (% of disease-free)')
    axs[0].legend()
    
    # Plot 2: Resistance Frequency
    axs[1].grid()
    axs[1].set_title("Resistance Frequency")
    axs[1].plot(tYield, Psi1, 'r', label=controlName1)
    axs[1].plot(tYield, Psi2, 'g', label=controlName2)
    axs[1].plot(tYield, Psi3, 'b', label=controlName3)
    axs[1].set_xlabel('Seasons')
    axs[1].set_ylabel('Resistance Frequency')
    axs[1].legend()
    
    # Plot 3: Lifetime Yield as a Bar Chart
    axs[2].bar(treatments, values, color=['r', 'g', 'b'], width=0.4)
    axs[2].set_title("Lifetime Yield")
    axs[2].set_ylabel("Lifetime Yield")
    
    # Adjust layout and show the plot
    plt.tight_layout()
    plt.show()
