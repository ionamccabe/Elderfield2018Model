import math

# transmission parameters
beta = 1.56*10**(-2) # infection rate param, degree-day^-1
gamma = 1/266 # 1/latent period, degree-day^-1
mu = 1/456 # # 1/infectious period, degree-day^-1

# fungicide paramenters
omega_h = 1; omega_l = 0.48 # fung max effect
theta_h = 9.6; theta_l = 9.9 # fung response curve
delta_h = 1.11*10**(-2); delta_l = 6.91*10**(-3) # fung decay param, degree-day^-1
Cl0 = 0; Ch0 = 0

# host parameters
r = 1.26 * 10**(-2) # host growth rate, degree-days^-1
k = 4.2 # max leaf area, unitless

# pathogen parameters
psi = 1*10**(-10) # initial resistance frequency
phi = 1.09*10**(-2) # initial incoculum density
v = 8.5*10**(-3) # decay rate of primary innoc, degree-days^-1

# initializing the simulation
S0 = 0.05; Er0 = 0; Es0 = 0; Ir0 = 0; Is0 = 0; R0 = 0 # initial populations
A0 = S0 + Er0 + Es0 + Ir0 + Is0 + R0
Pr0 = psi*phi # initial area of lower leaf infected by resistant strain
Ps0 = (1-psi)*phi # initial area of lower leaf infected by sensitive strain


# time points, degree days
t_Emerge = 1212; t_GS32 = 1456; t_GS39 = 1700; t_GS61 = 2066; t_GS87 = 2900

# length of a season
seasonLength = 3000 # degree-days

#%% functions

def fungicideEffect(concentration, omega, theta):
    fe = omega*(1-math.exp(-theta*concentration))
    return(fe)

def plantGrowth(time,totalTissueA):
    yearsSinceStart = time//seasonLength
    timeOfYear = time - yearsSinceStart*seasonLength
    if t_Emerge < timeOfYear < t_GS61: # assumed to stop growing when senescence begins
        g = r*(k-totalTissueA)
    else:
        g = 0
    return(g)

def plantSenescence(time):
    yearsSinceStart = time//seasonLength
    timeOfYear = time - yearsSinceStart*seasonLength
    if timeOfYear < t_GS61:
        s = 0
    else:
        s = 0.005*(timeOfYear-t_GS61)/(t_GS87-t_GS61) + 0.1*math.exp(-0.02*(t_GS87 - timeOfYear))
    return(s)