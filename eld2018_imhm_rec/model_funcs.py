import sys
sys.path.append(r"C:\Users\ionac\Documents\python\eld2018\Elderfield2018Model\eld2018_imhm_rec")

import numpy as np
from scipy import integrate
from scipy.integrate import simpson
from params_funcs import *
from control_funcs import *
from yield_funcs import calcYield_df

def elderfieldOdeSystem(time,X):
    S, Er, Es, Ir, Is, R, Pr, Ps, Ch, Cl, A = X
    dS = plantGrowth(time,A) - plantSenescence(time)*S - beta*S/A*(1-fungicideEffect(Cl, omega_l, theta_l))*((1-fungicideEffect(Ch, omega_h, theta_h))*(Is + Ps) + Ir + Pr)
    dEr = beta*S/A*(1-fungicideEffect(Cl, omega_l, theta_l))*(Ir + Pr) - plantSenescence(time)*Er - gamma*Er
    dEs = beta*S/A*(1-fungicideEffect(Cl, omega_l, theta_l))*(1-fungicideEffect(Ch, omega_h, theta_h))*(Is + Ps) - plantSenescence(time)*Es - gamma*(1-fungicideEffect(Ch, omega_h, theta_h))*Es
    dIr = gamma*Er - mu*Ir
    dIs = gamma*(1-fungicideEffect(Ch, omega_h, theta_h))*Es - mu*Is
    dR = mu*(Ir + Is) + plantSenescence(time)*(S + Er + Es)
    dPr = -v*Pr
    dPs = -v*Ps
    dCh = -delta_h*Ch
    dCl = -delta_l*Cl
    dA = dS + dEr + dEs + dIr + dIs + dR
    return np.array([dS, dEr, dEs, dIr, dIs, dR, dPr, dPs, dCh, dCl, dA])

# James' verion of above is Simulation.WheatSeptoria


def runChunk(X,time,controlFunction,wantYield, maxCl, maxCh): # integrates over t and adds results to existing array
    S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A = X
    Cl[-1],Ch[-1],controlName = controlFunction(time[0],Cl,Ch,maxCl,maxCh)
    X_startcon = S[-1],Er[-1],Es[-1],Ir[-1],Is[-1],R[-1],Pr[-1],Ps[-1],Ch[-1],Cl[-1],A[-1]
    for i, v in enumerate(X_startcon):
        if v < 0:
            X_startcon[i] = 0 # makes sure values are not negative
    
    soln = integrate.solve_ivp(fun = elderfieldOdeSystem, y0 = X_startcon, t_span = time, t_eval = np.arange(time[0],time[1],1), method = 'LSODA') # solve system of equations
    y = soln.y # extract vectors from the solution
    S_t, Er_t, Es_t, Ir_t, Is_t, R_t, Pr_t, Ps_t, Ch_t, Cl_t, A_t = y # split into vectors
    S = np.append(S,S_t)
    Er = np.append(Er,Er_t)
    Es = np.append(Es,Es_t)
    Ir = np.append(Ir,Ir_t)
    Is = np.append(Is,Is_t)
    R = np.append(R,R_t)
    Pr = np.append(Pr,Pr_t)
    Ps = np.append(Ps,Ps_t)
    Ch = np.append(Ch,Ch_t)
    Cl = np.append(Cl,Cl_t)
    A = np.append(A,A_t)
    X_ret = S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A

    soln_time = soln.t
    
    ### Calculating Yield
    if wantYield == True:
        yield_metrics = S_t + Er_t + Es_t
        totalYield = np.trapezoid(y = yield_metrics, x = soln.t)  # integral approximation
        dfYield = calcYield_df(time)
        percYield = (totalYield / dfYield) * 100
        # print(totalYield,percYield,dfYield)
    else:
        percYield = None
        totalYield = None
        dfYield = None
    ###
    
    return X_ret, percYield, totalYield, dfYield, controlName, soln_time

def systemWithControl(numSeasons,controlFunction, maxCl, maxCh, printDone = True, terminateEarly = False, includeStart = False):
    # year should be divided into five chunks:
    # t_Emerge - t_GS32, t_GS32 - t_GS_39, t_GS39 - t_GS61, t_GS61 - t_GS87, t_GS87 - onwards
    
    S = np.array([]); Er = np.array([]); Es = np.array([])
    Ir = np.array([]); Is = np.array([]); R = np.array([]); 
    Pr = np.array([]); Ps = np.array([])
    Ch = np.array([]); Cl = np.array([]); A = np.array([]) # initializing
    
    psi_array = np.array([])
    percYield_array = np.array([])
    totalYield_array = np.array([])
  
    psi_temp = psi  
    numSeasonsYield = 0
    time_arr = np.array([])
    
    
    for num in np.arange(1,numSeasons+1,1): 
        # print('Season', num, 'of', numSeasons)
        time_arr_temp = np.array([])
        if includeStart == True:
            time_temp = np.arange(1,t_Emerge,1)
            soln_time = np.arange(time_temp[0],time_temp[1],1)
            time_arr_temp = np.append(time_arr_temp,soln_time)
            initial_con_time = seasonLength + 1
        else:
            time_temp = [0]
            time_arr_temp = np.append(time_arr_temp,1211)
            initial_con_time = 1689 + 1 # 1 degree day after the end of the last season
            
        for t in time_temp:
            Pr_temp = psi_temp*phi # updating initial source of innoculum based on psi
            Ps_temp = (1-psi_temp)*phi
            
            S = np.append(S,S0)
            Er = np.append(Er,Er0)
            Es = np.append(Es,Es0)
            Ir = np.append(Ir,Ir0)
            Is = np.append(Is,Is0)
            R = np.append(R,R0)
            Pr = np.append(Pr,Pr_temp)
            Ps = np.append(Ps,Ps_temp)
            Ch = np.append(Ch,Ch0)
            Cl = np.append(Cl,Cl0)
            A = np.append(A,A0)
        
        X = S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A
        
        time_temp = [t_Emerge,t_GS32] # chunk one, already did t = t_Emerge
        X,pY,tY,dfY,controlName,soln_time = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh) # solve system of equations
        time_arr_temp = np.append(time_arr_temp,soln_time)

        time_temp = [t_GS32,t_GS39] # chunk two
        X,pY,tY,dfY,controlName,soln_time = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
        time_arr_temp = np.append(time_arr_temp,soln_time)
        
        time_temp = [t_GS39,t_GS61] # chunk three
        X,pY,tY,dfY,controlName,soln_time = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
        time_arr_temp = np.append(time_arr_temp,soln_time)

        if includeStart == True:
            endTime = t_GS87
        else:
            endTime = t_GS87 + 1

        time_temp = [t_GS61,endTime] # chunk four
        X,pY4,tY4,dfY4,controlName,soln_time = runChunk(X,time_temp,controlFunction,True, maxCl, maxCh)
        time_arr_temp = np.append(time_arr_temp,soln_time)
        
        ### Adding Yield to array
        percYield_array = np.append(percYield_array,pY4)
        totalYield_array = np.append(totalYield_array,tY4)
        dfYield = dfY4
        ###
        
        if includeStart == True:
            time_temp = [t_GS87,seasonLength+1] # chunk four
            X,pY,tY,dfY,controlName,soln_time = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
            time_arr_temp = np.append(time_arr_temp,soln_time)

        ### Updating psi
        S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A = X
        psi_temp = Ir[-1]/(Ir[-1] + Is[-1]) # updating psi (prop. of r in pop)
        psi_array = np.append(psi_array,psi_temp)
        ###

        ### Calculating Selection Ratio
        if num == 1:
            initialPropR = Pr0/(Pr0 + Ps0)
            SR = psi_temp/initialPropR
        ###

        numSeasonsYield += 1
        if pY4 < 95 and terminateEarly == True:
            break

        time_arr_temp = time_arr_temp[1:]
        time_arr = np.append(time_arr, time_arr_temp)

    if printDone == True:
        print('Done: ', numSeasonsYield, 'Seasons with', controlName, 'control')
    
    X_arr = np.array(X) 
    # indices_to_remove = [(initial_con_time * i) - 1 for i in range(1, numSeasons + 1)]
    indices_to_remove = [0]
    ict = initial_con_time + 1 # becuase the first element has not been removed yet
    for i in range(1, numSeasons): # removes initial conditions AND last element of each season (to match James)
        # Remove 1689th and 1690th elements for all but the last season
        indices_to_remove.extend([(ict * i) - 2 - (i-1), (ict * i) - 1 - (i-1)])

    # Step 2: Create a mask to retain only the desired columns
    mask = np.ones(X_arr.shape[1], dtype=bool)  # Start with all True
    mask[indices_to_remove] = False  # Set 1690th, 3380th, etc., to False

    # Step 3: Apply the mask to remove columns
    processed_X = X_arr[:, mask]
    # X_arr = np.array([state[1:] for state in X])
    X_tup = tuple(processed_X)

    return X_tup, percYield_array, totalYield_array, dfYield, psi_array, SR, controlName, time_arr
