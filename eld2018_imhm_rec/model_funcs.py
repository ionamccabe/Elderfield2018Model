import numpy as np
from scipy import integrate
from params_funcs import *
from control_funcs import *

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

def elderfieldOdeSystem_df(time,X):
    S_df, R_df = X
    dS_df = plantGrowth(time,S_df) - plantSenescence(time)*S_df
    dR_df = plantSenescence(time)*S_df
    return np.array([dS_df, dR_df])
    
def runChunk(X,time,controlFunction,wantYield, maxCl, maxCh): # integrates over t and adds results to existing array
    S,Er,Es,Ir,Is,R,Pr,Ps,Ch,Cl,A = X
    Cl[-1],Ch[-1],controlName = controlFunction(time[0],Cl,Ch,maxCl,maxCh)
    X_startcon = S[-1],Er[-1],Es[-1],Ir[-1],Is[-1],R[-1],Pr[-1],Ps[-1],Ch[-1],Cl[-1],A[-1]
    for i, v in enumerate(X_startcon):
        if v < 0:
            X_startcon[i] = 0 # makes sure values are not negative
    
    soln = integrate.solve_ivp(fun = elderfieldOdeSystem, y0 = X_startcon, t_span = time, t_eval = np.arange(time[0],time[1],1)) # solve system of equations
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
    
    ### Calculating Yield
    if wantYield == True:
        dt = np.diff(soln.t)  # Compute time step sizes
        totalYield = np.sum((S_t[:-1] + Er_t[:-1] + Es_t[:-1]) * dt)  # integral approximation
        dfYield = calcYield_df(time)
        percYield = (totalYield / dfYield) * 100
        # print(totalYield,percYield,dfYield)
    else:
        percYield = None
        totalYield = None
        dfYield = None
    ###
    
    return X_ret, percYield, totalYield, dfYield, controlName

def systemWithControl(numSeasons,controlFunction, maxCl, maxCh, printDone = True, terminateEarly = False):
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
    
    
    for num in np.arange(1,numSeasons+1,1): 
        # print('Season', num, 'of', numSeasons)
        time_temp = np.arange(1,t_Emerge,1)
        
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
        
        time_temp = [t_Emerge,t_GS32] # chunk one
        X,pY,tY,dfY,controlName = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh) # solve system of equations
        
        time_temp = [t_GS32,t_GS39] # chunk two
        X,pY,tY,dfY,controlName = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
        
        time_temp = [t_GS39,t_GS61] # chunk three
        X,pY,tY,dfY,controlName = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
        
        time_temp = [t_GS61,t_GS87] # chunk four
        X,pY4,tY4,dfY4,controlName = runChunk(X,time_temp,controlFunction,True, maxCl, maxCh)
        
        ### Adding Yield to array
        percYield_array = np.append(percYield_array,pY4)
        totalYield_array = np.append(totalYield_array,tY4)
        dfYield = dfY4
        ###
        
        time_temp = [t_GS87,seasonLength+1] # chunk four
        X,pY,tY,dfY,controlName = runChunk(X,time_temp,controlFunction,False, maxCl, maxCh)
        
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

    if printDone == True:
        print('Done: ', numSeasonsYield, 'Seasons with', controlName, 'control')
    
    return X, percYield_array, totalYield_array, dfYield, psi_array, SR, controlName

def calcYield_df(time):
    # integrating up until the start of the time
    y0_df = S0,R0 # starting conditions
    time_df_chunk1 = [1,time[0]]
    soln_df_chunk1 = integrate.solve_ivp(fun = elderfieldOdeSystem_df, y0 = y0_df, t_span = time_df_chunk1, t_eval = np.arange(time_df_chunk1[0],time_df_chunk1[1],1))
    y_df_chunk1 = soln_df_chunk1.y
    S_df_chunk1,R_df_chunk1 = y_df_chunk1
    
    # integrating the part we are calculating
    y0_df_chunk2 = S_df_chunk1[-1],R_df_chunk1[-1]
    time_df_chunk2 = [time[0],time[1]]
    soln_df_chunk2 = integrate.solve_ivp(fun = elderfieldOdeSystem_df, y0 = y0_df_chunk2, t_span = time_df_chunk2, t_eval = np.arange(time_df_chunk2[0],time_df_chunk2[1],1))
    y_df_chunk2 = soln_df_chunk2.y
    S_df_chunk2,R_df_chunk2 = y_df_chunk2
    
    # # calulating yield over time frame
    t_df_chunk2 = soln_df_chunk2.t
    dt = np.diff(t_df_chunk2)
    yield_df = np.sum(S_df_chunk2[:-1] * dt)  # Discrete integration
    
    return yield_df

def calcLifetimeYield(percYield_arr, totalYield_arr, dfYield):
    tY_sum = 0
    j = 0
    for i in np.arange(0,len(percYield_arr),1):
        if percYield_arr[i] >= 95:
            tY_sum += totalYield_arr[i]
            j = i + 1
            
        elif i == j and j != 1 and j != 0: ####################### MODIFIED to include the season where failure occurs
            
            tY_sum += totalYield_arr[i]
            j = i+1
    ltY = tY_sum / dfYield # defining lifetime yield in terms of multiples of df yield
    return ltY
            