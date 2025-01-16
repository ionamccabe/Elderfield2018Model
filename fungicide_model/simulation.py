#!usr/bin/python
from .constants import Model, StateIndices, SprayStrategies, FungicideEffectType
import copy
from functools import reduce
from .fungicide import Fungicide
from math import exp
from numpy import arange
from numpy import asarray
from numpy import trapz
from numpy import diff
from pandas import DataFrame
from scipy.integrate import ode
from .stepChange import StepChange
from .stepChange import StepChangeKind
from operator import attrgetter
from .fungicideTreatment import FungicideTreatment
import math


class Simulation(object):
    def __init__(self, params, stopOnYieldFailure = False):
        self.params = copy.deepcopy(params)
        self.selectionRatio = None
        self.yieldTilCriticalLoss = None
        self.yieldStop = stopOnYieldFailure
        self.output = None

    integrator = "lsoda" # setting up integrator to be 'lsoda'
    outputRate = 1
    maxRunTime = 50

    @staticmethod
    def _dictStateFromSolver(solver):
        d = {
            "t": solver.t,
            "S": solver.y[StateIndices.S],
            "ER": solver.y[StateIndices.ER],
            "ES": solver.y[StateIndices.ES],
            "IR": solver.y[StateIndices.IR],
            "IS": solver.y[StateIndices.IS],
            "R": solver.y[StateIndices.R],
            "PR": solver.y[StateIndices.PR],
            "PS": solver.y[StateIndices.PS],
            "high": solver.y[StateIndices.high],
            "low": solver.y[StateIndices.low],
            "O": solver.y[StateIndices.O]
        }
        return d

    def run(self):
        self._calculateDiseaseFreeYield()
        self._run()
        self._postProcess()

    def _run(self):
        if self.params.model == Model.WheatSeptoria:
            f = Simulation.wheatSeptoria
            nextYieldStart = self.params.GS61 - self.params.seasonStartTime
            nextYieldEnd = self.params.GS87 - self.params.seasonStartTime
            seasonLength = self.params.GS87 - self.params.seasonStartTime
        elif self.params.model == Model.GrapePowderyMildew:
            f = Simulation.grapePowderyMildew
            nextYieldStart = self.params.yieldStartTime - self.params.budBreakTime
            nextYieldEnd = self.params.yieldEndTime - self.params.budBreakTime
            seasonLength = self.params.seasonEndTime - self.params.budBreakTime
        else:
            raise RuntimeError("Unrecognised model: " + self.params.model)

        solver = ode(f).set_integrator(Simulation.integrator, nsteps = 1E3) # setting up ode() solver, integrator is 'lsoda' from above
        solver.set_f_params(self.params) # setting parameters

        state = [0 for x in range(11)] # setting all states to 0

        solver.set_initial_value(state, 0) # initial values: states (from above) and t=0

        lastOutput = None

        changes = self.getStepChanges()

        nextStepChangeIndex = 0
        nextStepTime = changes[0].time

        outputList = []
        self.yields = []
        thisYield = 0

        self.params.rFreq = self.params.psi

        stopSimulation = False # what does this do?

        # Needs to run for at least a season for selection ratio calculation
        if self.params.model == Model.WheatSeptoria:
            minimumRunTime = self.params.GS87 - self.params.seasonStartTime
        elif self.params.model == Model.GrapePowderyMildew:
            minimumRunTime = self.params.seasonEndTime - self.params.budBreakTime

        while solver.t < self.params.maxTime: # while the integrated time is less than the max time
            runTil = min( # run until next event, idk, max time, or idk
                nextStepTime,
                float("inf") if lastOutput is None else lastOutput + Simulation.outputRate, # lastOutput is solver.t, outputRate = 1, logging rate
                self.params.maxTime,
                solver.t + Simulation.maxRunTime # 
            )

            if runTil > solver.t: # run until time runTil
                solver.integrate(runTil)

            if not solver.successful():
                raise RuntimeError("Solver failed.")

            newState = solver.y # assigning the solution to be a new state

            for i, v in enumerate(newState): # making sure the states are never negative
                if v < 0:
                    newState[i] = 0

            solver.set_initial_value(newState, solver.t) # setting new solution and end time as initial conditions

            while solver.t == nextStepTime: # when the solver reaches a time step
                solver.set_initial_value(changes[nextStepChangeIndex].apply(solver.y), solver.t) 

                nextStepChangeIndex += 1 # next step change
                if nextStepChangeIndex >= len(changes): # if there are no more step changes, set nextStepTime to infitiy
                    nextStepTime = float("inf")
                else: # otherwise set to the next step time
                    nextStepTime = changes[nextStepChangeIndex].time
            if lastOutput is None or solver.t == lastOutput + Simulation.outputRate: # log the first one, then +1 after (repetedly)
                outputList.append(Simulation._dictStateFromSolver(solver))

                lastOutput = solver.t

                self.params.rFreq = self._resistanceFrequency(outputList[-1])

                if solver.t > nextYieldStart:
                    if solver.t < nextYieldEnd:
                        if self.params.model == Model.WheatSeptoria:
                            yieldCompartments = ["S", "ER", "ES"]

                            if len(outputList) > 1:
                                previousValues = outputList[-2]
                                tDiff = outputList[-1]["t"] - outputList[-2]["t"]
                            else:
                                previousValues = {
                                    "S": self.params.S0() if callable(self.params.S0) else self.params.S0,
                                    "ER": self.params.ER0() if callable(self.params.ER0) else self.params.ER0,
                                    "ES": self.params.ES0() if callable(self.params.ES0) else self.params.ES0,
                                    "O": self.params.O0() if callable(self.params.O0) else self.params.O0
                                }
                                tDiff = outputList[-1]["t"]

                            thisYield += reduce(
                                lambda accum, el: accum + tDiff / 4 * (previousValues[el] + outputList[-1][el]),
                                yieldCompartments,
                                0
                            )
                        else:
                            allCompartments = ["S", "ER", "ES", "IR", "IS", "O"]
                            infectiousCompartments = ["IR", "IS"]
                            totalArea = 0
                            totalInfection = 0

                            for c in allCompartments:
                                totalArea += outputList[-1][c]

                                if c in infectiousCompartments:
                                    totalInfection += outputList[-1][c]

                            severity = totalInfection / totalArea
                            thisYield = max(thisYield, severity)
                    if solver.t >= nextYieldEnd:
                        if self.params.model == Model.GrapePowderyMildew:
                            # For grape it's not really a yield...
                            thisYield = 1 - thisYield

                        if self.params.densityDependentInfection and self.yieldStop:
                            if thisYield < (1 - self.params.criticalYieldLoss) * self.diseaseFreeYield:
                                stopSimulation = True

                        self.yields.append(thisYield)

                        thisYield = 0
                        nextYieldStart += seasonLength
                        nextYieldEnd += seasonLength
            if stopSimulation and solver.t >= minimumRunTime:
                break

        self.output = DataFrame.from_records(outputList).set_index("t")

    def _calculateDiseaseFreeYield(self):
        diseaseFreeParams = copy.deepcopy(self.params)
        diseaseFreeParams.strategy = SprayStrategies.NoSpray

        if diseaseFreeParams.model == Model.WheatSeptoria:
            diseaseFreeParams.beta = 0
        else:
            diseaseFreeParams.beta = [0, 0]
            diseaseFreeParams.phi = 0

        diseaseFreeParams.maxSeasons = 1

        diseaseFreeSim = Simulation(diseaseFreeParams)
        diseaseFreeSim._run()

        self.diseaseFreeYield = diseaseFreeSim.yields[0]

    def getStepChanges(self):
        if isinstance(self.params.strategy, FungicideTreatment):
            changes = self.params.strategy.createStepChanges(self.params)
        elif self.params.model == Model.WheatSeptoria:
            t = 0
            changes = []

            if not self.params.fungicideDecay:
                highRiskLifetime = 1 / self.params.highRisk.decayRate
                lowRiskLifetime = 1 / self.params.lowRisk.decayRate

                maxLifeTime = self.params.GS39 - self.params.GS32
                largestLifeTime = max(highRiskLifetime, lowRiskLifetime)

                if largestLifeTime >= maxLifeTime:
                    multiplier = maxLifeTime / largestLifeTime
                    highRiskLifetime *= multiplier
                    lowRiskLifetime *= multiplier
                    highRiskLifetime = round(highRiskLifetime, 6)
                    lowRiskLifetime = round(lowRiskLifetime, 6)

            while t < self.params.maxTime:
                # Do start of season initialisation
                changes.append(StepChange(
                    StateIndices.S, t, self.params.S0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.ER, t, attrgetter("ER0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.ES, t, attrgetter("ES0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.IR, t, attrgetter("IR0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.IS, t, attrgetter("IS0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.R, t, self.params.R0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.PR, t, attrgetter("PR0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.PS, t, attrgetter("PS0"), StepChangeKind.SET, self.params
                ))

                # Add T1 spray
                t += self.params.GS32 - self.params.seasonStartTime
                if self.params.strategy == SprayStrategies.NoSpray:
                    pass
                elif self.params.strategy == SprayStrategies.Mixture: ### CHANGED: divided dose by 2
                    changes.append(StepChange(
                        StateIndices.high, t, self.params.highRiskDose/2, StepChangeKind.SET, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.low, t, self.params.lowRiskDose/2, StepChangeKind.SET, self.params
                    ))
                elif self.params.strategy == SprayStrategies.AltHiLo:
                    changes.append(StepChange(
                        StateIndices.high, t, self.params.highRiskDose, StepChangeKind.SET, self.params
                    ))
                elif self.params.strategy == SprayStrategies.AltLoHi:
                    changes.append(StepChange(
                        StateIndices.low, t, self.params.lowRiskDose, StepChangeKind.SET, self.params
                    ))

                if not self.params.fungicideDecay:
                    changes.append(StepChange(
                        StateIndices.high, t + highRiskLifetime, 0, StepChangeKind.SET, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.low, t + lowRiskLifetime, 0, StepChangeKind.SET, self.params
                    ))

                # Move time on to T2
                t += self.params.GS39 - self.params.GS32
                if self.params.strategy == SprayStrategies.NoSpray:
                    pass
                elif self.params.strategy == SprayStrategies.Mixture: ### CHANGED: divided dose by 2
                    changes.append(StepChange(
                        StateIndices.high, t, self.params.highRiskDose/2, StepChangeKind.ADD, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.low, t, self.params.lowRiskDose/2, StepChangeKind.ADD, self.params
                    ))
                elif self.params.strategy == SprayStrategies.AltHiLo:
                    changes.append(StepChange(
                        StateIndices.low, t, self.params.lowRiskDose, StepChangeKind.ADD, self.params
                    ))
                elif self.params.strategy == SprayStrategies.AltLoHi:
                    changes.append(StepChange(
                        StateIndices.high, t, self.params.highRiskDose, StepChangeKind.ADD, self.params
                    ))

                if not self.params.fungicideDecay:
                    changes.append(StepChange(
                        StateIndices.high, t + highRiskLifetime, 0, StepChangeKind.SET, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.low, t + lowRiskLifetime, 0, StepChangeKind.SET, self.params
                    ))

                # To next season!
                t += self.params.GS87 - self.params.GS39
        elif self.params.model == Model.GrapePowderyMildew:
            t = 0
            changes = []

            if not self.params.fungicideDecay:
                highRiskLifetime = 1 / self.params.highRisk.decayRate
                lowRiskLifetime = 1 / self.params.lowRisk.decayRate

                maxLifeTime = min(diff(asarray(self.params.sprayTimes)))
                largestLifeTime = max(highRiskLifetime, lowRiskLifetime)

                if largestLifeTime >= maxLifeTime:
                    multiplier = maxLifeTime / largestLifeTime
                    highRiskLifetime *= multiplier
                    lowRiskLifetime *= multiplier
                    highRiskLifetime = round(highRiskLifetime, 6)
                    lowRiskLifetime = round(lowRiskLifetime, 6)

            while t < self.params.maxTime:
                # Initialisation at start of season
                changes.append(StepChange(
                    StateIndices.S, t, self.params.S0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.ER, t, attrgetter("ER0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.ES, t, attrgetter("ES0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.IR, t, attrgetter("IR0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.IS, t, attrgetter("IS0"), StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.R, t, self.params.R0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.O, t, self.params.O0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.high, t, 0, StepChangeKind.SET, self.params
                ))
                changes.append(StepChange(
                    StateIndices.low, t, 0, StepChangeKind.SET, self.params
                ))

                # Apply shoot topping
                t += self.params.shootToppingTime - self.params.budBreakTime

                if self.params.seasonality:
                    changes.append(StepChange(
                        StateIndices.S, t, (1 - self.params.toppingLoss["S"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.ER, t, (1 - self.params.toppingLoss["E"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.ES, t, (1 - self.params.toppingLoss["E"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.IR, t, (1 - self.params.toppingLoss["I"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.IS, t, (1 - self.params.toppingLoss["I"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.R, t, (1 - self.params.toppingLoss["R"]), StepChangeKind.MULTIPLY, self.params
                    ))
                    changes.append(StepChange(
                        StateIndices.O, t, (1 - self.params.toppingLoss["O"]), StepChangeKind.MULTIPLY, self.params
                    ))

                # To next season
                t += self.params.seasonEndTime - self.params.shootToppingTime

            # Do fungicides sprays separately
            t = 0
            if len(self.params.sprayTimes) > 0:
                while t < self.params.maxTime:
                    for i, sprayTime in enumerate(self.params.sprayTimes):
                        thisT = t + sprayTime - self.params.budBreakTime

                        if thisT > self.params.maxTime:
                            break

                        if sprayTime in self.params.alternatingSprayTimes and self.params.strategy != SprayStrategies.NoSpray:
                            if self.params.strategy == SprayStrategies.Mixture:
                                sprayHigh = sprayLow = True
                            elif self.params.strategy == SprayStrategies.AltHiLo:
                                if i % 2 == 0:
                                    sprayHigh = True
                                    sprayLow = False
                                else:
                                    sprayHigh = False
                                    sprayLow = True
                            elif self.params.strategy == SprayStrategies.AltLoHi:
                                if i % 2 == 0:
                                    sprayHigh = False
                                    sprayLow = True
                                else:
                                    sprayHigh = True
                                    sprayLow = False
                            else:
                                sprayHigh = sprayLow = False
                            highDose = self.params.highRiskDose
                            lowDose = self.params.lowRiskDose
                        else:
                            if self.params.strategy == SprayStrategies.NoSpray:
                                sprayLow = sprayHigh = False
                            else:
                                sprayLow = True
                                sprayHigh = False

                            # Always apply these at full dose
                            lowDose = 1

                        if sprayHigh:
                            changes.append(StepChange(
                                StateIndices.high, thisT, highDose, StepChangeKind.ADD, self.params
                            ))
                            if not self.params.fungicideDecay:
                                changes.append(StepChange(
                                    StateIndices.high, thisT + highRiskLifetime, 0, StepChangeKind.SET, self.params
                                ))

                        if sprayLow:
                            changes.append(StepChange(
                                StateIndices.low, thisT, lowDose, StepChangeKind.ADD, self.params
                            ))
                            if not self.params.fungicideDecay:
                                changes.append(StepChange(
                                    StateIndices.low, thisT + lowRiskLifetime, 0, StepChangeKind.SET, self.params
                                ))

                    t += self.params.seasonEndTime - self.params.budBreakTime
        else:
            raise RuntimeError("Unknown model: " + params.model)

        # Make sure that the changes are happening in the correct order
        changes.sort(key = lambda x: x.time)

        return changes

    @staticmethod
    def wheatSeptoria(t, y, params):
        S = y[0]
        ER = y[1]
        ES = y[2]
        IR = y[3]
        IS = y[4]
        R = y[5]
        PR = y[6]
        PS = y[7]
        high = y[8]
        low = y[9]

        timeSinceSeasonStart = t % (params.GS87 - params.seasonStartTime) # remainder or the time/(season length), AKA time since season started

        A = S + ER + ES + IR + IS + R # total leaf area

        GS61 = params.GS61 - params.seasonStartTime # finding the time of GS61 relative to the start of the season
        GS87 = params.GS87 - params.seasonStartTime # finding the time of GS87 relative to the start of the season

        if params.seasonality and timeSinceSeasonStart >= GS61: # adding senescence if seasonality is used (it is by default)
            senescence = 0.005 * (timeSinceSeasonStart - GS61) / (GS87 - GS61) + 0.1 * exp(0.02 * (timeSinceSeasonStart - GS87))
        else:
            senescence = 0

        effectOnRInfection = 1 - params.lowRisk.effect(low, FungicideEffectType.Protectant)  # IDFU, uses the effect func from fungicide.py (exponential function)
        effectOnSInfection = 1 - Fungicide.combineEffects( # similar to above, but uses combineEffects, which is (presumably) both
            [params.highRisk, params.lowRisk],
            [high, low],
            FungicideEffectType.Protectant
        )

        infectionByR = params.beta * effectOnRInfection * (IR + PR) # infection term, missing S/A
        infectionByS = params.beta * effectOnSInfection * (IS + PS) # infection term, missing S/A

        if params.densityDependentInfection: # if density dependance is true (it is by default), multiple by S/A (as in the paper)
            infectionByR *= S / A
            infectionByS *= S / A

        deathOfR = params.mu * IR # death of Ir
        deathOfS = params.mu * IS # death of Is

        dydt = [0 for x in range(11)] # vector for the change in each state class

        if params.densityDependentInfection: # if density dependance is true (it is by default), have growth rate dependent on A (as in paper)
            hostGrowth = params.r * (params.k - A)
            dydt[0] = hostGrowth - senescence * S - (infectionByR + infectionByS) # updating change in state S

        if params.latentPeriod: # if latent period is true, find effects of fungicide on maturation from E to I (might be more complex in powdery mildew model, hence the weirdness)
            effectOnRMaturation = 1 - params.lowRisk.effect(low, FungicideEffectType.Eradicant) # IDFU, since 'low' and 'Eradicant', this evaluates to 1
            effectOnSMaturation = 1 - Fungicide.combineEffects( # IDFU, 'low' and 'eradicant' evaluates to 1, but 'high' is same as effectOnSInfection
                [params.highRisk, params.lowRisk],
                [high, low],
                FungicideEffectType.Eradicant
            )

            maturationOfR = params.gamma * ER * effectOnRMaturation # why is the fungicide here?
            maturationOfS = params.gamma * ES * effectOnSMaturation # why is Cl here?

            dydt[1] = infectionByR - senescence * ER - maturationOfR # updating change in state ER
            dydt[2] = infectionByS - senescence * ES - maturationOfS # updating change in state ES
            dydt[3] = maturationOfR - deathOfR # updating change in state IR
            dydt[4] = maturationOfS - deathOfS # updating change in state IS
        else:
            dydt[3] = infectionByR - deathOfR # updating change in state IR (if no latent period)
            dydt[4] = infectionByS - deathOfS # updating change in state IS (if no latent period)

        dydt[5] = deathOfS + deathOfR + senescence * (ER + ES + S) # updating change in state R

        if params.seasonality: # adding decay of initial innoculum if seasonality is used (it is by default)
            decayOfR = params.nu * PR # decay of PR
            decayOfS = params.nu * PS # decay of PS
            dydt[6] = -decayOfR # updating change in state PR
            dydt[7] = -decayOfS # updating change in state PS

        if params.fungicideDecay: # if decay of fungicides is true, adds it to the system
            highRiskDecay = params.highRisk.decayRate * high # 
            lowRiskDecay = params.lowRisk.decayRate * low
            dydt[8] = -highRiskDecay # updating change in fungicide concentration for high-risk
            dydt[9] = -lowRiskDecay # updating change in fungicide concentration for low-risk
        
        dydt[10] = dydt[0] + dydt[1] + dydt[2] + dydt[3] + dydt[4] + dydt[5] # area

        for i, v in enumerate(y):
            if v > 1E15:
                dydt[i] = 0

        return dydt # returns change in each class (same as mine)

    @staticmethod
    def grapePowderyMildew(t, y, params):
        S = y[StateIndices.S]
        ER = y[StateIndices.ER]
        ES = y[StateIndices.ES]
        IR = y[StateIndices.IR]
        IS = y[StateIndices.IS]
        R = y[StateIndices.R]
        O = y[StateIndices.O]
        high = y[StateIndices.high]
        low = y[StateIndices.low]

        tSinceBudBreak = t % (params.seasonEndTime - params.budBreakTime)

        if params.seasonality and tSinceBudBreak >= params.shootToppingTime - params.budBreakTime:
            r = params.r[1]
            k = params.k[1]
            beta = params.beta[1]

            hostGrowthT = tSinceBudBreak - (params.shootToppingTime - params.budBreakTime)
            N0 = 21700 * (1 - params.toppingLoss["S"])  # Assuming all equal
        else:
            r = params.r[0]
            k = params.k[0]
            beta = params.beta[0]
            hostGrowthT = tSinceBudBreak
            N0 = 42.47

        N = S + ER + ES + IR + IS + R + O

        effectOnRInfection = 1 - params.lowRisk.effect(low, FungicideEffectType.Protectant)
        effectOnSInfection = 1 - Fungicide.combineEffects(
            [params.highRisk, params.lowRisk],
            [high, low],
            FungicideEffectType.Protectant
        )
        infectionByR = beta * IR * effectOnRInfection
        infectionByS = beta * IS * effectOnSInfection

        if params.densityDependentInfection:
            infectionByR *= S / N
            infectionByS *= S / N

        deathOfR = params.mu * IR
        deathOfS = params.mu * IS

        dydt = [0 for x in range(11)]

        if params.densityDependentInfection:
            hostGrowth = (r * k * (k / N0 - 1) * math.exp(-r * hostGrowthT)) / math.pow(1 + (k / N0 - 1) * math.exp(-r * hostGrowthT), 2)

            dydt[StateIndices.S] = hostGrowth - (infectionByR + infectionByS)

            if params.seasonality:
                ontogenicResistance = params.m * S
                dydt[StateIndices.S] -= ontogenicResistance
                dydt[StateIndices.O] = ontogenicResistance

        dydt[StateIndices.IR] = -deathOfR
        dydt[StateIndices.IS] = -deathOfS

        if params.latentPeriod:
            effectOnRMaturation = 1 - params.lowRisk.effect(low, FungicideEffectType.Eradicant)
            effectOnSMaturation = 1 - Fungicide.combineEffects(
                [params.highRisk, params.lowRisk],
                [high, low],
                FungicideEffectType.Eradicant
            )
            maturationOfR = params.gamma * ER * effectOnRMaturation
            maturationOfS = params.gamma * ES * effectOnSMaturation

            dydt[StateIndices.ER] = infectionByR - maturationOfR
            dydt[StateIndices.ES] = infectionByS - maturationOfS
            dydt[StateIndices.IR] += maturationOfR
            dydt[StateIndices.IS] += maturationOfS
        else:
            dydt[StateIndices.IR] += infectionByR
            dydt[StateIndices.IS] += infectionByS

        dydt[StateIndices.R] = deathOfS + deathOfR

        if params.fungicideDecay:
            highRiskDecay = params.highRisk.decayRate * high
            lowRiskDecay = params.lowRisk.decayRate * low
            dydt[StateIndices.high] = -highRiskDecay
            dydt[StateIndices.low] = -lowRiskDecay

        return dydt

    def _resistanceFrequency(self, state):
        if state["IS"] == 0 and state["IR"] == 0:
            return self.params.psi
        return state["IR"] / (state["IR"] + state["IS"])

    def _postProcess(self):
        # Calculate the selection ratio first
        if self.params.model == Model.WheatSeptoria:
            endOfSeasonTime = self.params.GS87 - self.params.seasonStartTime
        else:
            endOfSeasonTime = self.params.seasonEndTime - self.params.budBreakTime

        try:
            initialFrequency = self.params.psi
            frequencyAtEnd = self._resistanceFrequency(self.output.loc[endOfSeasonTime - Simulation.outputRate, :])

            if initialFrequency != 0:
                self.selectionRatio = frequencyAtEnd / initialFrequency
        except KeyError:
            self.selectionRatio = None

        # Yield til critical loss
        self.yieldTilCriticalLoss = 0
        self.seasonsTilCriticalLoss = 0

        for thisYield in self.yields:
            self.yieldTilCriticalLoss += thisYield
            if thisYield < (1 - self.params.criticalYieldLoss) * self.diseaseFreeYield:
                break
            self.seasonsTilCriticalLoss += 1
