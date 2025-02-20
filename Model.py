#!/usr/bin/env python3

import numpy as np
import math
import random
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
import sympy as sy

from MutationLookUpTable import * # imports values for Vmax, kCatAct, KM, IC50, DeltaG, Names, etc

#*#*#*#*#*#*#*#*#*#*#*#*#*#*# NOTES #*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#                                                               #
# 1.    Asciminib regime 3 is once daily, and 4 is twice daily  #
# 2.    The timestep should be about 1/100th of a second        #
#                                                               #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#


################ Constant parameters ################

Temp = 310 # K
kB = 3.29982916 * 10**(-27) #in kcal/K
Beta = 1/(kB*Temp) # per kcal

SubConc = 10 * 10 ** (-6) # M
BaseInhConc = 1 * (10 ** (-9)) #changes nM to M

Etotal = 1*10**(-6) # M

WildTypeDeltaG = -1 # kcal/mol

MutationNumbersWanted = [0,1,2,3,4]
DrugsWanted = [0,1,2]
AsciminibCourseWanted = [3,4]
CombinationInhibitor = [0, 25, 50, 75, 100]
CombinationAsciminib = [0, 25, 50, 75, 100]

DrugName = ["Ponatinib","Imatinib","Nilotinib","Asciminib Single", "Asciminib Split"]
MutationNames = ["wild-type", "G250E","E255K", "E255V", "T315I"]

DayTotal = 10 # in days
TimeStep = 0.001 # in seconds
StepTotal = int(DayTotal * 24 * 3600 / TimeStep)
TimeTotal = DayTotal*24*3600

RateCatalysis = 0
RateTransActive = 0
RateOffSubstrate = 0
RateTransInactive = 0
RateOnSubstrate = 0
RateOnInhibitor = 0
RateOnAsciminib = 0

PrevActUnb = 0
PrevActBou = 0
PrevInaUnb = 0
PrevInaInh = 0
PrevActAsc = 0
PrevInaAsc = 0

PrevInhConc = 0
PrevAscConc = 0

SwitchDay = 0
Gamma = 0
Epsilon = 0
Alpha = 0
DosingInterval = 1

ThisDrug = 0

#plotter=0

################ Functions ################

#*# Convert per minute to per time step #*#

def convertRateToPerSec(RATE):
    return RATE/60

#*# Constants for concentration cycle equations #*#

def calculateGamma(FbIOAVAILABILITY,dOSE,aBSORPTIONrATE,eLIMINATIONrATE,vOLUMEoFdISTRIBUTION):
    return (FbIOAVAILABILITY * dOSE * aBSORPTIONrATE) / (vOLUMEoFdISTRIBUTION * (aBSORPTIONrATE - eLIMINATIONrATE))

def calculateEpsilon(eLIMINATIONrATE,dOSEiNTERVAL):
    return math.exp(eLIMINATIONrATE*dOSEiNTERVAL)

def calculateAlpha(aBSORPTIONrATE,dOSEiNTERVAL):
    return math.exp(aBSORPTIONrATE*dOSEiNTERVAL)

#*# Concentration cycle equations #*#

def findInhibitorConcentration(tIME,sWITCHdAY,gAMMA,ePSILON,aLPHA,aBSORPTIONrATE,eLIMINATIONrATE,dOSEiNTERVAL):

    DayCounter = math.floor(tIME/(3600*24))
    DoseCounter = math.floor(tIME/(3600*dOSEiNTERVAL))
    TimeSinceDose = (tIME / 3600) - (DoseCounter * dOSEiNTERVAL) #in hours
    TimeInHours = tIME / 3600

    if (DayCounter < sWITCHdAY):
        
        SumEpsilon = 0
        SumAlpha = 0

        for Dose in range(DoseCounter+1):

            SumEpsilon = SumEpsilon + (ePSILON ** Dose)  # no units
            SumAlpha = SumAlpha + (aLPHA ** Dose)        # no units

        return gAMMA * ((SumEpsilon * math.exp(-eLIMINATIONrATE * TimeInHours)) - (SumAlpha * math.exp(-aBSORPTIONrATE * TimeInHours)))

    else:
        return gAMMA * (((ePSILON * math.exp(-eLIMINATIONrATE * TimeSinceDose))/(ePSILON - 1)) - ((aLPHA * math.exp(-aBSORPTIONrATE * TimeSinceDose))/(aLPHA - 1)))

#*# Concentration of enzyme state equations #*#

#Decay#

def functionDecay(rate,time):
    return 1-rate*time

#Growth#

def functionGrowth(rate,time):
    return 1-functionDecay(rate,time)

#Concentration#

def findInhibitorBound():
    A = PrevInaInh * functionDecay((RateOffInhibitor + RateOnAsciminibInh * PrevAscConc),TimeStep)
    B = PrevInaUnb * functionGrowth((RateOnInhibitor * PrevInhConc),TimeStep)  
    C = PrevInhAsc * functionGrowth((RateOffAsciminib),TimeStep)
    return A + B + C

def findInactiveUnbound():
    A = PrevInaUnb * functionDecay((RateOnInhibitor * PrevInhConc + RateTransInactive + RateOnAsciminibIna * PrevAscConc),TimeStep)
    B = PrevActUnb * functionGrowth((RateTransActive),TimeStep)
    C = PrevInaInh * functionGrowth((RateOffInhibitor),TimeStep)
    D = PrevInaAsc * functionGrowth((RateOffAsciminib),TimeStep)
    return A + B + C + D

def findActiveUnbound():
    A = PrevActUnb * functionDecay((RateTransActive + RateOnSubstrate * SubConc + RateOnAsciminibAct * PrevAscConc),TimeStep)
    B = PrevInaUnb * functionGrowth((RateTransInactive),TimeStep)
    C = PrevActBou * functionGrowth((RateCatalysis + RateOffSubstrate),TimeStep)
    D = PrevActAsc * functionGrowth((RateOffAsciminib),TimeStep)
    return A + B + C + D

def findSubstrateBound():
    A = PrevActBou * functionDecay((RateCatalysis + RateOffSubstrate),TimeStep)
    B = PrevActUnb * functionGrowth((RateOnSubstrate * SubConc),TimeStep)
    return A + B

def findActiveAsciminib():
    A = PrevActAsc * functionDecay((RateOffAsciminib + RateTransAscActive),TimeStep)
    B = PrevActUnb * functionGrowth((RateOnAsciminibAct * PrevAscConc),TimeStep)
    C = PrevInaAsc * functionGrowth((RateTransAscInactive),TimeStep)
    return A + B + C

def findInactiveAsciminib():
    A = PrevInaAsc * functionDecay((RateOffAsciminib + RateTransAscInactive + RateOnInhibitorAsc * PrevInhConc),TimeStep)
    B = PrevInaUnb * functionGrowth((RateOnAsciminibIna * PrevAscConc),TimeStep)
    C = PrevActAsc * functionGrowth((RateTransAscActive),TimeStep)
    D = PrevInhAsc * functionGrowth((RateOffInhibitor),TimeStep)
    return A + B + C + D

def findInhibitedAsciminib():
    A = PrevInhAsc * functionDecay((RateOffAsciminib + RateOffInhibitor),TimeStep)
    B = PrevInaInh * functionGrowth((RateOnAsciminibInh * PrevAscConc),TimeStep)
    C = PrevInaAsc * functionGrowth((RateOnInhibitorAsc * PrevInhConc),TimeStep)
    return A + B + C


################ Calculations ################

plt.clf()

## Asciminib concentration values (ThisDrug=3)##

for ThisAscPercent in CombinationAsciminib:
    for ThisInhPercent in CombinationInhibitor:

        for ThisCourse in AsciminibCourseWanted:

            AscSwitchDay = DrugDynamicTimings[ThisCourse]
            AscFBioavailability = DrugDynamicsConstants[ThisCourse][0]                                  # no units
            AscDose = (ThisAscPercent/100) * DrugDynamicsConstants[ThisCourse][1] * 10 ** (-5)    # mol
            AscVolumeOfDistribution = DrugDynamicsConstants[ThisCourse][2]                              # litre
            AscEliminationRate = DrugDynamicsConstants[ThisCourse][3] * 10 ** (-2)                      # per hour
            AscAbsorptionRate = DrugDynamicsConstants[ThisCourse][4]                                    # per hour
            AscDosingInterval = DrugDynamicsConstants[ThisCourse][5]                                    # hours

            AscGamma = calculateGamma(AscFBioavailability,AscDose,AscAbsorptionRate,AscEliminationRate,AscVolumeOfDistribution) # M or mol/L
            AscEpsilon = calculateEpsilon(AscEliminationRate,AscDosingInterval)                                                 # no units
            AscAlpha = calculateAlpha(AscAbsorptionRate,AscDosingInterval)                                                      # no units


            for ThisDrug in DrugsWanted:

                SwitchDay = DrugDynamicTimings[ThisDrug]
                FBioavailability = DrugDynamicsConstants[ThisDrug][0]                           # no units
                Dose = (ThisInhPercent/100) * DrugDynamicsConstants[ThisDrug][1] * 10 ** (-5)   # mol
                VolumeOfDistribution = DrugDynamicsConstants[ThisDrug][2]                       # litre
                EliminationRate = DrugDynamicsConstants[ThisDrug][3] * 10 ** (-2)               # per hour
                AbsorptionRate = DrugDynamicsConstants[ThisDrug][4]                             # per hour
                DosingInterval = DrugDynamicsConstants[ThisDrug][5]                             # hours

                Gamma = calculateGamma(FBioavailability,Dose,AbsorptionRate,EliminationRate,VolumeOfDistribution)   # M or mol/L
                Epsilon = calculateEpsilon(EliminationRate,DosingInterval)                                          # no units
                Alpha = calculateAlpha(AbsorptionRate,DosingInterval)                                               # no units

                for ThisMutation in MutationNumbersWanted:

                    kCat = OriginalMutantData[ThisMutation][1]                                          # per minute
                    KM = OriginalMutantData[ThisMutation][2] * 10 ** (-6)                               # M
                    DeltaG = (WildTypeDeltaG + OriginalMutantData[ThisMutation][3]) /(6.022 * 10 ** 23) # kcal
                    IC50 = OriginalMutantData[ThisMutation][ThisDrug + 4] * 10 ** (-9)                  # M
                    AscIC50 = OriginalMutantData[ThisMutation][7] * 10 ** (-9)                          #M
                    
                    RateCatalysis = convertRateToPerSec(kCat)                                           # per second
                    RateTransActive = convertRateToPerSec(Rates[ThisMutation][0])                       # per second
                    RateTransAscActive = RateTransActive * math.exp(-Beta * 1.5 /(6.022 * 10 ** 23))    # per second
                    RateOffSubstrate = convertRateToPerSec(Rates[ThisMutation][1])                      # per second

                    RateOffInhibitor = convertRateToPerSec(Rates[ThisMutation][ThisDrug+2])             # per second
                    RateOffAsciminib = convertRateToPerSec(Rates[ThisMutation][ThisCourse])             # per second

                    RateOffInhibitorAsc = RateOffInhibitor
                    RateOffAsciminibInh = RateOffAsciminib
                  
                    ExpNegDelG = math.exp(-Beta * DeltaG)   # no unit
                    ExpPosDelG = math.exp(Beta * DeltaG)    # no unit
                    OnePlusSOverKM = 1 + (SubConc / KM)     # no unit

                    RD = IC50 / (1 + ExpNegDelG * OnePlusSOverKM)               # M
                    n = ExpNegDelG / 0.1                                        # no unit
                    AscRDI = 1.1 * AscIC50 / (1 + ExpNegDelG * OnePlusSOverKM)  # M
                    AscRDA = n * AscRDI                                         # M

                    RateTransInactive = RateTransActive * ExpNegDelG
                    RateTransAscInactive = RateTransInactive * math.exp(-Beta * DeltaGAsc[ThisMutation] /(6.022 * 10 ** 23))
                    RateOnSubstrate = (RateOffSubstrate + RateCatalysis) / KM
                    RateOnInhibitor = RateOffInhibitor / RD
                    RateOnAsciminibIna = RateOffAsciminib / AscRDI
                    RateOnAsciminibAct = RateOffAsciminib / AscRDA

                    RateOnAsciminibInh = RateOnAsciminibIna
                    RateOnInhibitorAsc = RateOnInhibitor

                    WeightInaUnb = 1                            # no unit
                    WeightActUnb = ExpNegDelG                   # no unit
                    WeightActBou = ExpNegDelG * SubConc / KM    # no unit

                    NonInhWeightTotal = WeightInaUnb + WeightActUnb + WeightActBou
                         
                    InitActUnb = Etotal * WeightActUnb / NonInhWeightTotal
                    InitActBou = Etotal * WeightActBou / NonInhWeightTotal
                    InitInaUnb = Etotal * WeightInaUnb / NonInhWeightTotal
                    InitInaInh = 0
                    InitActAsc = 0
                    InitInaAsc = 0
                    InitInhAsc = 0

                    PrevActUnb = InitActUnb
                    PrevActBou = InitActBou
                    PrevInaUnb = InitInaUnb
                    PrevInaInh = InitInaInh
                    PrevActAsc = InitActAsc
                    PrevInaAsc = InitInaAsc
                    PrevInhAsc = InitInhAsc

                    PrevInhConc = 0
                    PrecAscConc = 0

                    ActUnb = 0
                    ActBou = 0
                    InaUnb = 0
                    InaInh = 0
                    ActAsc = 0
                    InaAsc = 0
                    InhAsc = 0
                    
                    ProductInTimeStep = 0
                    TotalProduct = 0

                    Time=0
                    RecordedTime = []
                    InhConc = []
                    AscConc = []
                    ProductRate= []
                    SubBound = []
                    Inhibited = []
                    Asciminibed = []
                    InhibitedAsciminibed = []
                    Inactive = []
                    Active = []
                    Total = []

                    CurrentTotal=0

                    SecondsInMinute = 0
                    
                    ProductRateFile = 'ProductRateMutant{:02d}Drug{:02d}P{:03d}AndAsciminib{:02d}P{:03d}.dat'.format(ThisMutation,ThisDrug,ThisInhPercent,ThisCourse,ThisAscPercent)
                    OUTPUT = open(ProductRateFile,"w+")

                    while Time<TimeTotal :

                        
                        if Time==0:
                            ActBou = PrevActBou
                            ActUnb = PrevActUnb
                            InaUnb = PrevInaUnb
                            InaInh = PrevInaInh
                            ActAsc = PrevActAsc
                            InaAsc = PrevInaAsc
                            InhAsc = PrevInhAsc
                        else:
                            ActBou = findSubstrateBound()
                            ActUnb = findActiveUnbound()
                            InaUnb = findInactiveUnbound()
                            InaInh = findInhibitorBound()
                            ActAsc = findActiveAsciminib()
                            InaAsc = findInactiveAsciminib()
                            InhAsc = findInhibitedAsciminib()

                        PrevInhConc=findInhibitorConcentration(Time,SwitchDay,Gamma,Epsilon,Alpha,AbsorptionRate,EliminationRate,DosingInterval)
                        PrevAscConc=findInhibitorConcentration(Time,AscSwitchDay,AscGamma,AscEpsilon,AscAlpha,AscAbsorptionRate,AscEliminationRate,AscDosingInterval)

                        if SecondsInMinute>=60:
                            RecordedTime.append(Time/(24*3600))
                            SubBound.append(ActBou)
                            ProductRate.append(ActBou * kCat)
                            OUTPUT.write('%e  %e'%((Time/(24*3600)),(ActBou * kCat)))
                            Inhibited.append(InaInh)
                            Asciminibed.append(ActAsc+InaAsc)
                            InhibitedAsciminibed.append(InhAsc)
                            Inactive.append(InaUnb)
                            Active.append(ActUnb)
                            InhConc.append(PrevInhConc*10**9)
                            AscConc.append(PrevAscConc*10**9)
                            CurrentTotal=ActBou+ActUnb+InaUnb+InaInh+ActAsc+InaAsc+InhAsc
                            Total.append(CurrentTotal)
                            OUTPUT.write('  %e  %e  %e  %e  %e  %e  %e  %e  %e\n'%((PrevInhConc*10**9),(PrevAscConc*10**9),ActBou,ActUnb,InaUnb,InaInh,ActAsc,InaAsc,InhAsc))
                            SecondsInMinute=0
                        else:
                            SecondsInMinute = SecondsInMinute + TimeStep


                        if (CurrentTotal>=(Etotal*1.01)):
                            Time=TimeTotal
                        
                        PrevActUnb = ActUnb
                        PrevActBou = ActBou
                        PrevInaUnb = InaUnb
                        PrevInaInh = InaInh
                        PrevActAsc = ActAsc
                        PrevInaAsc = InaAsc
                        PrevInhAsc = InhAsc

                        Time=Time+TimeStep

                    if ThisMutation==0:
                        outputfilename = 'InhibitorConcentrationsDrug{:02d}P{:03d}AndAsciminib{:02d}P{:03d}.png'.format(ThisDrug,ThisInhPercent,ThisCourse,ThisAscPercent)
                        plt.plot(RecordedTime,InhConc, label='{drug} {percent}%%'.format(drug=DrugName[ThisDrug],percent=ThisInhPercent))
                        plt.plot(RecordedTime,AscConc, label='{drug} {percent}%%'.format(drug=DrugName[ThisCourse],percent=ThisAscPercent))
                        plt.xlabel('Time (days)', fontsize = 12)
                        plt.ylabel('Inhibitor Concentration (nM)', fontsize = 12)
                        plt.legend('',frameon=False)
                        plt.savefig(outputfilename)
                        plt.clf()

                    outputfilename = 'ProductRateMutant{:02d}Drug{:02d}P{:03d}AndAsciminib{:02d}P{:03d}.png'.format(ThisMutation,ThisDrug,ThisInhPercent,ThisCourse,ThisAscPercent)
                    plt.plot(RecordedTime,ProductRate)
                    plt.xlabel('Time (days)', fontsize = 12)
                    plt.ylabel('Product Rate', fontsize = 12)
                    plt.legend('',frameon=False)
                    plt.savefig(outputfilename)

                    plt.clf()

                    outputfilename = 'EnzymeStatesMutant{:02d}Drug{:02d}P{:03d}AndAsciminib{:02d}P{:03d}.png'.format(ThisMutation,ThisDrug,ThisInhPercent,ThisCourse,ThisAscPercent)
                    plt.plot(RecordedTime,SubBound, label='Substrate bound')
                    plt.plot(RecordedTime,Active, label='Active')
                    plt.plot(RecordedTime,Inactive, label='Inactive')
                    plt.plot(RecordedTime,Inhibited, label='Inhibited')
                    plt.plot(RecordedTime,Asciminibed, label='Asciminib Bound')
                    plt.plot(RecordedTime,InhibitedAsciminibed, label='Asciminib and inhibitor Bound')
                    plt.plot(RecordedTime,Total, label='Total')
                    plt.xlabel('Time (days)', fontsize = 12)
                    plt.ylabel('Enzyme Concentration', fontsize = 12)
                    plt.legend()
                    plt.savefig(outputfilename)

                    plt.clf()


