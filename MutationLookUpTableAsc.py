#!/usr/bin/env python3

import math
import numpy as np

#*#*#*#*#*#*#*#*#*#*#*#*#*#*# NOTES #*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#                                                               #
# 1. 	DeltaGs from file are assumed to be +ve for             #
#       DeltaG_inactive and -ve for DeltaG_active.              #
#                                                               #
# 2.    IC50s from "The specificity of asciminib, a potential   #
#       treatment for chronic myeloid leukemia, as a myristate- #
#       pocket binding ABL inhibitor and analysis of its        #
#       interactions with mutant forms of BCR-ABL1 kinase" by   #
#       Manley et al (2020)
#                                                               #
# 3.    Order and numbering of mutations:                       #
#       0.  Wild-type   1.  A337V                               #
#       2.  I502L       3.  P465S                               #
#							                                    #
# 4.    Drug numbering: 0. Ponatinib, 1. Imatinib,              #
#                       2. Nilotinib, 3. Asciminib Single       #
#                       4. Asciminib Split                      #
#							                                    #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#



#Vmax (pmol/min),kcatAct(/min),KM(microM),DeltaG(kcal/mol),Ponatinib-IC50(nanoM),Imatinib-IC50(nanoM),Nilotinib-IC50(nanoM),Asciminib-IC50(nanoM)
OriginalMutantData = np.array([[1.9,  66.0, 17.0,   0.0,    0.37,  90.5,  3.52,  0.61], \
                               [1.9,  66.0, 17.0,   0.0,    3.52,  82.1,  3.67,  453],\
                               [1.9,  66.0, 17.0,   0.0,    0.29,  56.2,  2.68,  30.2], \
                               [1.9,  66.0, 17.0,   0.0,    3.64,  92.2,  3.30,  369]],float)

OriginalDrugLabels = ["Ponatinib", "Imatinib", "Nilotinib","Asciminib Single Dose","Asciminib Split Dose"]
OriginalMutantLabels = ["Wild type", "A337V", "I502L", "P465S"]

#DayCountSteadyStateBegins - Ponatinib,Imatinib,Nilotinib,AsciminibSingle,AsciminibSplit
DrugDynamicTimings = [7,5,7,4,4] 

#F_Bioavailibility, Dose([10**(-5)]mol), VolumeOfDistribution(L), EliminationRate([10**(-2)/hour), AbsorptionRate(/hour), Dosing interval (hours)
DrugDynamicsConstants = np.array([  [   1, 8.449, 1223, 2.888, 1.302, 24],\
                                    [   1, 81.04,  435, 3.851, 0.940, 24],\
                                    [0.25, 75.54,  597, 4.332, 0.300, 12],\
                                    [   1, 17.78,  111, 8.664, 1.960, 24],\
                                    [   1, 8.892,  111, 8.664, 1.960, 12]],float)


#kATrans, koffSubstrate, Ponatinib koffInhibitorInactive, Imatinib koffInhibitorInactive, Nilotinib koffInhibitorActive, Asciminib koffInhibitorActive (all in /min)
Rates = np.array([[60, 33, 0.00488, 0.059, 0.02, 0.01],\
                  [60, 33, 0.00488, 0.059, 0.02, 0.01],\
                  [60, 33, 0.00488, 0.059, 0.02, 0.01],\
                  [60, 33, 0.00488, 0.059, 0.02, 0.01]],float)

# (kcal/mol)
DeltaGAsc = [3.7, -0.7, 2.7, 2.9, 3.9]

# (nM)


RDInhWhileAsc = np.array([[0.0147,0.0244,0.0095,0.0248],\
                       [3.61,2.66,2.30,2.78],\
                       [0.140,0.143,0.127,0.137]],float)

# (x 10^-3 nM)
RDIAscWhileInh = np.array([[6.381,825.272,258.965,658.403],\
                           [6.381,3851.177,324.821,2919.140],\
                           [6.381,4617.263,374.122,4017.456]],float)

