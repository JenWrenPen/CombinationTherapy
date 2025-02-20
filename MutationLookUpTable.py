#!/usr/bin/env python3

import math
import numpy as np

#*#*#*#*#*#*#*#*#*#*#*#*#*#*# NOTES #*#*#*#*#*#*#*#*#*#*#*#*#*#*#
#                                                               #
# 1. 	DeltaGs from file are assumed to be +ve for             #
#       DeltaG_inactive and -ve for DeltaG_active.              #
#                                                               #
# 2.    IC50s from two sources A & B.                           #
#       WT, G250E, E255K, E255V, T315I are from Source A        #
#       T315M & Y253H-E255V are from Source B (rest of the      #
#       mutations aren't used in the nuclear approch            #
#       simulation                                              #
#                                                               #
#       A:  Three novel patient-derived BCR/ABL mutants show    #
#           different sensitivity to second and third           #
#           generation tyrosine kinase inhibitors               #
#           - Redaelli et al                                    #
#                                                               #
#       B:  Supplementary material to: BCR-ABL1 Compound        #
#           Mutations Combining Key Kinase Domain Positions     #
#           Confer Clinical Resistance to Ponatinib in Ph       #
#           Chromosome-Positive Leukemia                        #
#           - Zabriskie et al                                   #
#                                                               #
# 3.    Order and numbering of mutations:                       #
#       0.  Wild-type   1.  G250E                               #
#       2.  E255K       3.  E255V                               #
#       4.  T315I                                               #
#							                                    #
# 4.    Drug numbering: 0. Ponatinib, 1. Imatinib,              #
#                       2. Nilotinib, 3. Asciminib Single       #
#                       4. Asciminib Split                      #
#							                                    #
#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#



#Vmax (pmol/min),kcatAct(/min),KM(microM),DeltaG(kcal/mol),Ponatinib-IC50(nanoM),Imatinib-IC50(nanoM),Nilotinib-IC50(nanoM),Asciminib-IC50(nanoM)
OriginalMutantData = np.array([[1.9,  66.0, 17.0,   0.0,    2.1,  527.0,  17.7,  3.8], \
                               [5.1, 175.2, 14.3,  -2.2,   12.5, 3613.0,  80.7,   22], \
                               [1.8,  63.6, 15.6,  -0.5,   17.6,   3174, 118.4,  9.5],\
                               [0.2,   6.8, 22.1,  -0.4,   27.2,   8953, 182.3, 13.3], \
                               [0.4,  12.2,  7.2,   0.1,    6.3,   9221, 697.1, 28.7]],float)

OriginalDrugLabels = ["Ponatinib", "Imatinib", "Nilotinib","Asciminib Single Dose","Asciminib Split Dose"]
OriginalMutantLabels = ["Wild type", "G250E", "E255K", "E255V", "T315I"]

#DayCountSteadyStateBegins - Ponatinib,Imatinib,Nilotinib,AsciminibSingle,AsciminibSplit
DrugDynamicTimings = [7,5,7,4,4] 

#F_Bioavailibility, Dose([10**(-5)]mol), VolumeOfDistribution(L), EliminationRate([10**(-2)/hour), AbsorptionRate(/hour), Dosing interval (hours)
DrugDynamicsConstants = np.array([  [   1, 8.449, 1223, 2.888, 1.302, 24],\
                                    [   1, 81.04,  435, 3.851, 0.940, 24],\
                                    [0.25, 75.54,  597, 4.332, 0.300, 12],\
                                    [   1, 17.78,  111, 8.664, 1.960, 24],\
                                    [   1, 8.892,  111, 8.664, 1.960, 12]],float)


#kATrans, koffSubstrate, Ponatinib koffInhibitorInactive, Imatinib koffInhibitorInactive, Nilotinib koffInhibitorActive, Asciminib koffInhibitorActive (all in /min)
Rates = np.array([[60,  33, 0.00488, 0.059, 0.02, 0.01],\
                  [55, 350, 0.00488, 0.059, 0.02, 0.01],\
                  [60, 127, 0.00488, 0.059, 0.02, 0.01],\
                  [60, 122, 0.00488, 0.059, 0.02, 0.01],\
                  [60,29.7, 0.00488, 0.059, 0.02, 0.01]],float)

# (kcal/mol)
DeltaGAsc = [3.7, -0.7, 2.7, 2.9, 3.9]

# (nM)
RDInhWhileAsc = np.array([[0.0371, 0.0374, 0.0592, 0.0955, 0.0776],\
                       [  7.31,   6.79,  18.74,  17.43,  24.56],\
                       [ 0.168,  0.126,  0.312,  0.383,  0.662]],float)

# (x 10^-3 nM)
RDIAscWhileInh = np.array([[5.766, 0.969, 3.079, 2.577, 10.436],\
                           [5.766, 1.109, 1.418, 2.921,  4.806],\
                           [5.766, 1.584, 2.412, 3.153,  3.717]],float)

