import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def n_plot(xlab, ylab, xs=14, ys=14):
    plt.minorticks_on()
    plt.tick_params(axis='both', which='major', labelsize=ys - 2, direction='in', length=6, width=2)
    plt.tick_params(axis='both', which='minor', labelsize=ys - 2, direction='in', length=4, width=2)
    plt.tick_params(axis='both', which='both', top=True, right=True)
    if xlab is not None:
        plt.xlabel(xlab, fontsize=xs)
    if ylab is not None:
        plt.ylabel(ylab, fontsize=ys)
    plt.tight_layout()
    return None

Beta = 1/(310.15*(1.380649*10**-23))

# Small Optimised Ground State Energy Results

CC_N_E = -1.833311543*1.609*10**-19
CG_N_E = -1.591993476*1.609*10**-19
GC_N_E = -1.556782136*1.609*10**-19
GG_N_E = -1.865596193*1.609*10**-19

TT_N_E = -1.57321388*1.609*10**-19
TA_N_E = -1.599604388*1.609*10**-19
AT_N_E = -1.602403587*1.609*10**-19
AA_N_E = -1.416247587*1.609*10**-19

Alpha_CC_N = np.exp(-Beta*CC_N_E)
Alpha_CG_N = np.exp(-Beta*CG_N_E)
Alpha_GC_N = np.exp(-Beta*GC_N_E)
Alpha_GG_N = np.exp(-Beta*GG_N_E)

Alpha_TT_N = np.exp(-Beta*TT_N_E)
Alpha_TA_N = np.exp(-Beta*TA_N_E)
Alpha_AT_N = np.exp(-Beta*AT_N_E)
Alpha_AA_N = np.exp(-Beta*AA_N_E)

part_func_N = Alpha_CC_N + Alpha_CG_N + Alpha_GC_N + Alpha_GG_N + Alpha_TT_N + Alpha_TA_N + Alpha_AT_N + Alpha_AA_N

Prob_CC_N = Alpha_CC_N/part_func_N
Prob_CG_N = Alpha_CG_N/part_func_N
Prob_GC_N = Alpha_GC_N/part_func_N
Prob_GG_N = Alpha_GG_N/part_func_N

Prob_TT_N = Alpha_TT_N/part_func_N
Prob_TA_N = Alpha_TA_N/part_func_N
Prob_AT_N = Alpha_AT_N/part_func_N
Prob_AA_N = Alpha_AA_N/part_func_N

# Small_BB Optimised Ground State Energy Results

CC_BB_E = -3.423114235*1.609*10**-19
CG_BB_E = -3.222166225*1.609*10**-19
GC_BB_E = -3.38863218*1.609*10**-19
GG_BB_E = -3.596662282*1.609*10**-19

TT_BB_E = -3.602782924*1.609*10**-19
TA_BB_E = -3.283056113*1.609*10**-19
AT_BB_E = -3.498333269*1.609*10**-19
AA_BB_E = -4.04622536*1.609*10**-19

Alpha_CC_BB = np.exp(-Beta*CC_BB_E)
Alpha_CG_BB = np.exp(-Beta*CG_BB_E)
Alpha_GC_BB = np.exp(-Beta*GC_BB_E)
Alpha_GG_BB = np.exp(-Beta*GG_BB_E)

Alpha_TT_BB = np.exp(-Beta*TT_BB_E)
Alpha_TA_BB = np.exp(-Beta*TA_BB_E)
Alpha_AT_BB = np.exp(-Beta*AT_BB_E)
Alpha_AA_BB = np.exp(-Beta*AA_BB_E)

part_func_BB = Alpha_CC_BB + Alpha_CG_BB + Alpha_GC_BB + Alpha_GG_BB + Alpha_TT_BB + Alpha_TA_BB + Alpha_AT_BB + Alpha_AA_BB

Prob_CC_BB = Alpha_CC_BB/part_func_BB
Prob_CG_BB = Alpha_CG_BB/part_func_BB
Prob_GC_BB = Alpha_GC_BB/part_func_BB
Prob_GG_BB = Alpha_GG_BB/part_func_BB

Prob_TT_BB = Alpha_TT_BB/part_func_BB
Prob_TA_BB = Alpha_TA_BB/part_func_BB
Prob_AT_BB = Alpha_AT_BB/part_func_BB
Prob_AA_BB = Alpha_AA_BB/part_func_BB

# Large_BB Optimised Ground State Energy Results

CC_4BB_E = -3.454068971*1.609*10**-19
CG_4BB_E = -4.372963242*1.609*10**-19
GC_4BB_E = -3.941490439*1.609*10**-19
GG_4BB_E = -3.481601376*1.609*10**-19

TT_4BB_E = -4.560740036*1.609*10**-19
TA_4BB_E = -4.415570863*1.609*10**-19
AT_4BB_E = -4.284213745*1.609*10**-19
AA_4BB_E = -4.51688244*1.609*10**-19

Alpha_CC_4BB = np.exp(-Beta*CC_4BB_E)
Alpha_CG_4BB = np.exp(-Beta*CG_4BB_E)
Alpha_GC_4BB = np.exp(-Beta*GC_4BB_E)
Alpha_GG_4BB = np.exp(-Beta*GG_4BB_E)

Alpha_TT_4BB = np.exp(-Beta*TT_4BB_E)
Alpha_TA_4BB = np.exp(-Beta*TA_4BB_E)
Alpha_AT_4BB = np.exp(-Beta*AT_4BB_E)
Alpha_AA_4BB = np.exp(-Beta*AA_4BB_E)

part_func_4BB = Alpha_CC_4BB + Alpha_CG_4BB + Alpha_GC_4BB + Alpha_GG_4BB + Alpha_TT_4BB + Alpha_TA_4BB + Alpha_AT_4BB + Alpha_AA_4BB

Prob_CC_4BB = Alpha_CC_4BB/part_func_4BB
Prob_CG_4BB = Alpha_CG_4BB/part_func_4BB
Prob_GC_4BB = Alpha_GC_4BB/part_func_4BB
Prob_GG_4BB = Alpha_GG_4BB/part_func_4BB

Prob_TT_4BB = Alpha_TT_4BB/part_func_4BB
Prob_TA_4BB = Alpha_TA_4BB/part_func_4BB
Prob_AT_4BB = Alpha_AT_4BB/part_func_4BB
Prob_AA_4BB = Alpha_AA_4BB/part_func_4BB

prob_2BP_NoBB_Opti = [Prob_CC_N, Prob_GG_N, Prob_GC_N, Prob_CG_N, Prob_TT_N, Prob_AA_N, Prob_TA_N, Prob_AT_N]
prob_2BP_BB_Opti = [Prob_CC_BB, Prob_CG_BB, Prob_GC_BB, Prob_GG_BB, Prob_TT_BB, Prob_AA_BB, Prob_TA_BB, Prob_AT_BB]
prob_4BP_BB_Opti = [Prob_CC_4BB, Prob_CG_4BB, Prob_GC_4BB, Prob_GG_4BB, Prob_TT_4BB, Prob_AA_4BB, Prob_TA_4BB, Prob_AT_4BB]

prob_2BP_NoBB_Opti_Fix = []
prob_2BP_BB_Opti_Fix = []
prob_4BP_BB_Opti_Fix = []

#Estimated Experimental Prediction
prob_Experi = [61.54, 23.08, 6.154, 9.231]

for j in range(4):
    prob_2BP_NoBB_Opti_Fix.append(prob_2BP_NoBB_Opti[j*2]*100 + prob_2BP_NoBB_Opti[j*2+1]*100)
    prob_2BP_BB_Opti_Fix.append(prob_2BP_BB_Opti[j*2]*100 + prob_2BP_BB_Opti[j*2+1]*100)
    prob_4BP_BB_Opti_Fix.append(prob_4BP_BB_Opti[j*2]*100 + prob_4BP_BB_Opti[j*2+1]*100)


index = ['CC&GG', 'CG&GC', 'TT&AA', 'TA&AT']
df = pd.DataFrame({
                    '2BP_NoBB_Opti': prob_2BP_NoBB_Opti_Fix,
                    '2BP_BB_Opti': prob_2BP_BB_Opti_Fix,
                    '4BP_BB_Opti': prob_4BP_BB_Opti_Fix,
                   'Experiment': prob_Experi}, index=index)
ax = df.plot.bar(rot=0, color={"2BP_NoBB_Opti": "red","Experiment": "black", "2BP_BB_Opti": "blue", "4BP_BB_Opti": "Purple"}, width=0.9)
n_plot(None, None, xs=18, ys=18)
plt.ylabel('Probability (%)', fontsize=18)
plt.savefig("Output.pdf", bbox_inches='tight')
exit()