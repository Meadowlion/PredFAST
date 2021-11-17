#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri March 26 10:19:44 2021

@author: sb069
"""

#makes the filepaths global variables

from collections import Counter
import sys
import textwrap
import math
import numpy as np
global iupfix
global nonproteinor
global pfamtxt
global prositetxt
global txt
import os
global g 
global Qtxt
g=0



from isol import *

#for iupred2a
PATH="Resources Used/iupred2a"

def DatabaseBuilder():
    Question()



def Question():
    global nonproteinor
    global pfamtxt
    global prositetxt
    global txt
    global csv3 
    global question1
    global Qtxt
    global g
    global finalin
    global PSpath
    global Inpath
    global Fullpath
    global Prositeq
    Fullpath=os.getcwd()
    PSpath="'" + Fullpath + "/Resources Used/PS_Scan/ps_scan_Linux/" + "'"
    Inpath= Fullpath + "/Input Files"
    Outpath=Fullpath + "/Output Files"
    os.environ['Input']= Inpath
    os.environ['PS']= PSpath 
    os.environ['Output']=Outpath
    os.environ['PSS']= "Resources Used/PS_Scan/ps_scan_Linux/"
    
    if g==0:
        
        HMMscanq=input("Do you want to run HMMscan?[Y/n]    ")
        Prositeq=input("Do you want to run prosite scan?[Y/n]    ")
    else:
        pass
        
    nonproteinor=""
    question1=input("Is this the enzyme , non-enzyme or an unlabeled file:" )
    question1=question1.strip()
    if question1=="enzyme":
        nonproteinor="1"
    if question1==float:
        question1=str(question1)
        
        
    elif question1=="non-enzyme" or "unlabeled":
        nonproteinor="0"
        
    else:
        print("Please enter 'enzyme' , 'non-enzyme' or 'unlabeled' ")
        return Question()
    txt=input("Please enter the name of the fasta file you wish to use:" )
    Qtxt="Input Files/" + txt
    csv3=input("Please enter what you wish to name the output file:" )
    csv3="Output Files/" + csv3 + ".csv" 
    pfamtxt="Machine Learning/PFAM/" + txt #+"' " 
    prositetxt="Machine Learning/Prosite/" + txt 
    print("The class number is " + nonproteinor)
  
    
     
    if Prositeq =="n" and HMMscanq == "n":
        count()
    elif Prositeq == "n" and HMMscanq=="Y":
        print("Remember to have the Prosite file in FASTA and in 'LM_ML/Machine Learning/Prosite' with the same name as the input file")
        HMMERSCAN()
    elif Prositeq== "Y"  and HMMscanq=="n":
        print("Remember to have the Pfam file in FASTA and in 'LM_ML/Machine Learning/Pfam' with the same name as the input file")
        PrositeScan()
    elif Prositeq=="Y" and HMMscanq=="Y":
        HMMERSCAN()
    else:
        print("Invalid Selection. Both programs will run.")
        HMMERSCAN()
def HMMERSCAN():
    global txt
    global Qtxt
    global hmmin
    global Fullpath
    global Prositeq
    hmmin="'" +Qtxt+ "'"
    print("The dataset will now be built Hmmer is running")
    hmmstr="hmmscan --domtblout 'Machine Learning/PFAM/" + txt +"'" + " 'Machine Learning/PFAM/Pfam-A.hmm' " + hmmin + " >/dev/null"
    
    os.system(hmmstr)
    if Prositeq=="Y":
        print("Hmmer has completed Prosite Scan is Running")
        PrositeScan()
    elif Prositeq=="n":
        print("HMMScan has finished please wait for the dataset to be fully built")

        count()
    else:
        print("Invalid Selection The program will now restart")
        PrositeScan()
def PrositeScan():
    global txt
    global Qtxt
    global hmmin
    import subprocess
    global Inpath
    #PSin= Inpath+ "/"  +txt
    prositein= "'" + "Machine Learning/Prosite/" + txt + "'" 
    psstr1="perl 'Resources Used/PS_Scan/ps_scan_Linux/ps_scan.pl' " + hmmin + " -o fasta > " + prositein
    psscan='Resources Used/PS_Scan/ps_scan_Linux/ps_scan.pl'
    # THE LINE MUST LOOK EXACTLY LIKE THIS IN ps_scan.pl line 690 my $PFSCAN  = "'/home/sb069/Documents/LM_ML/Resources Used/PS_Scan/ps_scan_Linux/pfscan'";
    exportpath= "'" + Fullpath + "/Resources Used/PS_Scan/ps_scan_Linux/pfscan'"
    exportpath='"'+exportpath+'"'
    exportpath="my $PFSCAN  =" + exportpath + ";"
    #print(exportpath)
    
    # This is how it appears in proscan we are going to manually replace path in PERL program: my $PFSCAN  = 'pfscan';
    
    with open(psscan, 'r') as prositefile:
        prositedata=prositefile.read()
        prositedata=prositedata.replace("my $PFSCAN  = 'pfscan';" , exportpath)
    prositefile.close()
    with open(psscan, 'w') as prositeoutfile:
        prositeoutfile.write(prositedata)
    prositeoutfile.close()
    #print("Export PATH", exportpath)
    #PSSTRING=exportpath +' ; ' + psstr1
    #print("THIS IS PSSTRING", PSSTRING)
    #PSPATH2=PSpath+":"
    #os.putenv("PATH", PSPATH2 +
    #os.getenv("PATH")
    #my_env=os.environ.copy()
    #my_env["PATH"] = PSPATH2 + my_env["PATH"]
    #print(my_env["PATH"])
    #subprocess.call(["perl", psscan, PSin],env=my_env)
    subprocess.call(psstr1, shell= True)
    #os.system(psstr1)
    #os.system(exportpath +' ; ' + psstr1)
    #termout.close()          
    #env={"PATH":PSpath,"2PATH":"/usr/bin/env perl"})
    #subprocess.call(psstr1,shell=True,env=myenv)
    #subprocess.call(["perl", "Resources Used/PS_Scan/ps_scan_Linux/ps_scan.pl", Qtxt, "-o", "fasta", ">" , prositetxt])#,
                    #env={"PATH":PSpath,"2PATH":"/usr/bin/env perl"}) #This seems to correctly change path but perl is not found? Perhaps need to set two PATHs?
    #print(proscanstr)
    #subprocess.call(psstr1,shell=True) 
    print("Prosite Scan has finished please wait for the dataset to be fully built")
    
    count()
#ignore this till you get to into() all iupred2a stuff
def avg(lst):
    return sum(lst) / len(lst)

def aa_freq(_seq):
    _freq = {}
    for _aa in _seq:
        if _aa in _freq:
            _freq[_aa] += 1
        else:
            _freq[_aa] = 1
    for _aa, _ins in _freq.items():
        _freq[_aa] = _ins / len(_seq)
    return _freq


def read_matrix(matrix_file):
    _mtx = {}
    with open(r"Resources Used/iupred2a/data/iupred2_long_energy_matrix", "r") as _fhm:
        for _line in _fhm:
            if _line.split()[0] in _mtx:
                _mtx[_line.split()[0]][_line.split()[1]] = float(_line.split()[2])
            else:
                _mtx[_line.split()[0]] = {}
                _mtx[_line.split()[0]][_line.split()[1]] = float(_line.split()[2])
    return _mtx


def read_histo(histo_file):
    hist = []
    h_min = float("inf")
    h_max = -float("inf")
    with open(histo_file, "r") as fnh:
        for _line in fnh:
            if _line.startswith("#"):
                continue
            if float(_line.split()[1]) < h_min:
                h_min = float(_line.split()[1])
            if float(_line.split()[1]) > h_max:
                h_max = float(_line.split()[1])
            hist.append(float(_line.split()[-1]))
    h_step = (h_max - h_min) / (len(hist))
    return hist, h_min, h_max, h_step


def smooth(energy_list, window):
    weighted_energy_score = [0] * len(energy_list)
    for idx in range(len(energy_list)):
        weighted_energy_score[idx] = avg(energy_list[max(0, idx - window):min(len(energy_list), idx + window + 1)])
    return weighted_energy_score


def iupred(seq, mode):
    #long 
    #can try short changing by changin to uc=25 and replace /long_ with /short_
    lc = 1
    uc = 100
    wc = 10
    mtx = read_matrix("{}/data/iupred2_long_energy_matrix".format(PATH))
    histo, histo_min, histo_max, histo_step = read_histo("{}/data/long_histogram".format(PATH))
                                                                                                                                                                                           
    unweighted_energy_score = [0] * len(seq)
    weighted_energy_score = [0] * len(seq)
    iupred_score = [0] * len(seq)

    for idx in range(len(seq)):
        freq_dct = aa_freq(seq[max(0, idx - uc):max(0, idx - lc)] + seq[idx + lc + 1:idx + uc + 1])
        for aa, freq in freq_dct.items():
            try:
                unweighted_energy_score[idx] += mtx[seq[idx]][aa] * freq
            except KeyError:
                unweighted_energy_score[idx] += 0

    if mode == 'short':
        for idx in range(len(seq)):
            for idx2 in range(idx - wc, idx + wc + 1):
                if idx2 < 0 or idx2 >= len(seq):
                    weighted_energy_score[idx] += -1.26
                else:
                    weighted_energy_score[idx] += unweighted_energy_score[idx2]
            weighted_energy_score[idx] /= len(range(idx - wc, idx + wc + 1))
    else:
        weighted_energy_score = smooth(unweighted_energy_score, wc)

    glob_text = ""
    if mode == 'glob':
        gr = []
        in_gr = False
        beg, end = 0, 0
        for idx, val in enumerate(weighted_energy_score):
            if in_gr and val <= 0.3:
                gr.append({0: beg, 1: end})
                in_gr = False
            elif in_gr:
                end += 1
            if val > 0.3 and not in_gr:
                beg = idx
                end = idx
                in_gr = True
        if in_gr:
            gr.append({0: beg, 1: end})
        mgr = []
        k = 0
        kk = k + 1
        if gr:
            beg = gr[0][0]
            end = gr[0][1]
        nr = len(gr)
        while k < nr:
            if kk < nr and gr[kk][0] - end < 45:
                beg = gr[k][0]
                end = gr[kk][1]
                kk += 1
            elif end - beg + 1 < 35:
                k += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
            else:
                mgr.append({0: beg, 1: end})
                k = kk
                kk += 1
                if k < nr:
                    beg = gr[k][0]
                    end = gr[k][1]
        seq = seq.lower()
        nr = 0
        res = ""
        for i in mgr:
            res += seq[nr:i[0]] + seq[i[0]:i[1] + 1].upper()
            nr = i[1] + 1
        res += seq[nr:]
        res = " ".join([res[i:i + 10] for i in range(0, len(res), 10)])
        glob_text += "Number of globular domains: {}\n".format(len(mgr))
        for n, i in enumerate(mgr):
            glob_text += "          globular domain   {}.\t{}-{}\n".format(n + 1, i[0] + 1, i[1] + 1)
        glob_text += "\n".join(textwrap.wrap(res, 70))

    for idx, val in enumerate(weighted_energy_score):
        if val <= histo_min + 2 * histo_step:
            iupred_score[idx] = 1
        elif val >= histo_max - 2 * histo_step:
            iupred_score[idx] = 0

        else:
            iupred_score[idx] = histo[int((weighted_energy_score[idx] - histo_min) * (1 / histo_step))]
    return iupred_score, glob_text


def anchor2(seq, iupred_scores):
    local_window_size = 41
    iupred_window_size = 30
    local_smoothing_window = 5
    par_a = 0.0013
    par_b = 0.26
    par_c = 0.43
    iupred_limit = par_c - (par_a / par_b)
    mtx = read_matrix('{}/data/anchor2_energy_matrix'.format(PATH))
    interface_comp = {}
    with open('{}/data/anchor2_interface_comp'.format(PATH)) as _fn:
        for line in _fn:
            interface_comp[line.split()[1]] = float(line.split()[2])
    local_energy_score = [0] * len(seq)
    interface_energy_score = [0] * len(seq)
    energy_gain = [0] * len(seq)
    for idx in range(len(seq)):
        freq_dct = aa_freq(seq[max(0, idx - local_window_size):max(0, idx - 1)] + seq[idx + 2:idx + local_window_size + 1])
        for aa, freq in freq_dct.items():
            try:
                local_energy_score[idx] += mtx[seq[idx]][aa] * freqprint(txt3)
            except KeyError:
                local_energy_score[idx] += 0
        for aa, freq in interface_comp.items():
            try:
                interface_energy_score[idx] += mtx[seq[idx]][aa] * freq
            except KeyError:
                interface_energy_score[idx] += 0
        energy_gain[idx] = local_energy_score[idx] - interface_energy_score[idx]
    iupred_scores = smooth(iupred_scores, iupred_window_size)
    energy_gain = smooth(smooth(energy_gain, local_smoothing_window), local_smoothing_window)
    anchor_score = [0] * len(seq)
    for idx in range(len(seq)):
        sign = 1
        if energy_gain[idx] < par_b and iupred_scores[idx] < par_c:
            sign = -1
        corr = 0
        if iupred_scores[idx] > iupred_limit and energy_gain[idx] < 0:
            corr = (par_a / (iupred_scores[idx] - par_c)) + par_b
        anchor_score[idx] = sign * (energy_gain[idx] + corr - par_b) * (iupred_scores[idx] - par_c)
        anchor_score[idx] = 1 / (1 + math.e ** (-22.97968 * (anchor_score[idx] - 0.0116)))
    return anchor_score

# all variables that go back to count need to be in here as global not at the top the strip stuff is due to how iupred2a outputs things
def into(iupfix):
    global mean
    global idp
    sequence = iupfix
    iupred2_result = iupred(sequence, sys.argv[-1])
    iupred2_result = iupred(sequence, sys.argv)
    iupred2_result=str(iupred2_result)
    iupred2_result=iupred2_result.strip("]")
    str_result=str(iupred2_result)
    str_result=str_result.strip("()")
    str_result=str_result.strip("]")
    str_result=str_result.strip("[")
    str_result=str_result.strip("''")
    str_result=str_result.split(",")
    str_result=str_result[:-2] 
    str_result=np.array(str_result).astype(np.float)
    l=len(str_result)
    H=sum(str_result)
    mean=H/l
    if mean>0.5:
        idp=1
    if mean<0.5:
        idp=0
    
    
    
def count():
    global nonproteinor
    global pfamtxt
    global prositetxt
    global txt
    global csv3 
    global question1 
    global g 
    global Qtxt

    PS=""
    
    with open(Qtxt, "r") as txt3:
    
    

        proteinout = []
#actual creation of data from sequence
        for line in txt3 :
            if ">" in line:
               
                PS="0"  #defines PS and PF as String
                PF="0"
                name=str(line[1:7]) #length of query name
                with open(prositetxt,"r") as r:
                    for line in r:
                        if name in line: 
                            PS=str(line.split(": PS",1)[1]) #allows us to get a string for Prosite
                            PS=PS[0:5]
                with open(pfamtxt,"r") as p:
                    x=0
                    for line in p:
                        if name in line and x<1:   #the x<1 just returns the top result setting it higher returns a different result IE 2 would return the 2nd result for that query name
                            PF=str(line.split("    PF",1)[1])
                            PF=PF[0:5]  #returns just the part of the list that is the PF family 
                            x+=1         #counter increases by one per run allows us to return the top sequence 
                        
                        
                
                    
                        
        
            elif ">" not in line:
                sequence=line #important for pI
                iupfix=sequence #important for iupred2a
                counted = Counter(line)
                lin= ((len(line))-1)
                lins= str(lin)
                #needed for iupred2a
                into(iupfix)
            
                input_pKa_set='IPC_protein' #this uses the same pka set they use
                alph= str((counted["A"]*1.420 + counted["R"]*0.980 + counted["N"]*0.670 + counted["D"]*1.010  + counted["C"]*0.700 + counted["Q"]*1.110  + counted["E"]*1.510  + counted["G"]*0.570 + counted["H"]*1.000  + counted["I"]*1.080 + counted["L"]*1.210 + counted["K"]*1.160 + counted["M"]*1.450 + counted["F"]*1.130 + counted["P"]*0.570 + counted["S"]*0.770 + counted["T"]*0.830  + counted["W"]*1.080 + counted["Y"]*0.690  + counted["V"]*1.060 ) / lin) #alpha helix %
                beta= str((counted["A"]*0.709 + counted["R"]*0.920 + counted["N"]*0.604 + counted["D"]*0.541 + counted["C"]*1.191  + counted["Q"]*0.840  + counted["E"]*0.567 + counted["G"]*0.657 + counted["H"]*0.863  + counted["I"]*1.799 + counted["L"]*1.261 + counted["K"]*0.721 + counted["M"]*1.210 + counted["F"]*1.393 + counted["P"]*0.354 + counted["S"]*0.928  + counted["T"]*1.221  + counted["W"]*1.306 + counted["Y"]*1.266  + counted["V"]*1.965 ) / lin) #beta strand %
                bulk= str((counted["A"]*11.500 + counted["R"]*14.280 + counted["N"]*12.820 + counted["D"]*11.680  + counted["C"]*13.460 + counted["Q"]*14.450  + counted["E"]*13.570  + counted["G"]*3.400 + counted["H"]*13.690  + counted["I"]*21.400 + counted["L"]*21.400 + counted["K"]*15.710 + counted["M"]*16.250 + counted["F"]*19.800 + counted["P"]*17.430 + counted["S"]*9.470 + counted["T"]*15.770 + counted["W"]*21.670 + counted["Y"]*18.030 + counted["V"]*21.570) / lin)#bulkiness
                ref= str((counted["A"]*4.340 + counted["R"]*26.660 + counted["N"]*13.280 + counted["D"]*12.000 + counted["C"]*35.770 + counted["Q"]*17.560 + counted["E"]*17.260 + counted["G"]*0.000 + counted["H"]*21.810 + counted["I"]*19.060 + counted["L"]*18.780 + counted["K"]*21.290 + counted["M"]*21.640 + counted["F"]*29.400 + counted["P"]*10.930 + counted["S"]*6.350 + counted["T"]*11.010 + counted["W"]*42.530 + counted["Y"]*31.530 + counted["V"]*13.920) / lin)#refractivity
                pI= str(predict_isoelectric_point(sequence, input_pKa_set)) #pI
                coil= str((counted["A"]* 0.824 + counted["R"]*0.893 + counted["N"]*1.167 + counted["D"]*1.197  + counted["C"]*0.953  + counted["Q"]*0.947 + counted["E"]*0.761 + counted["G"]* 1.251  + counted["H"]*1.068 + counted["I"]*0.886 + counted["L"]*0.810  + counted["K"]*0.897 + counted["M"]*0.810 + counted["F"]*0.797 + counted["P"]* 1.540 + counted["S"]*1.130 + counted["T"]*1.148 + counted["W"]*0.941 + counted["Y"]*1.109 + counted["V"]*0.772) / lin) # random coil
                buried= str((counted["A"]*11.200 + counted["R"]*0.500 + counted["N"]*2.900 + counted["D"]*2.900  + counted["C"]*4.100  + counted["Q"]*1.600 + counted["E"]*1.800 + counted["G"]* 11.800  + counted["H"]*2.000 + counted["I"]*8.600 + counted["L"]*11.700 + counted["K"]*0.500 + counted["M"]*1.900 + counted["F"]*5.100 + counted["P"]*2.700 + counted["S"]*8.000 + counted["T"]*2.200  + counted["W"]*0.941 + counted["Y"]*2.600 + counted["V"]*12.900) / lin) # buried residues
                hratio= str((counted["A"]*0.000 + counted["R"]*0.650 + counted["N"]*1.330 + counted["D"]*1.380  + counted["C"]*2.750   + counted["Q"]*0.890  + counted["E"]*0.920 + counted["G"]* 0.740  + counted["H"]* 0.580 + counted["I"]* 0.000+ counted["L"]*0.000 + counted["K"]*0.330 + counted["M"]*0.000 + counted["F"]*0.000 + counted["P"]*0.390 + counted["S"]*1.420 + counted["T"]*0.710   + counted["W"]*0.130 + counted["Y"]*0.200 + counted["V"]*0.000) / lin) #Hetero Ratio
                cdensity= str((counted["A"]*1.401 + counted["R"]*1.1 + counted["N"]*1.66 + counted["D"]*1.54  + counted["C"]*0.0  + counted["Q"]*0 + counted["E"]*1.460 + counted["G"]*1.607  + counted["H"]*0.0 + counted["I"]*0.0+ counted["L"]*1.191+ counted["K"]*0.0 + counted["M"]*1.340 + counted["F"]*0 + counted["P"]*0 + counted["S"]*1.537 + counted["T"]*0.0   + counted["W"]*0.0 + counted["Y"]*1.456 + counted["V"]*1.230) / lin) # crystal density prediction
                #for solubility some were listed as only "very soluble" at g/mL these values were given into the database as 100
                dalt= str(counted["A"]*89 + counted["R"]*174 + counted["N"]*132 + counted["D"]*133 + counted["C"]*121 + counted["Q"]*146 + counted["E"]*147 + counted["G"]*75 + counted["H"]*155 + counted["I"]*131 + counted["L"]*131 + counted["K"]*146 + counted["M"]*149 + counted["F"]*165 + counted["P"]*115 + counted["S"]*105 + counted["T"]*119 + counted["W"]*204 + counted["Y"]*181 + counted["V"]*117)#counts daltons
                polar=  str((counted["A"]*8.100 + counted["R"]*10.500 + counted["N"]*11.600  + counted["D"]*13.000 + counted["C"]*5.500 + counted["Q"]*10.500  + counted["E"]*12.300 + counted["G"]*9.000 + counted["H"]*10.400 + counted["I"]*5.200 + counted["L"]*4.900 + counted["K"]*11.300 + counted["M"]*5.700 + counted["F"]*5.200 + counted["P"]*8.000 + counted["S"]*9.200 + counted["T"]*8.600 + counted["W"]*5.400 + counted["Y"]*6.200  + counted["V"]*5.900) / lin ) #gets average polarity for the protein
                hydro= str((counted["A"]*-0.400 + counted["R"]*-0.590 + counted["N"]*-0.920 + counted["D"]*-1.310 + counted["C"]*0.170 + counted["Q"]*-0.910 + counted["E"]*-1.220 + counted["G"]*-0.670 + counted["H"]*-0.640 + counted["I"]*1.250 + counted["L"]*1.220 + counted["K"]*-0.670 + counted["M"]*1.020 + counted["F"]*1.920 + counted["P"]*-0.490 + counted["S"]*-0.550 + counted["T"]*-0.280 + counted["W"]*0.500 + counted["Y"]* 1.670 + counted["V"]*0.910) / lin)  # gets average hydrophobicity for protein    
                flex=  str((counted["A"]*0.360 + counted["R"]*0.530 + counted["N"]*0.460 + counted["D"]*0.510 + counted["C"]*0.350 + counted["Q"]*0.490 + counted["E"]* 0.500 + counted["G"]*0.540 + counted["H"]*0.320 + counted["I"]*0.460 + counted["L"]*0.370 + counted["K"]*0.470 + counted["M"]*0.300 + counted["F"]*0.310 + counted["P"]*0.510 + counted["S"]*0.510 + counted["T"]*0.440 + counted["W"]*0.310 + counted["Y"]*0.420 + counted["V"]*0.390) / lin) #gets average flexibility for the protein              
                proteinout.append(name + "," + str(counted["A"]/ lin)+ "," + str(counted["R"]/ lin)+ "," + str(counted["N"]/ lin)+ "," + str(counted["D"]/ lin)+ "," + str(counted["C"]/ lin)+ "," + str(counted["Q"]/ lin)+ "," + str(counted["E"]/ lin)+ "," + str(counted["G"]/ lin)+ "," + str(counted["H"]/ lin)+ "," + str(counted["I"]/ lin)+ "," + str(counted["L"]/ lin)+ "," + str(counted["K"]/ lin)+ "," + str(counted["M"]/ lin)+ "," + str(counted["F"]/ lin)+ "," + str(counted["P"]/ lin)+ "," + str(counted["S"]/ lin)+ "," + str(counted["T"]/ lin)+ "," + str(counted["W"]/ lin)+ "," + str(counted["Y"]/ lin)+ "," + str(counted["V"]/ lin)+ "," + lins + "," + dalt + "," + hydro + "," + polar + "," + flex + "," + pI + "," + ref + "," + bulk +"," + alph+ "," + beta + "," + coil+ "," + buried + "," + hratio + "," + cdensity + "," + str(idp) + "," + str(mean) + "," + str(PS) + "," +  str(PF) + "," + nonproteinor) #this is what actually makes the dataset   #this is what actually makes the dataset   
                
    
    
    proteinfinal='\n'.join(proteinout)
#removing final and last to remove lines
#closing file 
    txt3.close()
    global g
# Make a new file
    global csv3
    
    
    csv1= open(csv3, "w")
    
    csv1.write(proteinfinal)
    
    csv1.close()
    r.close()
    p.close()
    
    print("This has made the " + question1 + " part of the file have you run both enzymes and non-enzymes")
    print("Output Files located in the folder 'Output Files")
    if g<1 and question1!="unlabeled":
        g+=1
        DatabaseBuilder()
    elif g>1 or g==1:
        Finalizer()
    else:
        Finalizer()


def Finalizer():
    global csv3
    global txt
    final=""
    header="Machine Learning/Header.txt"
    Ftxt=csv3
    Ftxt2='no file'
    #Ftxt3=input("Please enter what you wish to name the database: ")
    Ftxt3=r"Databases/" + txt + ".csv"
    if Ftxt2=="no file":
         with open(header) as header1:
             header2=header1.read()
    
         with open(Ftxt) as read1:
             data1=read1.read()
         final+= header2
         final+= "\n"
         final+= data1
         with open (Ftxt3, 'w+') as output: 
            output.write(final) 
         print("The dataset is now ready to be used in an algorithm the file will be found in 'Databases' ")
    else:   
        Ftxt2="Output Files/" + Ftxt2
    
        with open(header) as header1:
            header2=header1.read()
    
        with open(Ftxt) as read1:
            data1=read1.read()
        
        with open(Ftxt2) as read2:
            data2=read2.read()
    
        final+= header2
        final+= "\n"
        final+= data1
        final+= "\n"
        final+=data2
        with open (Ftxt3, 'w+') as output: 
            output.write(final) 
        print("The full dataset is now ready to be used in an algorithm the file will be found in 'Databases'")


DatabaseBuilder()




          
        
    
    
    
