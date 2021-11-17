# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 22:29:37 2020

@author: meado
"""

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn import metrics
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import confusion_matrix 
import matplotlib.pyplot as plt
import seaborn as sns
import pickle 
from sklearn import model_selection
import os

csv0=input("What is the name of the test file?: ")
csv2=input("What model do you want to use?: ")
runout=csv0 
csv1="Databases/" + runout
proteindata = pd.read_csv(csv1, header=0)
print("The machine will now begin!")
X_test=proteindata[["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","Amino Acids","Daltons","Hydrophobicity","Polarity","Flex","pI","Refractivity","Bulk","Alpha Helix","Beta Sheet","Coil","Buried","Hetero Ratio","Crystal Density","Prosite","PFAM"]]

    #"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","Amino Acids","Daltons","Hydrophobicity","Polarity","Flex","pI","Refractivity","Bulk","Alpha Helix","Beta Sheet","Coil","Buried","Hetero Ratio","Crystal Density","Prosite","PFAM"
y_test=proteindata["Class"]
loaded_model = pickle.load(open("Machine Learning/Model/" + csv2, "rb"))
result = loaded_model.predict(X_test)
print("This will be found in Run Outputs/Predictions")
predout="Machine Learning/Run Outputs/Predictions/" + csv0 + " vs " + csv2 + ".csv"
out=open(predout,"x")
out.close()
ypred = pd.DataFrame(data=result)
df1 = X_test
df2=df1.insert(38,"Predicted Class",ypred,True)

df1.to_csv(predout)
