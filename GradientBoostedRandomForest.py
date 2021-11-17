# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 20:51:13 2019

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

n_jobs=10
s=0
z=0
csv1=input("What is the name of the Training Databse?: ")
labels=list(input("What are the classes of the database? (no commas) "))
print(labels)
runout=csv1 
csv1="Databases/" + csv1
print("The machine will now begin!")
while z<100:
    x=0
    acc=0.0
    mac=0.0
    proteindata = pd.read_csv(csv1, header=0)
    X=proteindata[["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","Amino Acids","Daltons","Hydrophobicity","Polarity","Flex","pI","Refractivity","Bulk","Alpha Helix","Beta Sheet","Coil","Buried","Hetero Ratio","Crystal Density","Prosite","PFAM"]]

    #"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","Amino Acids","Daltons","Hydrophobicity","Polarity","Flex","pI","Refractivity","Bulk","Alpha Helix","Beta Sheet","Coil","Buried","Hetero Ratio","Crystal Density","IDP","Average Disorder","Prosite","PFAM"
    y=proteindata["Class"]
            
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)
    while x<3:
        
                                                    
        #X1 and y1 define our test data's labels
        #gradient boosted similar to random forest
        clf = GradientBoostingClassifier(max_depth=50, n_estimators=500, learning_rate=0.1,min_samples_leaf=20)
       
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        errors = [mean_squared_error(y_test, y_pred) for y_pred in clf.staged_predict(X_test)]
        best_n_estimators = np.argmin(errors)
        #this is what "boosts" the random forest needs to be a int                       =best_n_estimators
        best_n_estimators=int(best_n_estimators)
        best_clf = GradientBoostingClassifier(max_depth=50, n_estimators=best_n_estimators, learning_rate=0.1,min_samples_leaf=20)
        best_clf.fit(X_train, y_train)
        y_pred = best_clf.predict(X_test)
        print("MAE:", mean_absolute_error(y_test, y_pred))
        print("MSE:",mean_squared_error(y_test,y_pred))
        mac+=mean_absolute_error(y_test, y_pred)
        accuracy=str(metrics.accuracy_score(y_test, y_pred))
        print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
        acc+=metrics.accuracy_score(y_test, y_pred)
        filename = str("Machine Learning/Model/" + runout + accuracy + ".dat" ) 
        out=open(filename,"x")
        out.close()
        pickle.dump(best_clf, open(filename, 'ab'))
         #view feature importance as a data set or a plot if graphed
        feature_importances = pd.DataFrame(clf.feature_importances_,index = X_train.columns,columns=['importance']).sort_values('importance',ascending=False)
        print(feature_importances)
        x+=1

    
    
    avg=(acc/x)
    avmac=(mac/x)
    strac=str(avmac)
    print("The Average is", avg)
    print("The Average Mean Error is", avmac)
    #makes a writable string of the above
    stravg=str(avg)
    #allows append instead of rewriting
    txt1 = open("Machine Learning/Run Outputs/" + runout ,"a")
    txt1.write("\n"+stravg+","+strac+"\n")
    txt1.close
   
    #This creates a confusion matrix and shows us what it is guessing labels are the classes defined in the csv
    cm=confusion_matrix(y_test, y_pred,labels=labels)
    #this plots out our confusion matrix as a heatmap arraw
    pl=pd.DataFrame(cm)

    sns.heatmap(pl, annot=True,fmt= "d")
    z+=1
    s+=1
    png="Machine Learning/Run Outputs/Plots/" + runout + str(s) + ".png"
    plt.savefig(png, dpi=1000)
    plt.show()
    
