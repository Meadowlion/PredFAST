#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 18:26:58 2021

@author: sb069
"""

import os 
#opening document using open 
In = input("What is the name of the file you wish to make buildable (ensure its in the 'Input Files' folder)?: ")

Out= input("What would you like to name the output file Warning: naming it the same as above will overwrite it: ")
In= "Input Files/" + In
In1=open(In, "r")

proteinlist=[]

for line in In1:
    if ">" in line:

        h="\n"+line
        proteinlist.append(h)
    else:
        proteinlist.append(line.strip("\n"))
proteinfinal="".join(proteinlist)
#removing final and last to remove lines
#print(proteinfinal.strip())
#closing file 


# Make a new file
Out1 = open("Input Files/" + Out + ".txt", "w+")
#write to file 
Out1.write(proteinfinal)
#close file
Out1.close()
In1.close()
Database=input("Do you want to launch database builder? [Y/n]:  ")
if Database=="Y"
	os.system('python3 '"'DatabaseBuilder.py'" )
if Database=="n"
	pass
