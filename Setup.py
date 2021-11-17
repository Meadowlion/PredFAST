#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri March 26 10:19:44 2021

@author: sb069
"""

from os.path import expanduser
import os
import errno
os.system("sudo apt-get install python 3.8 >/dev/null")
os.system("pip install scikit-learn >/dev/null")
os.system("pip install matplotlib >/dev/null")
os.system("pip install pickle-mixin >/dev/null")
os.system("pip install seaborn >/dev/null")
os.system("pip install isoelectric  >/dev/null")

Download=str(os.environ['_'])
Downloads=Download.replace('/Setup.py','')

hdir=expanduser('~')
try:
	os.system("sudo apt install hmmer")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Databases'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Output Files'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Input Files'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning/Model'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise

try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning/Prosite'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning/Run Outputs'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning/Run Outputs/Plots'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mkdir -p 'Documents/PredFast/Machine Learning/Run Outputs/Predictions'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
os.chdir(Downloads)
os.system("dir")
Dir=Downloads+"/"
try:
	os.system("mv 'Header.txt' '" + hdir + "/Documents/PredFast/Machine Learning/Header.txt'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'Resources Used' '" + hdir + "/Documents/PredFast'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'Model Output.py' '" + hdir + "/Documents/PredFast/Model Output.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'isol.py' '"+ hdir + "/Documents/PredFast/isol.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'GradientBoostedRandomForest.py' '" + hdir + "/Documents/PredFast/Gradient Boosted Random Forest.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'Database Builder For Unlabeled Dataset.py' '" + hdir + "/Documents/PredFast/Database Builder For Unlabeled Dataset.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'Line.py' '" + hdir + "/Documents/PredFast/Line.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'DatabaseBuilder.py' '" + hdir + "/Documents/PredFast/Database Builder.py'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise

try:
	os.system("mv 'prosite.dat' '" + hdir + "/Documents/PredFast/prosite.dat'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
try:
	os.system("mv 'PFAM' '" + hdir + "/Documents/PredFast/Machine Learning'")
except OSError as e:
	if e.errno != errno.EEXIST:
		raise
def main():
	print("SITE PACKAGE LOCATED HERE")
	os.system("python -m site")
	isoelectric=input("Drag and Drop your _init_.py file from site-packages/isoelectric/_init_.py:  ") 
	with open("config.txt", "w+") as config1:
		config1.write(isoelectric)
		config1.close()
	try:
		os.system("mv 'config.txt' '" + hdir + "/Documents/PredFast/config.txt'")
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise

	print("PredFast is now setup! Happy Protein Research!")
main() 

