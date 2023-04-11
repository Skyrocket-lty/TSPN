# -*- codeing = utf-8 -*-
# @Time : 2022/10/18 21:20
# @Author : 刘体耀
# @File : TNNRSMMA.py
# @Software: PyCharm

import numpy as np
import pandas as pd
import math
import numpy.matlib
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import copy

SM_miRNA_M = np.loadtxt(r"SM-miRNA association matrix.txt", dtype=int)
SM_miRNA_M_2 = np.loadtxt(r"SM-miRNA association matrix 2.txt", dtype=int)
SM = np.loadtxt(r"SM similarity matrix.txt", dtype=int)
miRNA = np.loadtxt(r"miRNA similarity matrix.txt", dtype=int)
SM_2 = np.loadtxt(r"SM similarity matrix 2.txt", dtype=int)
miRNA_2 = np.loadtxt(r"miRNA similarity matrix 2.txt", dtype=int)

#基于数据集1计算SM高斯轮廓核相似性
def Gaussian_SM():
    row=831
    sum=0
    SM1=np.matlib.zeros((row,row))
    for i in range(0,row):
        a=np.linalg.norm(SM_miRNA_M[i,])*np.linalg.norm(SM_miRNA_M[i,])
        sum=sum+a
    ps=row/sum
    for i in range(0,row):
        for j in range(0,row):
            SM1[i,j]=math.exp(-ps*np.linalg.norm(SM_miRNA_M[i,]-SM_miRNA_M[j,])*np.linalg.norm(SM_miRNA_M[i,]-SM_miRNA_M[j,]))
            if(SM[i,j]==0):
                SM[i,j]=SM1[i,j]
    return SM
#基于数据集2计算mirna高斯轮廓核相似性
def Gaussian_MM():
    column=541
    sum=0
    miRNA1=np.matlib.zeros((column,column))
    for i in range(0,column):
        a=np.linalg.norm(SM_miRNA_M[:,i])*np.linalg.norm(SM_miRNA_M[:,i])
        sum=sum+a
    ps=column/sum
    for i in range(0,column):
        for j in range(0,column):
            miRNA1[i,j]=math.exp(-ps*np.linalg.norm(SM_miRNA_M[:,i]-SM_miRNA_M[:,j])*np.linalg.norm(SM_miRNA_M[:,i]-SM_miRNA_M[:,j]))
            if(miRNA[i,j]==0):
                miRNA[i,j]=miRNA1[i,j]
    return miRNA




#基于数据集2计算SM高斯轮廓核相似性
def Gaussian_SM_2():
    row=39
    sum=0
    SM1=np.matlib.zeros((row,row))
    for i in range(0,row):
        a=np.linalg.norm(SM_miRNA_M_2[i,])*np.linalg.norm(SM_miRNA_M_2[i,])
        sum=sum+a
    ps=row/sum
    for i in range(0,row):
        for j in range(0,row):
            SM1[i,j]=math.exp(-ps*np.linalg.norm(SM_miRNA_M_2[i,]-SM_miRNA_M_2[j,])*np.linalg.norm(SM_miRNA_M_2[i,]-SM_miRNA_M_2[j,]))
            if(SM_2[i,j]==0):
                SM_2[i,j]=SM1[i,j]
    return SM_2
#基于数据集2计算mirna高斯轮廓核相似性
def Gaussian_MM_2():
    column=286
    sum=0
    miRNA1=np.matlib.zeros((column,column))
    for i in range(0,column):
        a=np.linalg.norm(SM_miRNA_M_2[:,i])*np.linalg.norm(SM_miRNA_M_2[:,i])
        sum=sum+a
    ps=column/sum
    for i in range(0,column):
        for j in range(0,column):
            miRNA1[i,j]=math.exp(-ps*np.linalg.norm(SM_miRNA_M_2[:,i]-SM_miRNA_M_2[:,j])*np.linalg.norm(SM_miRNA_M_2[:,i]-SM_miRNA_M_2[:,j]))
            if(miRNA_2[i,j]==0):
                miRNA_2[i,j]=miRNA1[i,j]
    return miRNA_2

if __name__ == "__main__":
    SM_GIPK = Gaussian_SM(SM_miRNA_M)
    miRNA_GIPK = Gaussian_MM(SM_miRNA_M)
    SM_GIPK_2 = Gaussian_SM_2(SM_miRNA_M_2)
    miRNA_GIPK_2 = Gaussian_MM_2(SM_miRNA_M_2)
