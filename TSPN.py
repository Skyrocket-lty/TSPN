# -*- codeing = utf-8 -*-
# @Time : 2022/10/18 21:20
# @Author : 刘体耀
# @File : TSPN.py
# @Software: PyCharm

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc
import copy


# SVT: singular value thresholding operator for matrix Y by thretholding parameter b and Weight sequence h
def SVT(Y, b, h):
    U,S,V= np.linalg.svd(Y)
    for index in range(0,S.size):
        if S[index] >= b*h[index]:
            S[index] = S[index] - b*h[index]
        else:
            S[index] = 0
    s = np.diag(S)
    row , col = Y.shape[0] , Y.shape[1]
    if row < col:
        s_n = np.column_stack((s, np.zeros((row, col - row))))
    else:
        s_n = np.row_stack((s, np.zeros((row-col, col))))
    Y_n = np.dot(U, np.dot(s_n, V))
    return Y_n

# TSPN: Iterative solution process of ADMM algorithm.
def TSPN(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r):
    X = t
    N = X
    Y = X
    iter0 = 1
    stop1 = 1
    stop2 = 1

    # the processing of computing W
    T=np.dot(B.transpose(),A)
    U,S,V= np.linalg.svd(T)
    u,s,v= np.linalg.svd(X)
    W=np.zeros_like(s)
    for i in range(0,S.size):
        W[i]=p*(1-S[i])*pow(s[i],p-1)
        if(W[i]<0):
            W[i]=0

    while stop1 > tol1 or stop2 > tol2:
        # the processing of computing N
        tran = (1/beta) * (Y+alpha*(t*omega))+X
        N = tran - (alpha/(alpha+beta))*omega*tran
        N[N < 0] = 0
        N[N > 1] = 1

        # the processing of computing X
        X_1 = SVT(N-(1/beta)*Y, 1/beta,W)

        # the processing of computing Y
        Y = Y + gamma*(X_1-N)

        stop1_0 = stop1
        if np.linalg.norm(X) != 0:
            stop1 = np.linalg.norm(X_1-X) / np.linalg.norm(X)
        else:
            stop1 = np.linalg.norm(X_1-X)
        stop2 = np.abs(stop1-stop1_0)/(max(1, np.abs(stop1_0)))
        X = X_1

        if iter0 >= maxiter:
            iter0 = maxiter
            print('reach maximum iteration,did not converge!')
            break
        iter0 = iter0 + 1
    T_recover = W
    return T_recover, iter0


def run_MC(t):
    # TSPN parameter
    maxiter = 300
    alpha = 1
    beta = 10
    gamma = 1
    p=0.8
    tol1 = 2 * 1e-3
    tol2 = 1 * 1e-5
    omega = np.zeros(t.shape)
    omega[t.nonzero()] = 1
    #插入外层循环，或者叫第一步循环
    for i in range(0,3):
        U, S, V = np.linalg.svd(t)
        r =1
        A = U[:r, :]
        B = V[:r, :]
        t, k = TSPN(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r)
    t_new = t
    M_1 = t_new[0:SM.shape[0], SM.shape[0]:T.shape[1]]
    return M_1

#M为异构图关联矩阵，左上为SM相似性矩阵，右上为SM—Mirna关联矩阵，右下为Mirna相似性矩阵，左下为SM—Mirna相似性矩阵的转置,维数为1372*1372
#SM为小分子相似性矩阵（831*831
M = np.loadtxt(r'adjaceny matrix T.txt', dtype=float)
SM = np.loadtxt(r"SM similarity matrix.txt",dtype=float)

if_name_ == "_main_":
    Scores_M = run_MC(M) #Scores_M is the prediction score matrix






