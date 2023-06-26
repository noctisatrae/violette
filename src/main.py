import numpy as np 
import matplotlib.pyplot as plt
import numpy.linalg as alg
from scipy.linalg import svd, sqrtm
from sklearn import linear_model

from rpy2.robjects.packages import importr
import rpy2.robjects as robjects
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri

rpy2.robjects.numpy2ri.activate()

package = robjects.r.source("./src/solvebeta.r", encoding="utf-8")

solvebeta = robjects.globalenv['solvebeta']
convcheck_bis = robjects.globalenv['convcheck']

fichier = open('pitprops_cor.txt','r')
X=[]
for ligne in fichier :
    X.append(ligne.split())
X=np.array(X[1:]).T[1:]
C=np.zeros((13,13))
for i in range(13):
    for j in range(13):
        C[i][j]=float(X[i][j])



def convcheck(beta1, beta2):
    a = np.apply_along_axis(np.max, 0, np.abs(beta1 + beta2))
    b = np.apply_along_axis(np.max, 0, np.abs(beta1 - beta2))
    d = len(a)
    x = np.ones(d)
    for i in range(d):
        x[i] = min(a[i], b[i])
    return np.max(x)

def SPCA_corr(C,k, lam, Mu, eps):
    n,p = 180,13
    X = sqrtm(C)
    U, S, V = alg.svd(X, full_matrices = False)
    A = V.T[:,:k]
    y = np.array(np.dot (X,A))
    XtX=X.T.dot(X)
    B = A.copy()
    vartot =  np.sum([S[i]**2 for i in range(min(n,p))])
    X_nr, X_nc =  np.shape(X)
    X_converted_array = ro.r.matrix(X, nrow=X_nr, ncol=X_nc)
    for i in range(k):
        _,y_nr = np.shape([y[:,i]])

        y_converted_array =  ro.r.matrix(np.array([y[:,i]]).T, nrow=y_nr, ncol=1)
        paras = ro.r.matrix(np.array([lam,Mu[i]]),nrow=1, ncol=2)
        B[:,i]=solvebeta(X_converted_array,y_converted_array,paras,650,'penalty')
    print('B',B)
    temp = B.copy()
    normtemp = np.sqrt(np.sum(temp ** 2, axis=0)) #tableau des normes de chaque ligne
    normtemp[normtemp == 0] = 1 #remplace les normes nulles par 1 pour Ã©viter la division par 0

    temp = (temp / normtemp)
    diff = 1
    compt = 0
    if n> p :
        while compt<200 and diff>eps:
            compt+=1
            U,_, V = alg.svd(np.dot(XtX,B), full_matrices = False)
            A = np.dot(U,V)#il semble que dans l un ce soit v' et l autre v
            y=np.dot(X,A)
            for i in range(k):
                _,y_nr = np.shape([y[:,i]])

                y_converted_array =  ro.r.matrix(np.array([y[:,i]]).T, nrow=y_nr, ncol=1)
                paras = ro.r.matrix(np.array([lam,Mu[i]]),nrow=1, ncol=2)
                B[:,i]=solvebeta(X_converted_array,y_converted_array,paras,650,'penalty')
            normbeta = np.sqrt(np.sum(B ** 2, axis=0))
            normbeta[normbeta == 0] = 1 #eviter la division par 0
            beta2 = (B / normbeta)
            diff = convcheck(temp,beta2)
            temp = beta2.copy()
    else :
        while compt<200 :
            compt+=1
            liste = [A[:,i].T*C for i in range(k)]
            B = [max(abs(liste[i]-Mu[i]*0.5),0)*np.sign(liste[i]) for i in range(k)]
            Ap = U*V.T
            if sum([np.linalg.norm(B[i],1) for i in range(k)]) <eps:
                print("arret")
                break
    normbeta = np.sqrt(np.sum(B ** 2, axis=0))
    normbeta[normbeta == 0] = 1 #eviter la division par 0
    B = (B / normbeta)
    _,R = np.linalg.qr(np.dot(X,B))
    varex=np.diag(R**2)/vartot
    return B,varex,compt
    
print(SPCA_corr(C,6,0,[0.06, 0.16, 0.1, 0.5, 0.5, 0.5],1e-3))