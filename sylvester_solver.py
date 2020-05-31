import numpy as np
import os
import math
from scipy.sparse import csr_matrix
import scipy
import scipy.linalg
"""
An implementation of the Generalized Sylvester Equation solver.

            AXB + CXD = E
using the Bartels-Stewart algorithm.

References
[1]  Chu, K.W.E., 1987. The solution of the matrix equations AXB - CXD= E AND (YA - DZ, YC -BZ)=(E, F). Linear Algebra and its applications, 93, pp.93-105.
[2]  Hern ́andez,  V.  and  Gass ́o,  M.,  1989.  Explicit  solution  of  the  matrix  equation  AXB  -  CXD=  E.Linear Algebra and its applications, 121, pp.333-344.
[3]  Bartels, R.H. and Stewart, G.W., 1972. Solution of the matrix equation AX + XB= C [F4]. Com-munications of the ACM, 15(9), pp.820-826.
[4]  Gardiner,  J.D.,  Laub,  A.J.,  Amato,  J.J. and Moler,  C.B.,  1992. Solution of the Sylvester matrixequationAXBT-CXDT=  E.  ACM  Transactions  on  Mathematical  Software  (TOMS),  18(2),pp.223-231.

This work is inspired from the MATLAB implementation of Lyapunov solver freeLYAP by Nick Hale, Alex Townsend, and Heather Wilber available at https://github.com/ajt60gaibb/freeLYAP.
"""
def sylvester_solver(A,B,C,D,E):

    tol = 1e-08
    
    #assuming uniqueness if RHS is zero then the solution is zero
    if(np.linalg.norm(E)<tol):
        X = np.zeros_like(E)
        return

    #checking for standard Lyapunov and Sylvester equation
    if not (len(B) and len(C)and A.all()==D.all()):
        X = scipy.linalg.solve_lyapunov(A,E)
        return X

    A = np.asarray(csr_matrix(A).todense())
    B = np.asarray(csr_matrix(B).todense())
    C = np.asarray(csr_matrix(C).todense())
    D = np.asarray(csr_matrix(D).todense())

    [r,c] = np.shape(E)
    Y = np.zeros_like(E)

    if (len(C)==0 and len(D)==0):
        AEQ = True
        C = B
    else:
        AEQ = False

    if len(C)==0:
        AA, Z1 = scipy.linalg.schur(A)
        Q1 = Z1
        BB = np.eye(r)

    else:
        AA, BB, Q1, Z1 = scipy.linalg.qz(A,C)
        Q1 = np.transpose(Q1)


    if AEQ:
        T = AA
        R = BB 
        Q2 = Q1
        Z2 = Z1

    elif len(B) ==0:
        T, Z2 = scipy.linalg.schur(D)
        Q2 = Z2
        R = np.eye(c)

    else:
        T, R, Q2, Z2 = scipy.linalg.qz(D,B)
        Q2 = np.transpose(Q2)

    F = np.linalg.multi_dot([Q1,E,np.transpose(Q2)])

    AAY = np.zeros_like(E)
    BBY = np.zeros_like(E)

    k = c

    while(k>=0):




        if(k==0 or T[k,k-1]==0):


            jj = np.arange(k+1,c)
            print('jj',jj)
            rhs = F[:,k] - np.multiply(AAY[:,jj],np.transpose(R[k,jj]))-np.multiply(BBY[:,jj],np.transpose(T[k,jj]))

            tmp = (AA + (T[k,k]/R[k,k])*BB)
            rhs = rhs/R[k,k]
            Y[:,k] = np.linalg.lstst(tmp,rhs)

            AAY[:,k] = np.multiply(AA,Y[:,k])
            BBY[:,k] = np.multiply(BB, Y[:,k])

            k -= 1

        else:

            jj = np.arange(k+1,c)
            rhs1 = F[:,k-1] - np.multiply(AAY[:,jj],(R[k-1,jj]))-np.multiply(BBY[:,jj],np.transpose(T[k-1,jj]))
            rhs2 = F[:,k] - np.multiply(AAY[:,jj],np.transpose(R[k,jj]))-np.multiply(BBY[:,jj],np.transpose(T[k,jj]))


            BBM = [[np.multiply(R[k-1,k-1],AA)+ np.multiply(T[k-1,k-1],BB), np.multiply(R[k-1,k],AA)+np.multiply(T[k-1,k],BB)],[np.multiply(T[k,k-1],S), np.multiply(R[k,k], AA)+np.multiply(T[k,k],BB)]]

            idx = np.reshape([np.arange(0,r),np.arange(r+1,2*r)],(2*r,1))
            rhs = [rhs1,rhs2]
            UM = np.linalg.lstsq(BBM[idx,idx],rhs[idx])
            UM[idx] = UM

            Y[:,k-1:k] = np.reshape(UM, (r,2))
            AAY[:,k-1:k] = np.multiply(AA, Y[:,k-1:k])
            BBY[:,k-1:k] = np.multiply(BB, Y[:,k-1:k])

            k -= 2


    X = np.linalg.multi_dot(Z1,Y,np.transpose(Z2))

    return X


def get_matrix(name):
    r = int(input('Enter number of rows for {}'.format(name)))
    c = int(input('Enter number of columns for {}'.format(name)))
    print("Enter the entries in a single line (separated by space):")
    entries = list(map(int, input().split()))
    matrix = np.array(entries).reshape(r,c)
    return matrix



if __name__=="__main__":
    print('\t---GENERALIZED  SYLVESTER EQUATION SOLVER---\t\n')

    A = get_matrix('A')
    B = get_matrix('B')
    C = get_matrix('C')
    D = get_matrix('D')
    E = get_matrix('E')
    """
    #Dummy Example
    A = [[1 ,-1 ,1],[1 ,1 ,-1],[1, 1, 1]]
    B = np.array([]).reshape(0,0) 
    C = np.array([]).reshape(0,0)
    D = [[8 ,1 ,6],[3 ,5 ,7],[4 ,9 ,2]]
    E = np.eye(3)
    """

    X = sylvester_solver(A,B,C,D,E)
    print("==> The solution for the given Sylvester equation is : \n", X)









