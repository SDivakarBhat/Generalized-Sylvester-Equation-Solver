### Generalized-Sylvester-Equation-Solver
A python based solver for Generalized Sylvester Equation.

An implementation of the Generalized Sylvester Equation solver.

            AXB + CXD = E
using the Bartels-Stewart algorithm.

## References

[1]  Chu, K.W.E., 1987. The solution of the matrix equations AXB - CXD= E AND (YA - DZ, YC -BZ)=(E, F). Linear Algebra and its applications, 93, pp.93-105.

[2]  Hern ́andez,  V.  and  Gass ́o,  M.,  1989.  Explicit  solution  of  the  matrix  equation  AXB  -  CXD=  E.Linear Algebra and its applications, 121, pp.333-344.

[3]  Bartels, R.H. and Stewart, G.W., 1972. Solution of the matrix equation AX + XB= C [F4]. Com-munications of the ACM, 15(9), pp.820-826.

[4]  Gardiner,  J.D.,  Laub,  A.J.,  Amato,  J.J. and Moler,  C.B.,  1992. Solution of the Sylvester matrixequationAXBT-CXDT=  E.  ACM  Transactions  on  Mathematical  Software  (TOMS),  18(2),pp.223-231.

This work is inspired from the MATLAB implementation of Lyapunov solver freeLYAP by Nick Hale, Alex Townsend, and Heather Wilber available at https://github.com/ajt60gaibb/freeLYAP.
Dummy Examples

# Example 1:
    A = [[1 ,-1 ,1],[1 ,1 ,-1],[1, 1, 1]]
    
    B = np.array([]).reshape(0,0) 
    
    C = np.array([]).reshape(0,0)
    
    D = [[8 ,1 ,6],[3 ,5 ,7],[4 ,9 ,2]]
    
    E = np.eye(3)
