#### Equivalenve of modules between pL and L that are principally polarized ####
##                                                                            ##
## The input was obtained using Script 1                                      ##
## The output is a pair [T1, alpha] such that                                 ##
##      T1 == T0                                                              ##
##      alphatp*TM*alphaTP.conjugate_transpose() == T0                        ##
##                                                                            ##
################################################################################

### Tools

def transformation(E,X,i,j,k,lst=None):
    E = E.with_added_multiple_of_row(i,j,k);
    E = E.with_added_multiple_of_column(i,j,k.conjugate());
    X = X.with_added_multiple_of_row(i,j,k);
    if lst != None:
        lst.append(['T',i,j,k])
        return E, X, lst
    return E, X
def swap(E,X,i,j,lst=None):
    E.swap_rows(i,j);
    E.swap_columns(i,j);
    X.swap_rows(i,j);
    if lst != None:
        lst.append(['S',i,j])
        return E, X, lst
    return E, X
def rescale(E,X,i,k,lst=None):
    E = E.with_rescaled_row(i,k)
    E = E.with_rescaled_col(i,k.conjugate())
    X = X.with_rescaled_row(i,k)
    if lst != None:
        lst.append(['R',i,k])
        return E, X, lst
    return E,X
def chain_of_trans(E,X,lst):
    for L in lst:
        if L[0] == 'T':
            E,X = transformation(E,X,L[1],L[2],L[3])
        elif L[0] == 'S':
            E,X = swap(E,X,L[1],L[2])
        elif L[0] == 'R':
            E,X = rescale(E,X,L[1],L[2])
        else:
            raise TypeError('Char should be T, S or R')
    return E, X


### Setting

k.<z> = CyclotomicField(5)

### Computation of T0 (see Example 3.4.2)

T = 1/5*matrix(k,[[z-z^4, 1-z^4, 1-z^4],[-1+z,z- z^4, 1-z^4],[-1+z,-1+z,z-z^4]])
X = identity_matrix(k,3)
(T0, X) = transformation(T, X, 1, 0, -(z-1)/(z - z^4))
(T0, X) = transformation(T0, X, 2, 0, -(z-1)/(z - z^4))
(T0, X) = transformation(T0, X, 2, 1, -1)
(T0, X) = transformation(T0, X, 2, 1, -z^2)
(T0, X) = rescale(T0, X, 2, (z+1))
(T0, X) = swap(T0, X, 1, 2)

### List of cases obtained with Script 1, and computation of equivalence with [M0,T0]

#Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [1 0 0]
# [0 1 0]
# [0 0 1]

# This is M0

 #############################################################
 
# Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5 -1/5*z^3 - 2/5*z^2 + 2/5*z - 9/5                                0]
# [                               0                            z - 1                                0]
# [                               0                                0                                1]

gammatp = matrix(k, [(-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5, -1/5*z^3 - 2/5*z^2 + 2/5*z - 9/5, 0), (0, z - 1, 0), (0, 0, 1)])
TM = gammatp*T0*gammatp.conjugate_transpose()
alphatp = identity_matrix(k,3)
(T1, alphatp) = transformation(TM, alphatp, 0, 1, -1)
(T1, alphatp) = transformation(T1, alphatp, 1, 0, -(-1/5*z^2 + 1/5)/(1/5*z^3 + 1/5*z^2 + 2/5*z + 1/5))
(T1, alphatp) = rescale(T1, alphatp, 1, (z + 1))
assert T1 == T0
[T1, alphatp]

 #############################################################
 
# Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5 -1/5*z^3 - 2/5*z^2 + 2/5*z + 6/5                                0]
# [                               0                            z - 1                                0]
# [                               0                                0                                1]

gammatp = matrix(k,[(-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5, -1/5*z^3 - 2/5*z^2 + 2/5*z + 6/5, 0),  (0, z - 1, 0), (0, 0, 1)])
TM = gammatp*T0*gammatp.conjugate_transpose()
alphatp = identity_matrix(k,3)
(T1, alphatp) = transformation(TM, alphatp, 0, 1, z+1)
(T1, alphatp) = transformation(T1, alphatp, 1, 0, -(1/5*z^3 + 2/5*z^2 + 1/5*z + 1/5)/(1/5*z^3 + 1/5*z^2 + 2/5*z + 1/5))
(T1, alphatp) = rescale(T1, alphatp, 1, (z + 1))
print T1 == T0
[T1, alphatp]
 
 #############################################################
 
# Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5 -1/5*z^3 - 2/5*z^2 + 2/5*z - 4/5                                0]
# [                               0                            z - 1                                0]
# [                               0                                0                                1]

gammatp = matrix(k,[(-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5, -1/5*z^3 - 2/5*z^2 + 2/5*z - 4/5, 0),  (0, z - 1, 0),  (0, 0, 1)])
TM = gammatp*T0*gammatp.conjugate_transpose()
alphatp = identity_matrix(k,3)
(T1, alphatp) = transformation(TM, alphatp, 1, 0, -1)
(T1, alphatp) = rescale(T1, alphatp, 0, (z + 1))
print T1 == T0
[T1, alphatp]

 
 #############################################################
 
# Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [ -3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5 -1/5*z^3 - 2/5*z^2 + 2/5*z - 14/5                                 0]
# [                                0                             z - 1                                 0]
# [                                0                                 0                                 1]

gammatp = matrix(k,[(-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5, -1/5*z^3 - 2/5*z^2 + 2/5*z - 14/5, 0),  (0, z - 1, 0),  (0, 0, 1)])
TM = gammatp*T0*gammatp.conjugate_transpose()
alphatp = identity_matrix(k,3)
(T1, alphatp) = transformation(TM, alphatp, 0, 1, -2)
(T1, alphatp) = transformation(T1, alphatp, 0, 1, z^4+1)
(T1, alphatp) = transformation(T1, alphatp, 1, 0, -(-1/5*z^3 + 1/5*z)/(1/5*z^3 + 1/5*z^2 + 2/5*z + 1/5))
(T1, alphatp) = rescale(T1, alphatp, 1, (z^4 + 1))
print T1 == T0
[T1, alphatp]
 
 #############################################################
 
# Free module of degree 3 and rank 3 over Maximal Order in Cyclotomic Field of order 5 and degree 4
# Echelon basis matrix:
# [-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5 -1/5*z^3 - 2/5*z^2 + 2/5*z + 1/5                                0]
# [                               0                            z - 1                                0]
# [                               0                                0                                1]

gammatp = matrix(k, [(-3/5*z^3 - 1/5*z^2 + 1/5*z - 2/5, -1/5*z^3 - 2/5*z^2 + 2/5*z + 1/5, 0),  (0, z - 1, 0),  (0, 0, 1)])
TM = gammatp*T0*gammatp.conjugate_transpose()
alphatp = identity_matrix(k,3)
(T1, alphatp) = transformation(TM, alphatp, 1, 0, -(1/5*z^3 - 1/5)/(1/5*z^3 + 1/5*z^2 + 2/5*z + 1/5))
(T1, alphatp) = rescale(T1, alphatp, 1, (z^4 + 1))
print T1 == T0
[T1, alphatp]
