K.<z> = CyclotomicField(5)
OK = K.OK()
p = z - 1
gr = (1 + K(5).sqrt())/2
S0 = diagonal_matrix(K,[gr,gr,-gr-1])
M0 = OK^3            
N_lst = [[p,0,0],[0,p,0],[0,0,p],[1,2,0]]
L_lst = [[1,0,0],[0,1,0],[0,0,1],[p^-1,2*p^-1,0]]
N = M0.span(N_lst)
L = M0.span(L_lst)

@parallel
def check_lst(lst):
    bs = [[el*p^-1 for el in N.linear_combination_of_basis(v.change_ring(ZZ)).list()] for v in lst]
    M_lst = [[p,0,0],[0,p,0],[0,0,p],[1,2,0]]+bs
    M = M0.span(M_lst)
    if M.matrix().det()*OK != 1*OK:
        return False
    mu0M = []
    M_bs = M.basis()
    for el1 in M_bs:
        for el2 in M_bs:
            v1 = (K^3)(el1)
            v2 = (K^3)(el2)
            mu0M.append(v1*S0*v2.conjugate())
    if OK.fractional_ideal(mu0M) != 1*OK:
        return False
    return True

final = sorted(list(check_lst(Combinations(GF(5)^3,3).list())))
mods = []
for el in final:
    if el[1]:
        bs = [[el*p^-1 for el in N.linear_combination_of_basis(v.change_ring(ZZ)).list()] for v in el[0][0][0]]
        M_lst = [[p,0,0],[0,p,0],[0,0,p],[1,2,0]]+bs
        M = M0.span(M_lst)
        if M not in mods:
            mods.append(M)
print mods.cardinality()
