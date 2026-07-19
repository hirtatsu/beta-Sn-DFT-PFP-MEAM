import re, numpy as np

HA_PER_BOHR3_TO_GPA = 29421.02648438959
delta = 0.005

def last_stress(logfile):
    with open(logfile) as f:
        txt = f.read()
    # find all stress tensor blocks
    pat = re.compile(r'Stress tensor \(Hartree/bohr\^3\).*?\n\*+\s*\n\s*\n((?:\s*-?\d+\.\d+\s+-?\d+\.\d+\s+-?\d+\.\d+\s*\n){3})', re.S)
    matches = pat.findall(txt)
    if not matches:
        raise RuntimeError(f'no stress in {logfile}')
    rows = [list(map(float, r.split())) for r in matches[-1].strip().split('\n')]
    return np.array(rows)  # 3x3 in Ha/Bohr^3

tags = ['xx_p','xx_m','yy_p','yy_m','zz_p','zz_m','yz_p','yz_m','xz_p','xz_m','xy_p','xy_m']
stress = {}
for t in tags:
    s = last_stress(f'beta_Sn_{t}.log') * HA_PER_BOHR3_TO_GPA  # GPa
    stress[t] = s

# Voigt indices: 1=xx,2=yy,3=zz,4=yz,5=xz,6=xy
def voigt(s):
    return np.array([s[0,0], s[1,1], s[2,2], s[1,2], s[0,2], s[0,1]])

# Engineering strain magnitudes used:
# normal (xx_p etc.): eps_voigt[i]=+0.005
# shear  (xy_p etc.): eps_voigt[6]= 2*(delta/2)= +0.005  (since F_ij=F_ji=delta/2)
# So all +/- have magnitude 0.005 in Voigt. Central difference denom = 2*0.005 = 0.01
denom = 2*delta  # 0.01

# Column j of C: (sigma(+j) - sigma(-j)) / denom   with sigma as Voigt vector
pairs = [('xx_p','xx_m'),('yy_p','yy_m'),('zz_p','zz_m'),('yz_p','yz_m'),('xz_p','xz_m'),('xy_p','xy_m')]
C = np.zeros((6,6))
for j,(pp,mm) in enumerate(pairs):
    dsig = (voigt(stress[pp]) - voigt(stress[mm])) / denom  # GPa
    C[:,j] = dsig

# Symmetrize (elastic tensor should be symmetric)
Csym = 0.5*(C + C.T)

np.set_printoptions(precision=3, suppress=True, linewidth=120)
print('=== Raw stresses (GPa), last MD step ===')
for t in tags:
    sv = voigt(stress[t])
    print(f'  {t:6s}  σ(xx,yy,zz,yz,xz,xy) = {sv}')

print('\n=== C_ij from central differences (GPa), column = strain index ===')
print('    order: 1=xx 2=yy 3=zz 4=yz 5=xz 6=xy')
print(C)

print('\n=== Symmetrized C_ij (GPa) ===')
print(Csym)

# Tetragonal 4/mmm independent constants
C11 = 0.5*(Csym[0,0]+Csym[1,1])
C33 = Csym[2,2]
C12 = Csym[0,1]
C13 = 0.5*(Csym[0,2]+Csym[1,2])
C44 = 0.5*(Csym[3,3]+Csym[4,4])
C66 = Csym[5,5]

print('\n=== Tetragonal (I4_1/amd, 4/mmm) independent elastic constants (GPa) ===')
print(f'  C11 = {C11:8.2f}')
print(f'  C33 = {C33:8.2f}')
print(f'  C12 = {C12:8.2f}')
print(f'  C13 = {C13:8.2f}')
print(f'  C44 = {C44:8.2f}')
print(f'  C66 = {C66:8.2f}')

# Bulk modulus (Voigt-Reuss-Hill average for tetragonal)
B_V = (2*(C11+C12) + 4*C13 + C33)/9
# shear modulus Voigt
G_V = (2*C11 - C12 - 2*C13 + C33 + 6*C44 + 3*C66)/15
# Reuss averages use compliance tensor S = C^{-1}
S = np.linalg.inv(Csym)
B_R = 1.0/(S[0,0]+S[1,1]+S[2,2] + 2*(S[0,1]+S[0,2]+S[1,2]))
G_R = 15.0/(4*(S[0,0]+S[1,1]+S[2,2]) - 4*(S[0,1]+S[0,2]+S[1,2]) + 3*(S[3,3]+S[4,4]+S[5,5]))
B_H = 0.5*(B_V+B_R); G_H = 0.5*(G_V+G_R)
E_H = 9*B_H*G_H/(3*B_H+G_H)
nu_H = (3*B_H - 2*G_H)/(2*(3*B_H + G_H))
print('\n=== Polycrystalline averages (Voigt / Reuss / Hill) (GPa) ===')
print(f'  Bulk K:   V={B_V:7.2f}  R={B_R:7.2f}  H={B_H:7.2f}')
print(f'  Shear G:  V={G_V:7.2f}  R={G_R:7.2f}  H={G_H:7.2f}')
print(f'  Young E (Hill) = {E_H:7.2f} GPa,   Poisson ν (Hill) = {nu_H:.3f}')

# Stability (Born) for tetragonal 4/mmm
c1 = C11 > abs(C12)
c2 = 2*C13**2 < C33*(C11+C12)
c3 = C44 > 0
c4 = C66 > 0
print(f'\nBorn stability checks: C11>|C12|:{c1}  2C13^2<C33(C11+C12):{c2}  C44>0:{c3}  C66>0:{c4}')
