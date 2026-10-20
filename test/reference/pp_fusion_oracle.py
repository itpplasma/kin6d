"""Optional independent ODE/BVP oracle, never imported by production.
Run from kin6d root: python test/reference/pp_fusion_oracle.py.
Uses scipy 1.18.1 / numpy 2.5.4 / mpmath 1.4.1; no FEM-DVR or generated leaves.
Writes a small frozen numeric table for the native Fortran test.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_bvp, solve_ivp, quad
import mpmath as m
m.mp.dps = 45
hc, mp, mn, me, alpha = 197.3269804, 938.27208943, 939.56542194, .51099895069, .0072973525643
kin = hc*hc/mp
R = 100.

def potential(r):
    x = .7*r
    e = np.exp(-x)
    singlet = (-10.463*e-1650.6*e**4+6484.2*e**7)/x
    vc = (-10.463*e+105.468*e**2-3187.8*e**4+9924.3*e**6)/x
    vt = -10.463*((1+3/x+3/x**2)*e-(12/x+3/x**2)*e**4)/x + (351.77*e**4-1673.5*e**6)/x
    ls = (708.91*e**4-2713.1*e**6)/x
    return singlet, vc, np.sqrt(8)*vt, vc-2*vt-3*ls+6*41.47/r**2

def rhs(r,y,p):
    _,ss,sd,dd = potential(r)
    return np.array([y[1], ((ss+p[0])*y[0]+sd*y[2])/41.47,
                     y[3], (sd*y[0]+(dd+p[0])*y[2])/41.47,
                     y[0]**2+y[2]**2])

def bc(a,b,p):
    return np.array([a[0],a[2],b[0],b[2],a[4],b[4]-1])

r = np.geomspace(.0001,R,1400)
u = (np.exp(-.2316*r)-np.exp(-1.5*r))*.9
w = .05*r**3*np.exp(-r)
y = np.array([u,np.gradient(u,r),w,np.gradient(w,r),1-np.exp(-r/3)])
sol = solve_bvp(rhs,bc,r,y,p=[2.2246],tol=3e-9,max_nodes=50000)
assert sol.success, sol.message
B = sol.p[0]
gamma = np.sqrt(B/41.47)
PD = quad(lambda r: sol.sol(r)[2]**2,.0001,R,epsabs=2e-10)[0]

def scatter(E):
    k = np.sqrt(E/kin)
    eta = alpha*hc/(2*kin*k)
    def rhs(r,y):
        v = potential(r)[0]+alpha*hc/r
        return [y[1],(v-E)/kin*y[0]]
    # Overall scale arbitrary; origin cutoff checked by halving separately.
    scat = solve_ivp(rhs,[.0001,R],[0.,1.],method='DOP853',
                     rtol=8e-12,atol=1e-12,dense_output=True)
    assert scat.success
    val, der = scat.y[:,-1]
    rho = m.mpf(str(k*R)); et = m.mpf(str(eta))
    F, G = m.coulombf(0,et,rho),m.coulombg(0,et,rho)
    Fp = m.diff(lambda x:m.coulombf(0,et,x),rho)
    Gp = m.diff(lambda x:m.coulombg(0,et,x),rho)
    A = float(G*der/k-Gp*val)
    C = float(Fp*val-F*der/k)
    norm = np.hypot(A,C)
    C2 = float(2*m.pi*et/m.expm1(2*m.pi*et))
    h = float(m.re(m.digamma(1+1j*et))-m.log(et))
    K = C2*k*A/C+2*k*eta*h
    overlap = quad(lambda r:sol.sol(r)[0]*scat.sol(r)[0]/(norm*np.sqrt(C2)*k),
                   .0001,R,epsabs=1e-12,epsrel=2e-10)[0]
    lam2 = gamma**3*overlap**2/2
    return K,lam2,np.arctan2(C,A),lambda r:scat.sol(r)[0]/norm

rows = [scatter(E) for E in [.0001,.0002,.0004]]
K,lam = (8*np.array(rows[0][:2])-6*np.array(rows[1][:2])+np.array(rows[2][:2]))/3
_,_,delta,upp = scatter(.001)
endpoint=(mp-mn+B)/me
f = quad(lambda x:x*np.sqrt(x*x-1)*(endpoint-x)**2 *
         (2*np.pi*alpha*x/np.sqrt(x*x-1))/np.expm1(2*np.pi*alpha*x/np.sqrt(x*x-1)),
         1,endpoint,epsabs=1e-12)[0]*1.0162
S = 6*np.pi**2*mp*alpha*np.log(2)*6.582119569e-22/hc*1.2695**2*lam*f/(gamma**3*3071.4*100)
out = Path('test/reference/pp_fusion_ode.dat')
with out.open('w') as stream:
    stream.write(' '.join(f'{v:.16e}' for v in [B,PD,-1/K,lam,S,delta,f])+'\n')
    for r in [.1,.5,1.,2.,5.,10.,20.,40.]:
        stream.write(' '.join(f'{v:.16e}' for v in
                             [r,sol.sol(r)[0],sol.sol(r)[2],upp(r)])+'\n')
print(out, B, PD, -1/K, lam, S, delta, f)
