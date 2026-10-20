"""Independent high-precision oracle; optional regeneration, never a runtime dependency.

Run from the repository root: python3 test/reference/pp_fusion_eft_oracle.py
Uses E1 (not production's Coulomb integral), unreduced C0/C2 matching, and
adaptive beta-decay phase-space quadrature (not production's Gauss rule).
Sources/equations and units: docs/pp-fusion-eft.md.
"""
from pathlib import Path
import mpmath as m

m.mp.dps = 70
hc, hb = m.mpf('197.3269804'), m.mpf('6.582119569e-22')
mp, mn, me = map(m.mpf, ['938.27208943', '939.56542194', '0.51099895069'])
alpha, ga, ft = map(m.mpf, ['0.0072973525643', '1.2695', '3071.4'])
bd, a, r, rho = map(m.mpf, ['2.224575', '-7.8063', '2.794', '1.764'])
gamma = m.sqrt(2*mp*mn*bd/(mp+mn))/hc
mass, mu = mp/hc, m.mpf(140)/hc
chi = alpha*mass/gamma
jc = 1-chi*m.exp(chi)*m.e1(chi)
c0p = 4*m.pi/mass/(1/a-mu+alpha*mass*(m.log(mu*m.sqrt(m.pi)/(alpha*mass))+1-3*m.euler/2))
c2p = mass*r*c0p**2/(8*m.pi)
c2d = 2*m.pi*rho/(mass*(mu-gamma)**2)
w0 = (mp-mn+bd)/me

def phase(w):
    if w <= 1 or w >= w0:
        return m.mpf(0)
    p = m.sqrt(w*w-1)
    z = 2*m.pi*alpha*w/p
    return w*p*(w0-w)**2*z/m.expm1(z)

f = m.mpf('1.0162')*m.quad(phase, [1, (1+w0)/2, w0])
prefactor = 6*m.pi**2*mp*alpha*m.log(2)*hb/hc*ga**2*f/(gamma**3*ft*100)
# Independent finite-energy spectral integral from Chen--Liu--Yu (7)-(11).
# A principal-value subtraction removes the real pole. Imaginary Coulomb
# penetration terms are flat at E=0 and do not contribute to Taylor slopes.
def spectral(p2):
    p = m.sqrt(p2)
    def weight(k):
        x = m.pi*alpha*mass/k
        return x/m.expm1(x)*m.exp(alpha*mass/k*m.atan(k/gamma))
    atpole = weight(p)*p2/(p2+gamma**2)
    def integrand(k):
        if k == p:
            return m.diff(lambda t: weight(t)*t*t/(t*t+gamma**2), p)/(2*p)
        return (weight(k)*k*k/(k*k+gamma**2)-atpole)/(k*k-p2)
    return 2*gamma/m.pi*m.quad(integrand, [0, p/2, 1.5*p, gamma, 1, m.inf])

p2step = m.mpf('2e-9')
spectra = [spectral(p2step), spectral(p2step/2)]

def energy_rates(p2, jenergy, lb):
    p = m.sqrt(p2)
    eta = alpha*mass/(2*p)
    h = m.re(m.digamma(1+1j*eta))-m.log(eta)
    den = 1/a+alpha*mass*h
    direct = gamma**2/(gamma**2+p2)*m.exp(alpha*mass/p*m.atan(p/gamma))
    aa = direct-gamma*jenergy/den
    bb = -gamma**2*lb/den-gamma*jenergy*r*p2/(2*den**2)
    cc = -gamma**2*lb*r*p2/(2*den**2)-gamma*jenergy*(r*p2/2)**2/den**3
    z = gamma*rho
    c0, c1, c2 = aa, aa*z/2+bb, 3*aa*z*z/8+bb*z/2+cc
    res = (direct-(gamma*jenergy+gamma**2*lb)/(den-r*p2/2))/m.sqrt(1-z)
    return [c0*c0, c0*c0+2*c0*c1, c0*c0+2*c0*c1+c1*c1+2*c0*c2, res*res]

def phase_energy(energy):
    endpoint = (mp-mn+bd+energy)/me
    def integrand(w):
        if w <= 1 or w >= endpoint:
            return m.mpf(0)
        p = m.sqrt(w*w-1)
        x = 2*m.pi*alpha*w/p
        return w*p*(endpoint-w)**2*x/m.expm1(x)
    return m.mpf('1.0162')*m.quad(integrand, [1, (1+endpoint)/2, endpoint])
phase_slope = m.diff(phase_energy, 0)/f
rows = []
for la in map(m.mpf, ['0', '3.9']):
    lb = -(mu-gamma)/(mass*c0p)*(la/ga-mass*(c2p+c2d)/2)
    aa, bb, z = m.exp(chi)-gamma*a*jc, -gamma**2*a*lb, gamma*rho
    c = [aa, aa*z/2+bb, 3*aa*z*z/8+bb*z/2]
    amps = [c[0], c[0]+c[1], sum(c), (aa+bb)/m.sqrt(1-z)]
    rates = [c[0]**2, c[0]**2+2*c[0]*c[1],
             c[0]**2+2*c[0]*c[1]+c[1]**2+2*c[0]*c[2], amps[3]**2]
    rates_at_step = energy_rates(p2step, spectra[0], lb)
    rates_at_half = energy_rates(p2step/2, spectra[1], lb)
    slopes = [mp/hc**2*((4*half-3*base-full)/p2step)/base+phase_slope
              for base, full, half in zip(rates, rates_at_step, rates_at_half)]
    rows.append(amps+rates+[prefactor*x for x in rates]+slopes)
path = Path('test/reference/pp_fusion_eft.dat')
with path.open('w') as out:
    for row in rows:
        out.write(' '.join(m.nstr(x, 40) for x in row)+'\n')
print(path)
for row in rows:
    print('Lambda2:', *(float(x) for x in row[4:8]))
    print('S [1e-25 MeV b]:', *(float(x/m.mpf('1e-25')) for x in row[8:12]))
    print('slope [MeV^-1]:', *(float(x) for x in row[12:]))
