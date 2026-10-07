import json
import numpy as np
R0=6.2;a=.62;F=32.86;q0=1.5;C=F/(8*q0*R0**3)
t=(np.arange(32768)+.5)*2*np.pi/32768
def flux(s):
 R=np.sqrt(R0**2+2*R0*a*np.sqrt(s)*np.cos(t))
 return 2*np.pi*F*R0**2*a*a*s*np.mean(np.sin(t)**2/R**3)
out=[]
for s in [.001,.1,.5,1.]:
 r=a*np.sqrt(s);R=np.sqrt(R0**2+2*R0*r*np.cos(t));Z=R0*r*np.sin(t)/R
 Rr=R0*np.cos(t)/R;Rt=-R0*r*np.sin(t)/R
 Zr=R0*np.sin(t)/R-R0*r*np.sin(t)*Rr/R**2
 Zt=R0*r*np.cos(t)/R-R0*r*np.sin(t)*Rt/R**2
 jac=Rr*Zt-Rt*Zr
 e_jac=np.max(np.abs(jac-R0**2*r/R**2))
 psi=C*((R*R-R0*R0)**2+4*R*R*Z*Z)
 e_psi=np.max(np.abs(psi-4*C*R0**2*r*r))
 q=F/(8*C)*np.mean(1/R**3)
 ds=1e-6
 q_flux=(flux(s+ds)-flux(s-ds))/(2*ds)/(4*C*R0**2*a*a)/(2*np.pi)
 out.append(dict(s_pol=s,jacobian_absolute_error_m=e_jac,psi_absolute_error_Wb_per_rad=e_psi,q=q,q_flux_derivative=q_flux,q_relative_error=abs(q_flux/q-1)))
assert max(x['jacobian_absolute_error_m'] for x in out)<1e-13
assert max(x['psi_absolute_error_Wb_per_rad'] for x in out)<1e-13
assert max(x['q_relative_error'] for x in out)<1e-9
print(json.dumps(out,indent=2))
