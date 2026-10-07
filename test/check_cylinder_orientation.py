#!/usr/bin/env python3
"""Independent Cartesian C0 orientation/curl diagnostic; not a CAS proof.

Straight axis points along physical +phi at phi=0 (Cartesian +y).
RH cylindrical angle rotates from +x towards -z; canonical poloidal angle
rotates from +x towards +z. No native outputs are modified.
"""
import json
import numpy as np

R0 = 6.2
errors = {key: 0.0 for key in ('basis', 'field_pullback', 'q', 'curl_relative')}
for phi in (0.0, 0.7, 2.1):
    eR = np.array([np.cos(phi), np.sin(phi), 0.0])
    ephi = np.array([-np.sin(phi), np.cos(phi), 0.0])
    eZ = np.array([0.0, 0.0, 1.0])
    for theta in (0.0, 0.4, 1.8, 3.2):
        er = np.cos(theta)*eR + np.sin(theta)*eZ
        etc = -np.sin(theta)*eR + np.cos(theta)*eZ
        etRH = np.cross(ephi, er)
        errors['basis'] = max(errors['basis'], np.linalg.norm(etRH+etc),
                              np.linalg.norm(np.cross(er, etc)+ephi))
        for qRH in (-1.5, 1.5):
            k = 1.0/(qRH*R0)
            for B0 in (-5.3, 5.3):
                r = 0.3
                bz = B0/(1.0+(k*r)**2)
                btRH = k*r*bz
                field_native = bz*ephi + btRH*etRH
                dpsi_c = -R0*btRH
                field_canonical = bz*ephi + dpsi_c*np.cross(er,ephi)/R0
                errors['field_pullback'] = max(errors['field_pullback'],
                                               np.linalg.norm(field_native-field_canonical))
                qc = r*bz/(R0*(-btRH))
                errors['q'] = max(errors['q'], abs(qc+qRH))

# Direct Cartesian field: (Bx,By,Bz)=(k*z*bz,bz,-k*x*bz).
# FD4 differentiates that field independently of cylindrical curl formulas.
for qRH in (-1.5, 1.5):
    k = 1.0/(qRH*R0)
    for B0 in (-5.3, 5.3):
        def field(xyz):
            x, _, z = xyz
            bz = B0/(1.0+k*k*(x*x+z*z))
            return np.array([k*z*bz, bz, -k*x*bz])
        for point in ([0.1,0.0,0.2], [0.4,0.3,-0.2]):
            point = np.array(point, dtype=float)
            h = 1.e-3
            jac = np.empty((3,3))
            for j in range(3):
                step = np.eye(3)[j]*h
                jac[:,j] = (field(point-2*step)-8*field(point-step)
                            +8*field(point+step)-field(point+2*step))/(12*h)
            curl = np.array([jac[2,1]-jac[1,2],jac[0,2]-jac[2,0],jac[1,0]-jac[0,1]])
            alpha = 2*k/(1+k*k*(point[0]**2+point[2]**2))
            exact = alpha*field(point)
            errors['curl_relative'] = max(errors['curl_relative'],
                                          np.linalg.norm(curl-exact)/np.linalg.norm(exact))
assert errors['basis'] < 2.e-15
assert errors['field_pullback'] < 2.e-15
assert errors['q'] < 2.e-15
assert errors['curl_relative'] < 1.e-9
print(json.dumps({'scope':'Cartesian numerical orientation/curl check, not proof',
                  'errors':errors}, indent=2))
