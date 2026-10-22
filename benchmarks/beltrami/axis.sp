&physicslist
 Igeometry=3, Istellsym=1, Lfreebound=0,
 phiedge=1.0, Nfp=2, Nvol=1, Mpol=4, Ntor=0, Lrad=8,
 tflux=1.0, pflux=0.0, pscale=0.0, pressure=0.0,
 mu=1.1, Lconstraint=0,
 Rac(0)=6.0, Rbc(0,0)=6.0, Rbc(0,1)=1.0, Zbs(0,1)=1.0,
/
&numericlist
 Linitialize=1, Ndiscrete=2, Nquad=-1, iMpol=-4, iNtor=-4,
 Lsparse=0, Lsvdiota=0,
/
&locallist
 LBeltrami=4, Linitgues=1,
/
&globallist
 Lfindzero=0,
/
&diagnosticslist
 nPpts=0, nPtrj=0, Ngrid=16, Lcheck=1,
/
&screenlist
 Wjo00aa=T,
/
