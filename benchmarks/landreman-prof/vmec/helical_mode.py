import numpy as np,sys
a=np.loadtxt(sys.argv[1]); e=a[:,2].reshape(3,24,24)[1]
np.set_printoptions(linewidth=250,precision=1,suppress=True)
print((e*1000)[::2,::2])
f=np.fft.fft2(e); print("top modes (pol,tor):",[ (i,j,round(abs(f[i,j])/576*1000,2)) for i,j in zip(*np.unravel_index(np.argsort(-abs(f).ravel())[:8],f.shape))])
