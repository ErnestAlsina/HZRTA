import numpy as np
import matplotlib.pyplot as plt
#
savfil = ["HZRTA_CaI_PRD_0.0g_90d_0d.npz","HZRTA_CaI_PRD_5.0g_90d_0d.npz",\
          "HZRTA_CaI_PRD_10.0g_90d_0d.npz","HZRTA_CaI_PRD_20.0g_90d_0d.npz",\
          "HZRTA_CaI_PRD_50.0g_90d_0d.npz","HZRTA_CaI_PRD_100.0g_90d_0d.npz"]
nfils = np.shape(savfil)[0]
tgfil = [0,0,0,0,0,0]
colrs = ["black","purple","blue","green","goldenrod","red"]
strlw = [2.25,2.0,1.75,1.5,1.25,1.0]
#
wlmin = 4225.9
wlmax = 4227.6
#
plt.subplot(221)
ax1 = plt.gca()
ax1.get_xaxis().get_major_formatter().set_useOffset(False)
#ax1.set_xlim(1205.66,1225.66)
ax1.set_xlim(wlmin,wlmax)
ax1.set_yscale("log")
#
plt.subplot(222)
ax2 = plt.gca()
ax2.get_xaxis().get_major_formatter().set_useOffset(False)
ax2.set_ylim(0.0,5.0)
ax2.set_xlim(wlmin,wlmax)
#
plt.subplot(223)
ax3 = plt.gca()
ax3.get_xaxis().get_major_formatter().set_useOffset(False)
ax3.set_ylim(0.0,2.5)
ax3.set_xlim(wlmin,wlmax)
#
plt.subplot(224)
ax4 = plt.gca()
ax4.get_xaxis().get_major_formatter().set_useOffset(False)
ax4.set_xlim(wlmin,wlmax)
#
for vv,stv in enumerate(savfil):
 tg = tgfil[vv]
 clr = colrs[vv]
 lwv = strlw[vv]
 sav = np.load(stv)
 wl = sav["wavl"]
 si = sav["stok"][0,:,tg]
 sq = 100.0*sav["stok"][1,:,tg]/sav["stok"][0,:,tg]
 su = 100.0*sav["stok"][2,:,tg]/sav["stok"][0,:,tg]
 sv = 100.0*sav["stok"][3,:,tg]/sav["stok"][0,:,tg]
 ax1.plot(wl,si,color=clr,lw=lwv)
 ax2.plot(wl,sq,color=clr,lw=lwv)
 ax3.plot(wl,su,color=clr,lw=lwv)
 ax4.plot(wl,sv,color=clr,lw=lwv)
#
plt.show()
