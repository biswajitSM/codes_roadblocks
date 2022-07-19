import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import os

savepath = r"C:\Users\romanbarth\Workspace\PhD\Writing\Pradhan (2021) Draft - Roadblock study\Cell Reports"
mpl.rcParams['font.size'] = 14
mpl.rcParams['font.family'] = 'Arial'


particleRadii_discrete = (np.array([14, 23, 29, 40, 51, 125, 182])+15.7)/2
flowSpeed = np.linspace(0, 150, 151) # in um/s
particleRadii = np.linspace(5, 100, 151) # in nm

X, Y = np.meshgrid(flowSpeed, particleRadii)
F = 6*np.pi*1e-3*np.multiply(Y*1e-9, X*1e-6) * 1e12 # in pN

# map
f, ax = plt.subplots(figsize=(5,3.5))
im = ax.imshow(F, cmap=plt.get_cmap('Reds'))
xticks = np.linspace(0,len(flowSpeed)-1,6).astype(int)
ax.set_xticks(xticks)
ax.set_xticklabels(flowSpeed[xticks].astype(int))

yticks = []
blues = mpl.cm.get_cmap('Blues')
cmap_values = np.linspace(0.33, 1, len(particleRadii_discrete))
count = 0
for radius in particleRadii_discrete:
    yticks.append( np.argmin(np.abs(radius-particleRadii)) )
    ax.plot([0, len(flowSpeed)], [yticks[-1], yticks[-1]], '--', c=blues(cmap_values[count]))
    count += 1
ax.set_yticks(yticks)
ax.set_yticklabels((2*particleRadii_discrete).astype(int))

ax.set_xlabel(r'Flow speed [$\mu m/s$]')
ax.set_ylabel(r'Particle diameter [nm]')
ax.set_xlim(left=0, right=len(flowSpeed))
ax.invert_yaxis()
plt.colorbar(im, label="Drag force [pN]")
plt.tight_layout()
plt.savefig(os.path.join(savepath, 'dragForce_vs_particleRadius_map.png'))
plt.savefig(os.path.join(savepath, 'dragForce_vs_particleRadius_map.svg'))

# specific lines
f, ax = plt.subplots(figsize=(5,3.5))
count = 0
for radius in particleRadii_discrete:
    ind = np.argmin(np.abs(radius-particleRadii))
    ax.plot(flowSpeed, F[ind,:], '-', c=blues(cmap_values[count]))
    rightValue = F[ind,-1]
    ax.annotate(str((2*particleRadii_discrete[count]).astype(int))+' nm', (flowSpeed.max()*1.015, rightValue), ha='left', va='center', c=blues(cmap_values[count]))
    count += 1


ax.set_xlabel(r'Flow speed [$\mu m/s$]')
ax.set_ylabel(r'Drag force [pN]')
ax.set_xlim(left=0, right=flowSpeed.max()*1.2)
ax.set_ylim(bottom=0, top=0.3)
ax.fill([79-26, 79+26, 79+26, 79-26, 79-26], [0, 0, 0.3, 0.3, 0], c=(0.8, 0.8, 0.8))
plt.tight_layout()
plt.savefig(os.path.join(savepath, 'dragForce_vs_particleRadius.png'))
plt.savefig(os.path.join(savepath, 'dragForce_vs_particleRadius.svg'))
# plt.show()






# drag force acting on DNA
# from : https://backend.orbit.dtu.dk/ws/portalfiles/portal/123372399/PedersenPRE2016.pdf
gamma = 1e-3 # pN*s/(um^2) parallel and/or orthogonal drag cofficient per unit length
dgamma = 1e-3
v_min = 79-26
v_max = 79+26
v = 79 # um/s
dv = 26
kb = np.linspace(0, 48.5, 100)
L = kb * 1000 * 0.34 / 1000 # bps times nm/bp, expressed in microns
F_drag = gamma * v * L
F_drag_min = gamma * v_min * L
F_drag_max = gamma * v_max * L
F_drag_std = np.sqrt((v*L*dgamma)**2 + (gamma*L*dv)**2)


f, ax = plt.subplots(figsize=(5,3.5))
ax.set_xlabel(r'Loop length [kb]')
ax.set_ylabel(r'Drag force [pN]')
ax.plot(kb, F_drag, 'k', lw=2)
ax.fill(list(kb)+list(kb[::-1]), list(F_drag+F_drag_std/2)+list(F_drag[::-1]-F_drag_std[::-1]/2), c=(0.8, 0.8, 0.8))
ax.set_xlim(left=0, right=kb.max())
ax.set_ylim(bottom=0, top=1.3)

plt.tight_layout()
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength.png'))
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength.svg'))
ax.set_xlim(left=0, right=10)
ax.set_ylim(bottom=0, top=0.4)
plt.tight_layout()
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength_zoom.png'))
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength_zoom.svg'))
# plt.show()

def set_xaxis_color(ax, color):
    ax.xaxis.label.set_color(color)          #setting up Y-axis label color to blue
    ax.tick_params(axis='x', colors=color)  #setting up Y-axis tick color to black
    # ax.spines['bottom'].set_color(color) 

# plot both together
particleRadii = np.linspace(0, 100, 151) # in nm
F_drag_particle = 6*np.pi*1e-3*np.multiply(particleRadii*1e-9, v*1e-6) * 1e12 # in pN
F_drag_particle_std = 6*np.pi*1e-3*np.multiply(particleRadii*1e-9, dv*1e-6) * 1e12 # in pN

f, ax = plt.subplots(figsize=(5,3.5))
ax.set_xlabel(r'Loop length [kb]')
ax.set_ylabel(r'Drag force [pN]')
ax.plot(kb, F_drag, 'b', lw=2)
ax.fill(list(kb)+list(kb[::-1]), list(F_drag+F_drag_std/2)+list(F_drag[::-1]-F_drag_std[::-1]/2), \
    c=(0, 0, 1), alpha=0.25)
ax.set_xlim(left=0, right=10)
ax.set_ylim(bottom=0, top=0.4)
set_xaxis_color(ax, (0, 0, 1))

ax2 = ax.twiny()
ax2.set_xlabel(r'Particle diameter [nm]')
ax2.set_ylabel(r'Drag force [pN]')
ax2.plot(2*particleRadii, F_drag_particle, 'r', lw=2)
ax2.fill(list(2*particleRadii)+list(2*particleRadii[::-1]), \
    list(F_drag_particle+F_drag_particle_std/2)+list(F_drag_particle[::-1]-F_drag_particle_std[::-1]/2), \
        c=(1, 0, 0), alpha=0.25)
ax2.set_xlim(left=0, right=2*particleRadii.max())
# ax2.set_ylim(bottom=0, top=1.3)
set_xaxis_color(ax2, (1, 0, 0))

plt.tight_layout()
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength_particleSize.png'))
plt.savefig(os.path.join(savepath, 'dragForce_vs_loopLength_particleSize.svg'))

plt.show()

# experimental loop sizes
loopSizes = np.array( [22.8, 13.72, 21] ) # kb
particleRadius = 35.7/2
F_drag_particle = 6*np.pi*1e-3*np.multiply(particleRadius*1e-9, v*1e-6) * 1e12 # in pN
F_drag_particle_std = 6*np.pi*1e-3*np.multiply(particleRadius*1e-9, dv*1e-6) * 1e12 # in pN
gamma = 1e-3 # pN*s/(um^2) parallel and/or orthogonal drag cofficient per unit length
dgamma = 1e-3
v_min = 79-26
v_max = 79+26
v = 79 # um/s
dv = 26
# kb = np.linspace(0, 48.5, 100)
L = loopSizes * 1000 * 0.34 / 1000 # bps times nm/bp, expressed in microns
F_drag = gamma * v * L
F_drag_std = np.sqrt((v*L*dgamma)**2 + (gamma*L*dv)**2)
print('Force on loops: '+str(F_drag)+' +/- '+str(F_drag_std))
print('Force on particle: '+str(F_drag_particle)+' +/- '+str(F_drag_particle_std))
a=1



# particleRadii_discrete = (np.array([14, 23, 29, 40, 51, 125, 182])+15.7)/2
# for particleRadius in particleRadii_discrete:
#     f, ax = plt.subplots(figsize=(3,4))
#     # forces on loop 2
#     F_drag_particle = 6*np.pi*1e-3*np.multiply(particleRadius*1e-9, v*1e-6) * 1e12 # in pN
#     F_drag_particle_std = 6*np.pi*1e-3*np.multiply(particleRadius*1e-9, dv*1e-6) * 1e12 # in pN
#     kb = np.linspace(0, 10, 50)
#     F_drag_loop2 = L = loopSizes * 1000 * 0.34 / 1000 # bps times nm/bp, expressed in microns
#     F_drag = gamma * v * L
#     F_drag_loop2_std = np.sqrt((v*L*dgamma)**2 + (gamma*L*dv)**2)
#     F_loop2 = F_drag_loop2 + F_drag_particle
#     F_loop2_std = np.sqrt(F_drag_loop2_std**2+F_drag_particle_std**2)



