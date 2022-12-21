# %%
from colossus.cosmology import cosmology
from colossus.lss import mass_function

# %%
from __future__ import print_function 
import matplotlib.pyplot as plt
import numpy as np
%matplotlib inline

# %%
from colossus.cosmology import cosmology
cosmology.setCosmology('planck15');


# %%
from colossus.lss import mass_function

cosmology.setCosmology('planck15');

z = [0.0, 1.0, 2.0, 4.0]
M = 10**np.arange(11.0, 15.5, 0.1)

plt.figure()
plt.xlabel('M200m')
plt.ylabel('Tinker dn/dln(M)')
plt.loglog()
plt.xlim(1E11, 4E15)
plt.ylim(1E-7, 1E-1)
for i in range(len(z)):
    mfunc = mass_function.massFunction(M, z[i], mdef = '200m', 
        model = 'tinker08', q_out = 'dndlnM')
    plt.plot(M, mfunc, '-', label = 'P15, z = %.1f' % (z[i]))




cosmology.setCosmology('planck18');

z = [0.0, 1.0, 2.0, 4.0]
M = 10**np.arange(11.0, 15.5, 0.1)

for i in range(len(z)):
    mfunc = mass_function.massFunction(M, z[i], mdef = '200c', 
        model = 'tinker08', q_out = 'dndlnM')
    plt.plot(M, mfunc, '--', label = 'P18, z = %.1f' % (z[i]))

plt.legend(ncol=2);

# %%
from colossus.halo import mass_adv
Mvir = 1E12
z = 2.0
M200m, R200m, c200m = mass_adv.changeMassDefinitionCModel(Mvir, z, '200c', '200m')

print(M200m / Mvir)

# %%
