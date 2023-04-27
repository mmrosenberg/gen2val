
import numpy as np
import matplotlib.pyplot as plt

from selection_output.optimize_selection_cuts_output_arrays import pCuts as pCuts_2d
from selection_output.optimize_selection_cuts_output_arrays import sCuts as sCuts_2d
from selection_output.optimize_selection_cuts_output_arrays import purity as purity_2d
from selection_output.optimize_selection_cuts_output_arrays import efficiency as efficiency_2d
from selection_output.optimize_selection_cuts_output_arrays import purXeff as purXeff_2d

from selection_output.optimize_selection_cuts_output_arrays_noPCut import pCuts as pCuts_1d
from selection_output.optimize_selection_cuts_output_arrays_noPCut import sCuts as sCuts_1d
from selection_output.optimize_selection_cuts_output_arrays_noPCut import purity as purity_1d
from selection_output.optimize_selection_cuts_output_arrays_noPCut import efficiency as efficiency_1d
from selection_output.optimize_selection_cuts_output_arrays_noPCut import purXeff as purXeff_1d

from selection_output.optimize_selection_cuts_output_arrays_piPhCuts import piCuts as piCuts_s2d
from selection_output.optimize_selection_cuts_output_arrays_piPhCuts import phCuts as phCuts_s2d
from selection_output.optimize_selection_cuts_output_arrays_piPhCuts import purity as purity_s2d
from selection_output.optimize_selection_cuts_output_arrays_piPhCuts import efficiency as efficiency_s2d
from selection_output.optimize_selection_cuts_output_arrays_piPhCuts import purXeff as purXeff_s2d

pCuts_2d_optP = []
sCuts_2d_optP = []
purity_2d_optP = []
efficiency_2d_optP = []
purXeff_2d_optP = []

pCuts_2d_optS = []
sCuts_2d_optS = []
purity_2d_optS = []
efficiency_2d_optS = []
purXeff_2d_optS = []

for i, pCut in enumerate(pCuts_2d):
  if abs(pCut - 0.77) < 1e-6:
    pCuts_2d_optP.append(pCuts_2d[i])
    sCuts_2d_optP.append(sCuts_2d[i])
    purity_2d_optP.append(purity_2d[i])
    efficiency_2d_optP.append(efficiency_2d[i])
    purXeff_2d_optP.append(purXeff_2d[i])

for i, sCut in enumerate(sCuts_2d):
  if abs(sCut - 7.1) < 1e-6:
    pCuts_2d_optS.append(pCuts_2d[i])
    sCuts_2d_optS.append(sCuts_2d[i])
    purity_2d_optS.append(purity_2d[i])
    efficiency_2d_optS.append(efficiency_2d[i])
    purXeff_2d_optS.append(purXeff_2d[i])


piCuts_s2d_optPi = []
phCuts_s2d_optPi = []
purity_s2d_optPi = []
efficiency_s2d_optPi = []
purXeff_s2d_optPi = []

piCuts_s2d_optPh = []
phCuts_s2d_optPh = []
purity_s2d_optPh = []
efficiency_s2d_optPh = []
purXeff_s2d_optPh = []

for i, piCut in enumerate(piCuts_s2d):
  if abs(piCut - 8.9) < 1e-6:
    piCuts_s2d_optPi.append(piCuts_s2d[i])
    phCuts_s2d_optPi.append(phCuts_s2d[i])
    purity_s2d_optPi.append(purity_s2d[i])
    efficiency_s2d_optPi.append(efficiency_s2d[i])
    purXeff_s2d_optPi.append(purXeff_s2d[i])

for i, phCut in enumerate(phCuts_s2d):
  if abs(phCut - 3.1) < 1e-6:
    piCuts_s2d_optPh.append(piCuts_s2d[i])
    phCuts_s2d_optPh.append(phCuts_s2d[i])
    purity_s2d_optPh.append(purity_s2d[i])
    efficiency_s2d_optPh.append(efficiency_s2d[i])
    purXeff_s2d_optPh.append(purXeff_s2d[i])


fig = plt.figure(0)
plt.plot(sCuts_2d_optP, purXeff_2d_optP, 'k-', label='selection purity*efficiency with prong purity cut')
plt.plot(sCuts_1d, purXeff_1d, 'k--', label='selection purity*efficiency without prong purity cut')
plt.plot(sCuts_2d_optP, purity_2d_optP, 'r-', label='selection purity with prong purity cut')
plt.plot(sCuts_1d, purity_1d, 'r--', label='selection purity without prong purity cut')
plt.plot(sCuts_2d_optP, efficiency_2d_optP, 'b-', label='selection efficiency with prong purity cut')
plt.plot(sCuts_1d, efficiency_1d, 'b--', label='selection efficiency without prong purity cut')
plt.title("selection results with prong purity > 0.77 and electron confidence > X cuts", fontsize=10)
plt.xlabel("X")
plt.ylim(0.3, 1.25)
#plt.yticks(np.arange(0.35, 1.01, 0.05))
yticks = np.arange(0.3, 1.01, 0.05)
ytick_labels = ['' if i % 2 != 0 else '%.2f'%tick for i, tick in enumerate(yticks)]
plt.yticks(yticks, ytick_labels)
plt.legend(loc='upper left', fontsize=8)
plt.grid(True)

fig = plt.figure(1)
plt.plot(sCuts_1d, purXeff_1d, 'b-', label='selection purity*efficiency')
plt.plot(sCuts_1d, purity_1d, 'g-', label='selection purity')
plt.plot(sCuts_1d, efficiency_1d, 'k-', label='selection efficiency')
plt.title("selection results with electron confidence > X cut", fontsize=12)
plt.xlabel("X")
plt.ylim(0.3, 1.)
#plt.yticks(np.arange(0.35, 1.01, 0.05))
yticks = np.arange(0.3, 1.01, 0.05)
ytick_labels = ['' if i % 2 != 0 else '%.2f'%tick for i, tick in enumerate(yticks)]
plt.yticks(yticks, ytick_labels)
plt.legend(loc='upper left', fontsize=8)
plt.grid(True)

fig = plt.figure(2)
plt.plot(pCuts_2d_optS, purXeff_2d_optS, 'k-', label='selection purity*efficiency')
plt.plot(pCuts_2d_optS, purity_2d_optS, 'r-', label='selection purity')
plt.plot(pCuts_2d_optS, efficiency_2d_optS, 'b-', label='selection efficiency')
plt.title("selection results with electron confidence > 7.1 and prong purity > X cuts", fontsize=10)
plt.xlabel("X")
plt.ylim(0.4, 1.)
#plt.yticks(np.arange(0.35, 1.01, 0.05))
yticks = np.arange(0.4, 1.01, 0.05)
ytick_labels = ['' if i % 2 != 0 else '%.2f'%tick for i, tick in enumerate(yticks)]
plt.yticks(yticks, ytick_labels)
#plt.legend(loc='upper left', fontsize=8)
plt.legend(fontsize=8)
plt.grid(True)

fig = plt.figure(3)
plt.plot(phCuts_s2d_optPi, purXeff_s2d_optPi, 'k-', label='selection purity*efficiency')
plt.plot(phCuts_s2d_optPi, purity_s2d_optPi, 'r-', label='selection purity')
plt.plot(phCuts_s2d_optPi, efficiency_s2d_optPi, 'b-', label='selection efficiency')
plt.title("selection results with pion difference > 8.9 and photon difference > X cuts", fontsize=10)
plt.xlabel("X")
plt.ylim(0.4, 1.)
plt.yticks(np.arange(0.4, 1.01, 0.05))
ytick_labels = ['' if i % 2 != 0 else '%.2f'%tick for i, tick in enumerate(yticks)]
plt.yticks(yticks, ytick_labels)
plt.legend(fontsize=8)
plt.grid(True)

fig = plt.figure(4)
plt.plot(piCuts_s2d_optPh, purXeff_s2d_optPh, 'k-', label='selection purity*efficiency')
plt.plot(piCuts_s2d_optPh, purity_s2d_optPh, 'r-', label='selection purity')
plt.plot(piCuts_s2d_optPh, efficiency_s2d_optPh, 'b-', label='selection efficiency')
plt.title("selection results with photon difference > 3.1 and pion difference > X cuts", fontsize=10)
plt.xlabel("X")
plt.ylim(0.35, 1.)
plt.yticks(np.arange(0.35, 1.01, 0.05))
ytick_labels = ['' if i % 2 != 0 else '%.2f'%tick for i, tick in enumerate(yticks)]
plt.yticks(yticks, ytick_labels)
plt.legend(fontsize=8)
plt.grid(True)

plt.show()
