
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt

pcaxis = [0.121072, -0.120848, 0.985260]
pccenter = [73.676590, -62.803082, 637.737976]
vertex = [65.89543151855469, -54.615806579589844, 573.2275390625]
shw1 = [65.951630, -54.576836, 573.931580]
shw2 = [66.622658, -56.417892, 572.640503]

xs = []
ys = []
zs = []

with open("PCAtest_points.txt") as points_list:
  for line in points_list:
    points = line.replace("\n","").split()
    xs.append(float(points[0]))
    ys.append(float(points[1]))
    zs.append(float(points[2]))

xs_proj = []
ys_proj = []
zs_proj = []

for i in range(len(xs)):
  p = [xs[i], ys[i], zs[i]]
  lp = [0.,0.,0.]
  ldist = 0.
  for c in range(3):
    ldist += (p[c] - pccenter[c])*pcaxis[c]
  for c in range(3):
    lp[c] = pccenter[c] + ldist*pcaxis[c]
  xs_proj.append(lp[0])
  ys_proj.append(lp[1])
  zs_proj.append(lp[2])

specialPs = [vertex, shw1, shw2]
specialPs_proj = []

for p in specialPs:  
  lp = [0.,0.,0.]
  ldist = 0.
  for c in range(3):
    ldist += (p[c] - pccenter[c])*pcaxis[c]
  for c in range(3):
    lp[c] = pccenter[c] + ldist*pcaxis[c]
  specialPs_proj.append(lp)

vtxProj = specialPs_proj[0]
oneDimProj = []
for i in range(len(xs)):
  oneDimProj.append(sqrt((xs_proj[i] - vtxProj[0])**2 + (ys_proj[i] - vtxProj[1])**2 + (zs_proj[i] - vtxProj[2])**2))

oneDimProj_sorted = sorted(oneDimProj)
maxGap = -1.
for i in range(1, len(oneDimProj_sorted)):
  gap = oneDimProj_sorted[i] - oneDimProj_sorted[i-1]
  if gap > maxGap:
    maxGap = gap

print("largest gap:", maxGap)

oneDimProj = np.array(oneDimProj)
oneDimProjY = np.zeros(len(oneDimProj))


fig = plt.figure(0)
ax = fig.add_subplot(projection='3d')
#ax.scatter(xs, ys, zs)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
ax.scatter(zs, xs, ys, marker='.', s=1)
ax.scatter([specialPs[0][2]], [specialPs[0][0]], [specialPs[0][1]], marker='P', s=80)
ax.scatter([specialPs[1][2]], [specialPs[1][0]], [specialPs[1][1]], marker='s', s=20, c='r')
ax.scatter([specialPs[2][2]], [specialPs[2][0]], [specialPs[2][1]], marker='s', s=20, c='r')
ax.set_xlabel('z')
ax.set_ylabel('x')
ax.set_zlabel('y')

fig = plt.figure(1)
ax = fig.add_subplot(projection='3d')
#ax.scatter(xs_proj, ys_proj, zs_proj)
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
ax.scatter(zs_proj, xs_proj, ys_proj)
ax.scatter([specialPs_proj[0][2]], [specialPs_proj[0][0]], [specialPs_proj[0][1]], marker='P', s=80)
ax.scatter([specialPs_proj[1][2]], [specialPs_proj[1][0]], [specialPs_proj[1][1]], marker='s', s=20, c='r')
ax.scatter([specialPs_proj[2][2]], [specialPs_proj[2][0]], [specialPs_proj[2][1]], marker='s', s=20, c='r')
ax.set_xlabel('z')
ax.set_ylabel('x')
ax.set_zlabel('y')

fig = plt.figure(2)
plt.plot(oneDimProj, oneDimProjY, '.')

plt.show()



