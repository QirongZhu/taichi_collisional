import numpy as np
import matplotlib.pyplot as plt
import matplotlib 
import matplotlib as mpl

matplotlib.use('pdf')

#print(mpl.rcParams['agg.path.chunksize'])
mpl.rcParams['agg.path.chunksize'] = 10000

import h5py as h5py

fileName = "snapshot_0_ref.hdf5"

file = h5py.File(fileName, "r")

accx0 = np.array(file['Accx'])
accy0 = np.array(file['Accy'])
accz0 = np.array(file['Accz'])
acc   = np.sqrt(accx0*accx0 + accy0*accy0 + accz0*accz0)

jerkx0 = np.array(file['Jerkx'])
jerky0 = np.array(file['Jerky'])
jerkz0 = np.array(file['Jerkz'])
jerk   = np.sqrt(jerkx0*jerkx0 + jerky0*jerky0 + jerkz0*jerkz0)

ax = plt.subplot(111)

fileName = "snapshot_0.hdf5"
file = h5py.File(fileName, "r")
daccx = np.array(file['Accx']) - accx0
daccy = np.array(file['Accy']) - accy0
daccz = np.array(file['Accz']) - accz0
djerkx= np.array(file['Jerkx']) - jerkx0
djerky= np.array(file['Jerky']) - jerky0
djerkz= np.array(file['Jerkz']) - jerkz0
error_acc = np.sqrt(daccx * daccx + daccy * daccy + daccz * daccz)/acc
error_jerk= np.sqrt(djerkx * djerkx + djerky * djerky + djerkz * djerkz)/jerk
sorted_array = np.sort(error_acc)
ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, '--', lw=4, label="test")
sorted_array = np.sort(error_jerk)
#ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, '-', lw=4, label="test")

fileName = "snapshot_0_ref_single_pre.hdf5"
file = h5py.File(fileName, "r")
daccx = np.array(file['Accx']) - accx0
daccy = np.array(file['Accy']) - accy0
daccz = np.array(file['Accz']) - accz0
djerkx= np.array(file['Jerkx']) - jerkx0
djerky= np.array(file['Jerky']) - jerky0
djerkz= np.array(file['Jerkz']) - jerkz0
error_acc = np.sqrt(daccx * daccx + daccy * daccy + daccz * daccz)/acc
error_jerk= np.sqrt(djerkx * djerkx + djerky * djerky + djerkz * djerkz)/jerk
sorted_array = np.sort(error_acc)
ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, 'k--', lw=4, label="test")



ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

ax.set_xlim(1e-12, 1e-0)
#ax.text(2e-4, 2e-5, r'$\epsilon=1e-7$', fontsize='x-large')
plt.xlabel(r"$|\Delta \mathbf{a}|/|\mathbf{a}|$")
plt.ylabel(r"$P\ (>|\Delta \mathbf{a}|/|\mathbf{a}|)$")
plt.title("accelerations")
#plt.legend()
plt.tight_layout()
plt.grid(True)
plt.savefig("acc_errors.pdf")

#print accx
#print jerkx


plt.close()

fileName = "snapshot_0_ref.hdf5"

file = h5py.File(fileName, "r")

accx0 = np.array(file['Accx'])
accy0 = np.array(file['Accy'])
accz0 = np.array(file['Accz'])
acc   = np.sqrt(accx0*accx0 + accy0*accy0 + accz0*accz0)

jerkx0 = np.array(file['Jerkx'])
jerky0 = np.array(file['Jerky'])
jerkz0 = np.array(file['Jerkz'])
jerk   = np.sqrt(jerkx0*jerkx0 + jerky0*jerky0 + jerkz0*jerkz0)

ax = plt.subplot(111)

fileName = "snapshot_0.hdf5"
file = h5py.File(fileName, "r")

daccx = np.array(file['Accx']) - accx0
daccy = np.array(file['Accy']) - accy0
daccz = np.array(file['Accz']) - accz0
djerkx= np.array(file['Jerkx']) - jerkx0
djerky= np.array(file['Jerky']) - jerky0
djerkz= np.array(file['Jerkz']) - jerkz0

error_acc = np.sqrt(daccx * daccx + daccy * daccy + daccz * daccz)/acc
error_jerk= np.sqrt(djerkx * djerkx + djerky * djerky + djerkz * djerkz)/jerk

sorted_array = np.sort(error_acc)
#ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, '-', lw=4, label="test")
sorted_array = np.sort(error_jerk)
ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, '--', lw=4, label="test")



fileName = "snapshot_0_ref_single_pre.hdf5"
file = h5py.File(fileName, "r")
daccx = np.array(file['Accx']) - accx0
daccy = np.array(file['Accy']) - accy0
daccz = np.array(file['Accz']) - accz0
djerkx= np.array(file['Jerkx']) - jerkx0
djerky= np.array(file['Jerky']) - jerky0
djerkz= np.array(file['Jerkz']) - jerkz0
error_acc = np.sqrt(daccx * daccx + daccy * daccy + daccz * daccz)/acc
error_jerk= np.sqrt(djerkx * djerkx + djerky * djerky + djerkz * djerkz)/jerk
sorted_array = np.sort(error_jerk)
ax.plot(sorted_array, 1-np.arange(0, 1e5)/1e5, 'k--', lw=4, label="test")

ax.set_xscale("log", nonposx='clip')
ax.set_yscale("log", nonposy='clip')

ax.set_xlim(1e-12, 1e-0)
#ax.text(2e-4, 2e-5, r'$\epsilon=1e-7$', fontsize='x-large')
plt.xlabel(r"$|\Delta \mathbf{j}|/|\mathbf{j}}|$")
plt.ylabel(r"$P\ (>|\Delta \mathbf{j}|/|\mathbf{j}|)$")
plt.title("jerks")
#plt.legend()
plt.tight_layout()
plt.grid(True)
plt.savefig("jerk_errors.pdf")
