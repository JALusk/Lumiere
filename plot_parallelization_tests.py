import matplotlib.pyplot as plt
import numpy as np

P = [1, 2, 4, 6, 8]

#Data from Lusk's computer
control_l = [3.213221073, None, None, None, None]
mpi_l = [None, 5.16646719, 1.873862028, 1.57524085, 1.479753494]
thread_l = [None, 1290.144, 2280.8804, 4288.0869, 7304.2]

#Plot Lusk's data
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(P, control_l, label='Linear')
ax1.scatter(P, mpi_l, label='MPI')
ax1.scatter(P, thread_l, label='Threading')
plt.legend()
plt.xlabel('# processes or threads, P')
plt.ylabel('time (s)')
plt.title('MPI vs. Threading: Computer 1')
plt.savefig('../Desktop/mpi_vs_threading_lusk.png')

#Data from Erik's computer
control_e = [1.929160357, None, None, None, None]
mpi_e = [None, 4.474668264, 4.762196779, 5.299091339, None]
thread_e = [None, 856.5893123, 1235.359966, 1725.284964, 3094.525664]

#Plot Erik's data
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(P, control_e, label='Linear')
ax1.scatter(P, mpi_e, label='MPI')
ax1.scatter(P, thread_e, label='Threading')
plt.legend()
plt.xlabel('# processes or threads, P')
plt.ylabel('time (s)')
plt.title('MPI vs. Threading: Computer 2')
plt.savefig('../Desktop/mpi_vs_threading_erik.png')

#Data from Hypatia's computer
control_h = [2.868185043, None, None, None, None]
mpi_h = [5.697077036, 2.646511555, 4.498752832]
thread_h = [None, 3234.731062, 748.9400649, 930.9155359, 1847.476208]

#Plot Hypatia's data
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.scatter(P, control_h, label='Linear')
ax1.scatter([2,6,8], mpi_h, label='MPI')
ax1.scatter(P, thread_h, label='Threading')
plt.legend()
plt.xlabel('# processes or threads, P')
plt.ylabel('time (s)')
plt.title('MPI vs. Threading: Computer 3')
plt.savefig('../Desktop/mpi_vs_threading_hypatia.png')

#Plot all MPI results together
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(P, mpi_l, '-o', label='Computer 1')
ax1.plot(P, mpi_e, '-o', label='Computer 2')
ax1.plot([2,6,8], mpi_h, '-o', label='Computer 3')

plt.legend()
plt.title('Message Passing Interface (MPI)')
plt.xlabel('# processes, P')
plt.ylabel('time (s)')
plt.savefig('../Desktop/mpi_all.png')


#Plot all threading results together
fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.plot(P, thread_l, '-o', label='Computer 1')
ax1.plot(P, thread_e, '-o', label='Computer 2')
ax1.plot(P, thread_h, '-o', label='Computer 3')

plt.legend()
plt.title('Threading')
plt.xlabel('# processes, P')
plt.ylabel('time (s)')
plt.savefig('../Desktop/threading_all.png')

