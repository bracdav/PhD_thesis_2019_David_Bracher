from thesis_db2019 import *

path = 'D:\\Beamtime_SLS\\Beamtime 1805_RHEED\\David\\180508_goe_etoh_rheed_002\\analysis\\data\\'
files = glob(path + '*etoh*.dat')

path_sim = 'D:\\Beamtime_SLS\\Beamtime 1805_RHEED\\David\\180508_goe_etoh_rheed_002\\analysis\\data\\'
files_sim = glob(path_sim + '*mo*.dat')
sim_names = ['Magnetite', 'Goethite', 'Hematite']


rheed = RheedData(files[0])
data = rheed.load_data(delimiter= '\t')
rheed.set_zero(9)
rheed.set_temp()
rheed.px_to_2theta()
rheed.set_energy(30_000)
rheed.get_wavelength()
rheed.set_q()


plt.figure()
plt.plot(rheed.data.q, rheed.data.Y)
plt.xlim([3.5, 4.5])
plt.show()


sim = RheedSim(files_sim, sim_names)

