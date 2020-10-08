from thesis_db2019 import *

path = 'D:\\Beamtime_SLS\\Beamtime 1706\\David\\Analysis\\feooh'
en = path + '\\s170612_005.dat'
stem = 's170612'
suffix = '_raw_small.dat'

pol1 = XASpectrum()
pol1.construct_nums([5, 7, 8, 9, 10])
pol1.load_energy(en)
pol1.xas_load(path = path, stem = stem, suffix = suffix)
pol1.xas_normalize()
pol1.xas_mean()

pol2 = XASpectrum()
pol2.construct_nums([11, 12, 13, 14, 15])
pol2.load_energy(en)
pol2.xas_load(path = path, stem = stem, suffix = suffix)
pol2.xas_normalize()
pol2.xas_mean()

Dichroism = pd.DataFrame()
dichr = pd.DataFrame()
l3 = pd.DataFrame()
l32 = pd.DataFrame()
l2 = pd.DataFrame()
for np in pol1.avg:
    dichroism = FeOOH(pol1.avg[np],pol2.avg[np])
    dichroism.select_feature(inter = [0,24])
    dichroism.norm_on_feature()
    dichroism.xas_match()
    
    Dichroism[np+'_pol1'] = dichroism.match1
    Dichroism[np+'_pol2'] = dichroism.match2
    dichroism.calc_dichroism()
    Dichroism[np+'_dichroism'] =dichroism.dichroism
    dichr[np+'_dichroism'] = dichroism.dichroism
    dichroism.calc_dichroism_norm()
    Dichroism[np+'_dichroism_norm'] = dichroism.dichroism_norm
    
    try:
        dichroism.det_features(dichr, [70, 90])
        dichroism.fit_features_l32(710.8, [0.6, 1.6, 2.4, 13.6, 15.2])
        l32[np] = dichroism.fit_l32
        dichroism.fit_features_l2(724.4, [0.6])
        l2[np] = dichroism.fit_l2
        dichroism.fit_features_l3(710.8, [0.6, 1.6, 2.4])
        l3[np] = dichroism.fit_l3
        print('fit successfull for {}'.format(np))
    except:
        print('no fit found for {}'.format(np))
        
l32.index = l2.index = l3.index = Dichroism.index