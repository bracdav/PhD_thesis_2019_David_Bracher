from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import scipy.constants as nat
import seaborn as sns
from Stoner.Fit import langevin as lv


"""Analysis Functions as used in thesis"""
def ScanPhi90(X,I01,I2,THETA,PHI):
    """
    function to model f(x) = I0 + I2 * cos(x)**2
    Parameters
    ----------
    X : (n,) array of floats
        array of floats used to calculate the function
    I01 : float
        Offset parameter
    I2 : float
        Amplitude parameter
    THETA : float
        Out of plane orientation of the  spin axis
    PHI : float
        In plane orientation of the spin axis

    Returns
    -------
    TYPE
        (n,) array of floats

    """
   
    X=X*np.pi/180
    THETAS1 = 90 *np.pi / 180
    PHI = PHI * np.pi / 180
    THETA = THETA * np.pi / 180
    return I01 - I2 * ((np.cos(X)*np.cos(THETAS1)*np.sin(16*np.pi/180) + np.sin(X)*np.sin(THETAS1))*np.sin(THETA)*np.cos(PHI) - (-np.sin(X)*np.cos(THETAS1)*np.sin(16*np.pi/180) + np.cos(X)*np.sin(THETAS1))*np.sin(THETA)*np.sin(PHI) + np.cos(THETAS1)*np.cos(16*np.pi/180)*np.cos(THETA))**2  
    
def ScanPhi0(X,I02,I2,THETA,PHI):
    """
    function to model f(x) = I0 + I2 * cos(x)**2
    Parameters
    ----------
    X : (n,) array of floats
        array of floats used to calculate the function
    I02 : float
        Offset parameter
    I2 : float
        Amplitude parameter
    THETA : float
        Out of plane orientation of the  spin axis
    PHI : float
        In plane orientation of the spin axis

    Returns
    -------
    TYPE
        (n,) array of floats

    """
    
    X=X*np.pi/180
    THETAS2 = 0
    PHI = PHI * np.pi / 180
    THETA = THETA * np.pi / 180
    return I02 - I2 * ((np.cos(X)*np.cos(THETAS2)*np.sin(16*np.pi/180) + np.sin(X)*np.sin(THETAS2))*np.sin(THETA)*np.cos(PHI) - (-np.sin(X)*np.cos(THETAS2)*np.sin(16*np.pi/180) + np.cos(X)*np.sin(THETAS2))*np.sin(THETA)*np.sin(PHI) + np.cos(THETAS2)*np.cos(16*np.pi/180)*np.cos(THETA))**2

def ScanTheta(X,I03,I2,THETA,PHI):
    """
    function to model f(x) = I0 + I2 * cos(x)**2
    Parameters
    ----------
    X : (n,) array of floats
        array of floats used to calculate the function
    I03 : float
        Offset parameter
    I2 : float
        Amplitude parameter
    THETA : float
        Out of plane orientation of the  spin axis
    PHI : float
        In plane orientation of the spin axis

    Returns
    -------
    float
        (n,) array of floats

    """

    X=X*np.pi/180
    PHI = PHI * np.pi / 180
    THETA = THETA * np.pi / 180
    PHIS = 20 * np.pi / 180
    return I03 - I2 * ((np.cos(PHIS)*np.cos(X)*np.sin(16*np.pi/180) + np.sin(PHIS)*np.sin(X))*np.sin(THETA)*np.cos(PHI) - (-np.sin(20*np.pi/180)*np.cos(X)*np.sin(16*np.pi/180) + np.cos(PHIS)*np.sin(X))*np.sin(THETA)*np.sin(PHI) + np.cos(X)*np.cos(16*np.pi/180)*np.cos(THETA))**2

def Find_S(ANGLES,I01,I02,I03,I2,THETA,PHI):
    """
    Parameters
    ----------
    ANGLES : (n,) array of floats
        array of floats used to calculate the function.
    I01 : float
        Offset parameter
    I02 : float
        Offset parameter.
    I03 : float
        Offset parameter.
    I2 : float
    Amplitude parameter
        THETA : float
        Out of plane orientation of the  spin axis.
    PHI : float
        In plane orientation of the spin axis.

    Returns
    -------
    out : float
        (n,) array of floats.

    """

    ext1 = ANGLES[:6]
    ext2 = ANGLES[6:13]
    ext3 = ANGLES[13:]
    
    res1 = ScanPhi90(ext1,I01,I2,THETA,PHI)
    res2 = ScanPhi0(ext2,I02,I2,THETA,PHI)
    res3 = ScanTheta(ext3,I03,I2,THETA,PHI)
    out = np.append(res1,res2)
    out = np.append(out,res3)
    return out



"""Analysis classes as used in thesis"""
class XASpectrum:
    def __init__(self):
        pass
    
    def construct_nums(self,nums = None):
        """
        Parameters
        ----------
        nums : list of integers
            DESCRIPTION. The default is None.

        Returns
        -------
        List of strings
            The input numbers are converted to strings of numbers between
            '0001' - '999' in self.string_nums
        """
        strings = [] 
        for num in nums:
            if num < 10:
                strings.append('00'+str(num))
            elif num >= 100:
                strings.append(str(num))
            else:
                strings.append('0'+str(num))
        self.string_nums = strings
    
    def load_energy(self, location, sep = '\t'):
        """
        loads the energy scale from a file
        Parameters
        ----------
        location : String
            should contain the location of the energy file (output spectrum from PEEM)
        sep : String, optional
            The default is '\t'.
            sep is used to determine the seperator in the text file
        
        Returns
        -------
        No return

        """
        energy = pd.read_csv(location, delimiter = sep)
        self.energy = energy['Energy'].values        
    
    def xas_load(self,path = None, stem = None, suffix = '.csv', sep = '\t'):
        """
        loads raw XA spectra
        Parameters
        ----------
        path : str, optional
            Location where spectroscopy files are saved
        stem : str, optional
            stem of files
        suffix : str, optional
            The default is '.csv'.
            Suffix of the text files
        sep : str, optional
            The default is '\t'.
            sep defines the seperator in the text files

        Returns
        -------
        No return
        """
        self.path = path
        self.stem = stem
        self.file_loc = file_loc = path + '\\' + stem + '_'
        try:
            df = pd.DataFrame()
            for num in self.string_nums:
                print(file_loc + num + suffix)
                single = pd.read_csv(file_loc + num + suffix, delimiter = sep)
                self.single = single
                single['supl'] = num
                single.index = np.round(self.energy,1)
                df = pd.concat((df, single))
            first = df.columns[0]
            df.drop(first, axis = 1, inplace = True)
            self.raw_data = df
        except:
            print('Method .construct_nums() and .load_energy() need to be run first.')
            
    def xas_normalize(self):
        """
        normalizes PEEM XA spectra.
        Returns
        -------
        No return

        """
        try:
            data = self.raw_data.values[:,:-1]
            norm = data[:,::2] / data[:,1::2]
            cols = []
            self.norm = pd.DataFrame(norm).astype('float')
            self.norm.index = self.raw_data.index
            for i in range(1,self.norm.shape[1]+1):
                cols.append('np_'+str(i))
            self.norm.columns = cols
            #self.norm['supl'] = self.raw_data['supl']
            
        except:
            print('Error : Run .xas_load() properly first')
        return self
            
    def xas_mean(self):
        """
        Calculate the average of a number of spectra

        Returns
        -------
        None.

        """
        self.avg = self.norm.reset_index()
        self.avg = self.avg.groupby('index').mean()    
    
class DichroicSpectrum:
    
    def __init__(self,pol1,pol2):
        self.pol1 = pol1
        self.pol2 = pol2
    
    def select_feature(self, inter = None, plot = 0):
        """
        selets a feature of a spectrum to normalize.
        Parameters
        ----------
        inter : TYPE, optional
            Intervall in which the features are located.
        plot : TYPE, optional
            1 will plot the selected intervall and 0 supresses plotting.
            0 is the standard value

        Returns
        -------
        No return

        """
        if inter == None:
            inter = [0, len(self.pol1.index)]
        if plot == 1:
            plt.plot(self.pol1.values)
            plt.xlim(inter)
            plt.show()
        
        self.inter = inter
    
    def norm_on_feature(self):
        """
        calculates the normalized spectrum
        Returns
        -------
        No return
        """
        start = self.inter[0]
        end = self.inter[1]
        nom = self.pol1 - self.pol1.iloc[:10].mean()
        denom = self.pol1.iloc[start:end].max() - self.pol1.iloc[:10].mean()
        self.pol1 = nom/denom
        nom = self.pol2 - self.pol2.iloc[:10].mean()
        denom = self.pol2.iloc[start:end].max() - self.pol2.iloc[:10].mean()
        self.pol2 = nom/denom                                                 
    
    def xas_match(self, intervall = [], int_given = False):
        """
        match the spectrra of two polarization using a least-square fit

        Parameters
        ----------
        intervall : list of integers, optional
            determines the points to consider for the fit
        int_given : Bool, optional
            Boolean (True/False) which determines if a list of integers is given.
            The default is False, meaning no list is entered.
        """
        
        if int_given == False:
            
            def match(x,a,b,c):
                return self.pol1[intervall] * a + x * b + c
            fit, cov = curve_fit(match, self.pol1.index, self.pol2.values)
        else:
            def match(x,a,b,c):
                return self.pol1[intervall] * a + x * b + c
            fit, cov = curve_fit(match,self.pol1.index[intervall], 
                                 self.pol2.values[intervall])
        
        self.match1 = self.pol1 * fit[0] + self.pol1.index * fit[1] +fit[2]
        self.match2 = self.pol2
    
    def calc_dichroism(self):
        """
        calculate the dichroism spectrum as pol1 - pol2
        Returns
        -------
        None.
        """
        self.dichroism = self.match1 - self.match2
    
    def calc_dichroism_norm(self):
        """
        calculate the dichroism spectrum as (pol1 - pol2)/(pol1 + pol2)
        Returns
        -------
        TYPE
            DESCRIPTION.
        """
        self.dichroism_norm = (self.match1 - self.match2)/(2+self.match1 - self.match2)

class FeOOH(DichroicSpectrum):
    
    def det_features(self, dichr, lim):
        """
        Parameters
        ----------
        dichr : TYPE
            XA spectrum
        lim : list of integer
            Image numbers which should be displayed.
        Returns
        -------
        None.
        """
        plt.plot(dichr.mean(axis = 1).values)
        plt.xlim(lim)
        
    def fit_features_l3(self, cent, shift):
        """
        calculate Gaussian fit to spectral signature L3 edge
        Parameters
        ----------
        cent : int
            The location where the first Gaussian is centered.
        shift : list of integers
            List of shifts to the remaining Gaussian features.
        Returns
        -------
        (n,) array
            The fitted line

        """
        def fun(x, ampl, cent, width):
            return ampl * np.exp(-(x - cent)**2/(2 * width**2))
        
        def fit_full(x, a1, a2, a3, a4, w1, w2, w3, w4):
            return fun(x, a1, cent, w1) + fun(x, a2, cent+shift[0], w2) + fun(x, a3, cent+shift[1], w3) + fun(x,a4, cent+shift[2], w4)
        
        fit, cov = curve_fit(fit_full, self.dichroism.index, self.dichroism.values)
        fit_line = fit_full(self.dichroism.index, fit[0], fit[1], fit[2], 
                            fit[3], fit[4], fit[5], fit[6], fit[7])
        self.fit_l3 = fit_line
        
    def fit_features_l32(self, cent, shift):
        """
        calculate Gaussian fit to spectral signature L2/3 edges
        Parameters
        ----------
        cent : int
            The location where the first Gaussian is centered.
        shift : list of integers
            List of shifts to the remaining Gaussian features.
        Returns
        -------
        (n,) array
            The fitted line

        """
        def fun(x, ampl, cent, width):
            return ampl * np.exp(-(x - cent)**2/(2 * width**2))
        
        def fit_full(x, a1, a2, a3, a4, a5, a6, w1, w2, w3, w4, w5, w6):
            return fun(x, a1, cent, w1) + fun(x, a2, cent+shift[0], w2) + fun(x, a3, cent+shift[1], w3) + fun(x,a4, cent+shift[2], w4) + fun(x,a5, cent+shift[3], w5) + fun(x,a6, cent+shift[4], w6)
        
        fit, cov = curve_fit(fit_full, self.dichroism.index, self.dichroism.values)
        
        fit_line = fit_full(self.dichroism.index, fit[0], fit[1], fit[2], 
                            fit[3], fit[4], fit[5], fit[6], fit[7], fit[8], 
                            fit[9], fit[10], fit[11])
        self.fit_l32 = fit_line
        return self
        
    def fit_features_l2(self, cent, shift):
        """
        calculate Gaussian fit to spectral signature L2 edge
        Parameters
        ----------
        cent : int
            The location where the first Gaussian is centered.
        shift : list of integers
            List of shifts to the remaining Gaussian features.
        Returns
        -------
        (n,) array
            The fitted line

        """
        def fun(x, ampl, cent, width):
            return ampl * np.exp(-(x - cent)**2/(2 * width**2))
        
        def fit_full(x, a1, a2, w1, w2):
            return fun(x, a1, cent, w1) + fun(x, a2, cent+shift[0], w2)
        
        fit, cov = curve_fit(fit_full, self.dichroism.index,
                             self.dichroism.values)
        fit_line = fit_full(self.dichroism.index, fit[0], fit[1], fit[2], 
                            fit[3])
        self.fit_l2 = fit_line

class CoO(DichroicSpectrum):
    """
    Empty class inherits from DichroicSpectrum. Specific methods for CoO analysis
    should be appended in this class.
    
    (Please remove disclaimer as soon as a method is introduced.)
    """
    pass


class SpinAxis:
    
    def __init__(self, files):
        """
        Parameters
        ----------
        files : List
            List of files.

        Returns
        -------
        None.
        """
        self.files = files
        
    def load_files(self, names):
        """
        Load the full angular dependency data of nanoparticles.
        Parameters
        ----------
        names : List of names
            Names for files.

        Returns
        -------
        (n,m) arrays
            angular dependent XMLD data of nanoparticles
        """
        self.names = names
        zipped = zip(self.files, names)
        self.file_dict = {}
        self.data = {}
        for file, name in zipped:
            self.file_dict[name] = file
            self.data[name] = pd.read_csv(file, header = None)
            self.data[name].set_index(self.data[name].columns[0],
                                      inplace = True)
            particles = int(len(self.data[name].columns)/2)
            cols = []
            self.np_names = []
            for np in range(particles):
                cols.append('np_'+str(np+1))
                cols.append('np_'+str(np+1)+'err')
                self.np_names.append('np_'+str(np+1))
            self.data[name].columns = cols
        return self.data
            
    def construct_fit(self, np):
        """
        Parameters
        ----------
        np : String
            Name of the nanoparticle which should be analyzed.

        Returns
        -------
        None.
        """
        self.combined = pd.DataFrame()
        for name in self.names:
            self.combined = pd.concat((self.combined, self.data[name][np]))
        self.scale = self.combined.index
        self.meas = self.combined.values.transpose()[0]
        
    def execute_fit(self):
        """
        perform the fits to the XMLD data
        Returns
        -------
        self.fit: List
            Fit parameters of global fit.
            
        self.cov: ndimArray
            Covariance matrix of the global fit.
        """
        init = np.array([0.7, 0.7, 0.7, 0.5, 75, 225])
        bounds = ((0, 0, 0, 0, 0, 0), (3, 3, 3, 2, 180, 360))
        self.fit, self.cov = curve_fit(Find_S,self.scale,self.meas,
                                       p0 = init, bounds = bounds)
        return self.fit, self.cov
    
    def sim_lines(self, intervall =[0,100], points = 1000, dims = 4):
        """
        calculate the simulated fit lines of the XMLD fits according to the fi parameters
        Parameters
        ----------
        intervall : List, optional
            The default is [0,100]. Intervall which is simulated.
        points : TYPE, optional
            The default is 1000. Number of point for simulated curves.
        dims : TYPE, optional
            The default is 4. Amount of lines to simulate.

        Returns
        -------
        self.sim pandas dataframe
            All the simulated lines.

        """
        self.sim = {'scale' : np.linspace(intervall[0], intervall[1], points)}
        params = self.fit
        self.sim[self.names[0]] = ScanPhi0(self.sim['scale'],params[0],
                                           params[3], params[4], params[5])
        self.sim[self.names[1]] = ScanPhi90(self.sim['scale'],params[1],
                                            params[3], params[4], params[5])
        self.sim[self.names[2]] = ScanTheta(self.sim['scale'],params[2],
                                            params[3], params[4], params[5])
        self.sim = pd.DataFrame(self.sim)
        self.sim.set_index('scale', inplace = True)
        return self.sim
    
    
    
class XASExtreme:
    
    def __init__(self, path):
        """
        Parameters
        ----------
        path : str
            String of the path where the XA spectra files are located

        Returns
        -------
        None.

        """
        self.path = path
        
    def list_files(self, pattern):
        """
        list files in path
        Parameters
        ----------
        pattern : str
            pattern to match in order to list all files within a locataction

        Returns
        -------
        None.

        """
        self.files = glob(self.path+'\\'+pattern)
    
    def seperate_polarization_circ(self, delimiter):
        """
        seperate the names of files of recorded with c+ and c-

        Parameters
        ----------
        delimiter : str
            delimiter which is used in the spectrum file ('\t', ',')

        Returns
        -------
        None.

        """
        self.delimiter = delimiter
        self.pol1_list = []
        self.pol2_list = []
        for file in self.files:
            data = pd.read_csv(file, delimiter = delimiter, usecols = ['FieldX', 'Pol'],
                               nrows = 2).drop(0)
            data['FieldX'] = data['FieldX'].astype('float').round().astype('int')
            data = data.values
            data = data[0]
            if (float(data[0]) > 0) & (data[1] == 'CIRC +'):
                self.pol1_list.append(file)
            if (float(data[0]) < 0) & (data[1] == 'CIRC -'):
                self.pol1_list.append(file)
            if (float(data[0]) < 0) & (data[1] == 'CIRC +'):
                self.pol2_list.append(file)
            if (float(data[0]) > 0) & (data[1] == 'CIRC -'):
                self.pol2_list.append(file)
            if (float(data[0]) == 0) & (data[1] == 'CIRC +'):
                self.pol1_list.append(file)
            if (float(data[0]) == 0) & (data[1] == 'CIRC -'):
                self.pol2_list.append(file)
                
                
    def load_data(self,cols = ['#Ecrbk', 'FieldX', 'FieldZ', 'Pol', 'Temperature', 'NORMtey']):
        """
        load all data within a certain path
        Parameters
        ----------
        cols : list of str
            all columns which are used from file. standard is:
                ['#Ecrbk', 'FieldX', 'FieldZ', 'Pol', 'Temperature', 'NORMtey'] 

        Returns
        -------
        None

        """
        def process_df(df, mapper):
            """
            preprocessing data within the load data method.
            Parameters
            ----------
            df : pandas DataFrame
                raw data
            mapper : dictionary
                dictionary is used

            Returns
            -------
            df : pandas DataFrame
                processed data.

            """
            df.drop(0, inplace = True)
            df['#Ecrbk'] = df['#Ecrbk'].astype('float').round(1)
            cols = ['FieldX', 'FieldZ', 'Temperature']
            df[cols] = df[cols].astype('float').round(0)
            df['NORMtey'] = df['NORMtey'].astype('float')
            cols = df.columns
            df = df.groupby('#Ecrbk', as_index = False).agg(mapper)
            df.columns = cols
            return df
        self.plus = pd.DataFrame()
        self.minus = pd.DataFrame()
        mapper = {'FieldX' : ['mean'], 'FieldZ' : ['mean'], 'Pol' : ['first'], 
                  'Temperature' : ['first'], 'NORMtey' : ['mean']}
        for p, m in zip(self.pol1_list, self.pol2_list):
             P = pd.read_csv(p, delimiter = self.delimiter, usecols = cols)
             P = process_df(P, mapper)
             self.plus = pd.concat((self.plus, P), axis = 0)
             M = pd.read_csv(m, delimiter = self.delimiter, usecols = cols)
             M = process_df(M, mapper)
             self.minus = pd.concat((self.minus, M),axis = 0)
    
            
    def aggregate_data_circ(self):
        """
        Average data of all measurements.
        Parameters
        ----------
        None

        Returns
        -------
        self.raw_spectrua : DataFrame
            averaged XA spectra

        """
        plus = self.pol1 = self.plus.groupby('#Ecrbk')['NORMtey'].mean()
        minus = self.pol2 = self.minus.groupby('#Ecrbk')['NORMtey'].mean()
        self.raw_spectra = pd.DataFrame()
        self.raw_spectra['c+'] = plus
        self.raw_spectra['c-'] = minus
        self.raw_spectra['xmcd'] = plus - minus
        self.raw_spectra.columns = ['c+', 'c-', 'xmcd']
        return self.raw_spectra
    
    def step_norm(self, interval):
        """
        subtracting linear background.
        Parameters
        ----------
        interval : list of ints
            list of intervals

        Returns
        -------
        None
        """
        interval = np.append(np.arange(interval[0],interval[1],1), np.arange(interval[2], interval[3], 1))
        def lin(x, a, b):
            return a*x + b
        fit, cov = curve_fit(lin, self.pol2.index[interval], self.pol2.values[interval])
        self.pol2 = self.pol2 - lin(self.pol2.index, fit[0], fit[1])
        self.pol2 = self.pol2/self.pol2.max()

    def xas_match(self, intervall = [], int_given = False):
        """
        match XA spectra in order to get dichroic spectrum
        Parameters
        ----------
        intervall : (n,) array of ints, optional
            DESCRIPTION. The default is [].
        int_given : Bool, optional
            The default is False.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        if int_given == False:        
            def match(x,a,b,c):
                return self.pol1 * a + x * b + c                    
            fit, cov = curve_fit(match,self.pol2.index, self.pol2.values)
        else:
            pol1 = self.pol1.iloc[intervall]
            pol2 = self.pol2.iloc[intervall]  
            def match(x,a,b,c):
                return pol1.values * a + x * b + c
            fit, cov = curve_fit(match,pol2.index, pol2.values)    
        self.match1 = self.pol1 * fit[0] + self.pol1.index * fit[1] +fit[2]
        self.match2 = self.pol2
    
    def construct_data(self):
        """
        construct the final data sets
        Returns
        -------
        DataFrame
            Results

        """
        self.final = pd.DataFrame()
        self.final['c+'] = self.match1
        self.final['c-'] = self.match2
        self.final['xmcd'] = self.final['c+'] - self.final['c-']
        self.final['temperature'] = self.plus.groupby('#Ecrbk')['Temperature'].agg('first')
        self.final['fieldx'] = np.abs(self.plus.groupby('#Ecrbk')['FieldX'].agg('first'))
        return self.final
        
    def construct_interval(self,ranges, full_xas = True, sample = 5):
        """
        generate intervals for method xas_match()
        Parameters
        ----------
        ranges : list of integers
            
        full_xas : Bool, optional
            decides if only L3 or L2/3 edge (standard) is used. The default is True.
        sample : TYPE, optional
            DESCRIPTION. The default is 5.

        Returns
        -------
        None.

        """
        if full_xas == True:
            preedge = np.arange(ranges[0][0], ranges[0][1], 1)
            l3 = np.arange(ranges[1][0], ranges[1][1], sample)
            between = np.arange(ranges[2][0], ranges[2][1], 1)
            l2 = np.arange(ranges[3][0], ranges[3][1], sample)
            postedge = np.arange(ranges[4][0], ranges[4][1], 1)
            self.intervall = preedge
            self.intervall = np.append(self.intervall, l3)
            self.intervall = np.append(self.intervall, between)
            self.intervall = np.append(self.intervall, l2)
            self.intervall = np.append(self.intervall, postedge)
            
        else:
            preedge = np.arange(ranges[0][0], ranges[0][1], 1)
            l3 = np.arange(ranges[1][0], ranges[1][1], undersample)
            postedge = np.arange(ranges[2][0], ranges[2][1], 1)
            self.intervall = preedge
            self.intervall = np.append(self.intervall, l3)
            self.intervall = np.append(self.intervall, postedge)
            
            
class SumRules:
    
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            string of destination of file

        Returns
        -------
        None.

        """
        self.filename = filename
        
    def load_data(self):
        """
        load data with spectra
        Returns
        -------
        None.

        """
        self.data = pd.read_csv(self.filename)
        self.data['comb'] = self.data.temperature.astype('str') + ',' + self.data.fieldx.astype('str')
        
    def get_param_pairs(self):
        """
        get the pairs of temperature and magnetic field
        Returns
        -------
        None.

        """
        params = pd.Series(self.data.comb.unique()).str.split(',').values
        self.params = []
        for i, j in params:
            self.params.append([float(i), float(j)])

        
    def use_parm(self, temp, h):
        """
        get data according to parameters temperature and magnetic field
        Parameters
        ----------
        temp : int
            temperature
        h : int
            magnetic field

        Returns
        -------
        None.

        """
        self.P1 = self.data.loc[(self.data.temperature == temp) & (self.data.fieldx == h), 'c+']
        self.P2 = self.data.loc[(self.data.temperature == temp) & (self.data.fieldx == h), 'c-']
                
    def calculate_A(self, interval, off = 0):
        """
        Calculate the A for sum rule analysis.
        Parameters
        ----------
        interval : list of integer
            list conatins [start, end]
        off : float, optional
            offset parameter. The default is 0.

        Returns
        -------
        A : float
            DESCRIPTION.

        """
        A = [0]
        A_int = 0
        EN = self.data.index
        P1 = self.P1.values
        P2 = self.P2.values +off
        
        for i in range(interval[0], interval[1] + 1):
            a_int = (EN[i+1] - EN[i]) * ((P1[i] - P2[i]) + (P1[i+1] - P2[i+1]))/2
            A_int += a_int
            A.append(A_int)
        self.A = A
        self.A_int = A_int
        return A
    
    
    def calculate_B(self, interval, off = 0):
        """
        Calculate the B for sum rule analysis.
        Parameters
        ----------
        interval : list of integer
            list conatins [start, end]
        off : float, optional
            offset parameter. The default is 0.

        Returns
        -------
        B : float
            DESCRIPTION.

        """
        B = [0]
        B_int = 0
        EN = self.data.index
        P1 = self.data['c+']
        P2 = self.data['c-'] + off
        
        for i in range(interval[0], interval[1] + 1):
            b_int = (EN[i+1] - EN[i]) * ((P1[i] - P2[i]) + (P1[i+1] - P2[i+1]))/2
            B_int += b_int
            B.append(B_int)
        self.B = B
        self.B_int = B_int
        return B
    
    
    def calculate_C(self, interval, off = 0):
        """
        Calculate the C for sum rule analysis.
        Parameters
        ----------
        interval : list of integer
            list conatins [start, end]
        off : float, optional
            offset parameter. The default is 0.

        Returns
        -------
        C : TYPE
            DESCRIPTION.

        """
        C = [0]
        C_int = 0
        EN = self.data.index
        P1 = self.data['c+']
        P2 = self.data['c-'] + off
        
        for i in range(interval[0], interval[1] + 1):
            c_int = (EN[i+1] - EN[i]) * ((P1[i] + P2[i]) + (P1[i+1] + P2[i+1]))/4
            C_int += c_int
            C.append(C_int)
        self.C = C
        self.C_int = C_int
        return C
    
    def calculate_orb(self, Nd = 6):
        """
        Calculate the orbital moments.
        Parameters
        ----------
        Nd : int, optional
            number of electron holes of absorber atom. The default is 6.

        Returns
        -------
        float
            DESCRIPTION.

        """
        self.orb = -2*(10-Nd)*(self.A_int + self.B_int)/self.C_int
        return self.orb
    
    def calculate_spin(self, Nd = 6):
        """
        Calculate the spin moments.
        Parameters
        ----------
        Nd : int, optional
            number of electron holes of absorber atom. The default is 6.

        Returns
        -------
        float
            DESCRIPTION.

        """
        self.spin = -(10-Nd)*3 * (self.A_int - 2*self.B_int)/(2 * self.C_int)
        return self.spin
        
        
class MagnetizationCurves:
    
    def __init__(self, filename):
        """
        Parameters
        ----------
        filename : str
            DESCRIPTION.

        Returns
        -------
        None.
        """        
        self.filename = filename
        
    def load_data(self):
        """
        Load the file containing the calculated, spin, orbital, and total moments
        Returns
        -------
        DataFrame
            Contains the spin, orbita, and total moments with respect to temperature and magnetic field

        """
        self.moments = pd.read_csv(self.filename)
        return self.moments
    
    def set_zero(self):
        """
        Moments ar set to 0 ifmagnetic field is 0.
        Returns
        -------
        None.

        """
        self.moments.loc[self.moments.H == 0, ['20K', '20K_err', '50K',
                                               '50K_err', '300K',
                                               '300K_err']] = 0
        
    def plot_total(self, tp = 'total'):
        """
        plot the data of the magnetic moments.
        Parameters
        ----------
        tp : str, optional
            string which describes which moments are plotted (total, orbit, spin). The default is 'total'.

        Returns
        -------
        None.

        """
        red = self.moments.loc[self.moments.type == tp+'_moment']
        red.drop('type', inplace = True, axis = 1)
        plt.figure()
        plt.errorbar(red.H, red['20K'], red['20K_err'], color = 'blue',
                     capsize = 5, marker = 'o', linestyle = '')
        plt.errorbar(red.H, red['50K'], red['50K_err'], color = 'black',
                     capsize = 5, marker = 'o', linestyle = '')
        plt.errorbar(red.H, red['300K'], red['300K_err'], color = 'red',
                     capsize = 5, marker = 'o', linestyle = '')
        
    def plot_sims(self):
        """
        Plot simulated lines.
        Returns
        -------
        None.

        """
        plt.plot(self.sim_lines.index, self.sim_lines['20K'].values, color = 'blue')
        plt.plot(self.sim_lines.index, self.sim_lines['50K'].values, color = 'black')
        plt.plot(self.sim_lines.index, self.sim_lines['300K'].values, color = 'red')

        
    def fit_langevin(self, tp = 'total'): 
        """
        Langevin fit to magnetic moments data with respect to the magnetic field.
        Parameters
        ----------
        tp : str, optional
            defines the moments to fir (total, spin, orbit). The default is 'total'.

        Returns
        -------
        None.

        """
        df = self.moments.loc[self.moments.type == tp+'_moment']
        df.drop('type', inplace = True, axis = 1)
        self.mag300, self.cov300 = curve_fit(lv, df.H + 1e-6, df['300K'], p0 = (0.2, 1e-20, 300), 
                             bounds = ((0, 1e-24, 300 - 1e-6), (3, 1e-16, 300 +1e-6)))
        
        self.mag50, self.cov50 = curve_fit(lv, df.H + 1e-6, df['50K'], p0 = (0.2, 1e-20, 50), 
                             bounds = ((0, 1e-24, 50 - 1e-6), (3, 1e-16, 50 +1e-6)))
        
        self.mag20, self.cov20 = curve_fit(lv, df.H + 1e-6, df['20K'], p0 = (0.2, 1e-20, 20), 
                             bounds = ((0, 1e-24, 20 - 1e-6), (3, 1e-16, 20 +1e-6)))
        
    def construct_fit_lines(self, rng = [0, 6], points = 1000):
        """
        generates the simulated Langvin lines.
        Parameters
        ----------
        rng : list of floats, optional
            The default is [0, 6].
        points : integer, optional
            The default is 1000.

        Returns
        -------
        None.

        """
        scale = np.linspace(rng[0], rng[1], points)
        
        lines_dict = {}
        lines_dict['20K'] = lv(scale+1e-6, self.mag20[0], self.mag20[1], self.mag20[2])
        lines_dict['50K'] = lv(scale+1e-6, self.mag50[0], self.mag50[1], self.mag50[2])
        lines_dict['300K'] = lv(scale+1e-6, self.mag300[0], self.mag300[1], self.mag300[2])
        
        self.sim_lines = pd.DataFrame(index = scale, columns = ['20K', '50K', '300K'])
        for col in self.sim_lines:
            self.sim_lines[col] = lines_dict[col]

        
class MomentsPerNP:       
    
    def __init__(self, size):
        """
        Parameters
        ----------
        size : float
            edge length of nanoparticle

        Returns
        -------
        None.

        """
        self.size = size
        self.edge = size * 1e-9
        
    def set_volume(self):
        """
        Calculate the volume of the nanoparticle.
        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        self.volume = np.sqrt(2)/3 * (self.edge)**3
        return self.volume
    
    def set_unit_cell(self, a):
        """
        Calculate volume of unit cell.
        Parameters
        ----------
        a : float
            Lattice parameter in Angstrom.

        Returns
        -------
        float
            DESCRIPTION.

        """
        self.uc = (a * 1e-10)**3
        return self.uc
    
    def set_interface(self, thickness = 1):
        """
        Calculate the atoms within the interface.
        Parameters
        ----------
        thickness : float, optional
            Thickness in nanometer. The default is 1.

        Returns
        -------
        
            DESCRIPTION.

        """
        self.shell = thickness * 1e-9       
        if self.shell == 0:
            self.vol_interface = self.volume
        else:
            self.vol_interface = np.sqrt(2)/3 * (self.edge**3 - (self.edge**3 - self.shell)**3)
            
        self.n_atoms = self.vol_interface / self.uc
        return np.round(self.n_atoms)
    
    def get_moments(self, ms):
        """
        Calculate the magnetic moments of a nanoparticles.
        Parameters
        ----------
        ms : float
            number of magnetic moments per absorber atom.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        self.ms = ms
        self.moments_per_np = self.ms * self.n_atoms
        return self.moments_per_np
    
    
class RheedData:
    
    def __init__(self, file):
        self.file = file
        
    def load_data(self, delimiter):
        self.data = pd.read_csv(self.file, delimiter = delimiter, index_col = 'X')
        return self.data
    
    def set_zero(self, shift):
        self.shift = shift
        self.data.index = self.data.index - self.shift
    
    def set_temp(self):
        string = self.file.split('.')[0].split('_')[-1]
        self.data['temperature (K)'] = string
        
    def px_to_2theta(self):
        twotheta = np.arctan((self.data.index * 0.2645) / 341.6)
        self.data['2theta'] = twotheta
        
    def set_energy(self, energy):
        self.energy = energy 
        
        
    def get_wavelength(self):
        c = nat.physical_constants['speed of light in vacuum'][0]
        m = nat.physical_constants['electron mass'][0]
        e = nat.physical_constants['elementary charge'][0]
        h = nat.physical_constants['Planck constant'][0]
        
        self.wavelength = h / (2 * m * e * self.energy * (1 + e * self.energy / (2 * m * c**2)))**0.5
        
    def set_q(self):
        self.q = 4 * np.pi * np.sin(self.data['2theta']/2) / (self.wavelength / 1e-10)
        self.data['q'] = self.q
        
        
        
class RheedSim:
    
    def __init__(self, sim_files, names):
        self.sim_files = {}
        
        for name, file in zip(names, sim_files):
            self.sim_files[name] = pd.read_csv(file, delimiter = '\s+')
            
        
        
        
        
        
        
        
        
        
        
        
        