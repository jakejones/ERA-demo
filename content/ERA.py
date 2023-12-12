import nbformat
import ipywidgets as widgets
import os
import numpy as np
import astropy.io.fits as fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import pandas
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

# -- Constants --
c = 299792458.0


class ERAData:

    def __init__(self,file_path):
        
        # -- Load FITS File --
        hdu=fits.open(file_path, memmap=False)
        
        # -- FITS Header Convenience Functions --
        NAXIS  = hdu[0].header["NAXIS"]
        PCOUNT = hdu[0].header["PCOUNT"]
        CTYPE  = [""]*NAXIS
        for i in range(len(CTYPE)):
            try:
                CTYPE[i] = hdu[0].header["CTYPE"+str(i+1)]
            except:
                pass
            
        def CRVAL(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["CRVAL"+str(i+1)]
        
        def CDELT(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["CDELT"+str(i+1)]
        
        def AxisSize(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return hdu[0].header["NAXIS"+str(i+1)]
        
        def AxisIndex(name):
            for i in range(len(CTYPE)):
                if name == CTYPE[i]:
                    return int( hdu[0].header["NAXIS"] - (i+1) )
        
        # -- Telescope Info --
        Nant = int( hdu[1].header["NAXIS2"] )   # Number of antenna's
        Nbls = int( Nant * (Nant+1) / 2 )       # Number of baselines
        
        # -- Observation Info & Data --
        Npol  = AxisSize("STOKES")                        # Number of polarisations
        Nfreq = AxisSize("FREQ")                          # Number of frequency bins
        freq  = np.arange(Nfreq)*CDELT("FREQ") + CRVAL("FREQ") # Frequency Axis
        freq_center = freq[int(Nfreq/2)]                  # Center freq to use after averaging
        Ntimes = len(np.unique(hdu[0].data[:]['DATE']))   # Number of unique timesteps in the data
        Nblts  = len(hdu[0].data)                         # Total number of baselines for all timesteps
        if Nblts / Nbls != Ntimes:
            print("Beware: Not all baselines exist for every timestep")
        
        OBSRA   = CRVAL("RA")                             # Phase Center (Right Ascension)
        OBSDEC  = CRVAL("DEC")                            # Phase Center (Declination)
        
        UVWs = np.zeros([Nblts,3])
        UVWs[:,0] = hdu[0].data[:].par('UU') * freq_center # Mult by freq to get UVW's in (Number of wavelengths)
        UVWs[:,1] = hdu[0].data[:].par('VV') * freq_center # Mult by freq to get UVW's in (Number of wavelengths)
        UVdist = np.sqrt(np.sum( UVWs[:,0:2]**2 , axis=1)) # Compute baseline lengths (Used in figure & to compute image resolution)
                
        print(Nant)
        print(Nbls)
        print(Nfreq)
        print(freq)
        print(freq_center)
        print(Ntimes)
        print(Nblts)
        print(OBSRA)
        print(OBSDEC)

class FileUploader:
    def __init__(self):
        self.uploader = widgets.FileUpload(
            accept='',
            multiple=False
        )
        display(self.uploader)
    
    def write_file(self,file_path="telescope-data.uvfits"):
        uploaded_file = self.uploader.value[0]
        with open(file_path, 'wb') as f: 
            f.write(uploaded_file.content)


