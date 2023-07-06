import numpy as np
import os.path
import sys
from astropy.io import fits
from astropy.table import Table


class convert_data():
    """
    Converte #.rvo.dat and #.par.dat file to fits file.
    Combine both data files to one fits file- HDU[1] contains rvo data, HDU[2] par data.
    Use CPL for writing to full fill ESO standards. 
    """
    
    def __init__(self, filename, dat=1, fits=0, cpl=0, final=0):
        
        self.filename = filename     
        if fits or cpl:
            self.data_rvo = np.genfromtxt(filename+'.rvo.dat', dtype=None, names=True,
                deletechars='', encoding=None).view(np.recarray)
            self.header_rvo = self.data_rvo.dtype.names

            self.data_par = np.genfromtxt(filename+'.par.dat', dtype=None, names=True,
                deletechars='', encoding=None).view(np.recarray)
            self.header_par = self.data_par.dtype.names
        
        if fits:
            self.write_fits()
        if cpl:
            self.write_cpl()
        if final:
            self.write_finalRV()            
        if not dat:
            os.remove(filename+'.rvo.dat') 
            os.remove(filename+'.par.dat')      
        
    def write_fits(self):
        # convert data to fits files using astropy       
        hdu0 = fits.PrimaryHDU()
        
        # create rvo table
        if self.data_rvo.size == 1:
            self.data_rvo = np.array([self.data_rvo])
                    
        table = []   
        for h, col in enumerate(self.header_rvo):
            if col == 'filename':
                c = fits.Column(name=str(col), array=self.data_rvo[col], format='50A')
            else:    
                c = fits.Column(name=str(col), array=self.data_rvo[col], format='F')
            table.append(c)
        hdr = fits.Header()
        hdr.set('TYPE', 'rvo', 'RV data')
        hdr.set('UNIT', 'm/s', 'Unit of RV data')
        table_rvo = fits.BinTableHDU.from_columns(table)
        
        # create par table
        table = []   
        for h, col in enumerate(self.header_par):
            c = fits.Column(name=str(col), array=self.data_par[col], format='F')
            table.append(c)
        hdr = fits.Header()
        hdr.set('TYPE', 'par', 'parameter data')
        hdr.set('UNIT', 'm/s', 'Unit of  RV data')
        table_par = fits.BinTableHDU.from_columns(table)

        # combine and write to file
        hdul = fits.HDUList([hdu0, table_rvo, table_par])
        hdul.writeto(self.filename+'_rvo_par.fits', overwrite=True)
    
    def write_cpl(self):
        # convert data to fits files using ESO PyCPL libaries
        import cpl
        from cpl.core import Table
        from cpl.core import PropertyList, Property

        hdr = cpl.core.PropertyList()
        tbl = cpl.core.Table(self.data_rvo.size)

        if self.data_rvo.size == 1:
            self.data_rvo = np.array([self.data_rvo])
        
        hdr.append(Property('TYPE', 'rvo', 'RV data'))
        hdr.append(Property('UNIT', 'm/s', 'Unit of RV data'))

        # save rvo data
        for h, col in enumerate(self.header_rvo):
           if str(col) == 'filename':
               tbl.new_column(str(col), cpl.core.Type.STRING)
           else:     
               tbl.new_column(str(col), cpl.core.Type.DOUBLE)
           tbl[str(col)] = self.data_rvo[col]

        Table.save(tbl, None, hdr, self.filename+'_rvo_par.fits', cpl.core.io.CREATE)
 
        # save par data       
        hdr = cpl.core.PropertyList()
        tbl = cpl.core.Table(self.data_par.size)
        
        hdr.append(Property('TYPE', 'par', 'parameter data'))
        hdr.append(Property('UNIT', 'm/s', 'Unit of RV data'))

        for h, col in enumerate(self.header_par):
           tbl.new_column(str(col), cpl.core.Type.DOUBLE)
           tbl[str(col)] = self.data_par[col]

        Table.save(tbl, None, hdr, self.filename+'_rvo_par.fits', cpl.core.io.EXTEND) 

    def write_finalRV(self):
        # convert tmp.dat to fits
                
        data = np.genfromtxt(self.filename+'.dat', dtype=None, names=True,
            deletechars='', encoding=None).view(np.recarray)

        tab = Table(data, names=data.dtype.names)
        tab.write(self.filename+'.fits', format='fits', overwrite=True)
       
