#! /usr/bin/env python3
# Licensed under a GPLv3 style license - see LICENSE


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
    
    def __init__(self, filename, args, dat=1, fits=0, cpl=0, final=0):
        
        self.filename = filename     
        if fits or cpl:
            self.data_rvo = np.genfromtxt(filename+'.rvo.dat', dtype=None, names=True,
                deletechars='', encoding=None).view(np.recarray)
            self.header_rvo = self.data_rvo.dtype.names

            self.data_par = np.genfromtxt(filename+'.par.dat', dtype=None, names=True,
                deletechars='', encoding=None).view(np.recarray)
            self.header_par = self.data_par.dtype.names
            
        self.args = vars(args)
        self.args.pop('config_file', None)
        self.args.pop('ftsname', None)
        self.args = {key: self.args[key] for key in self.args if 'look' not in key}           
        if self.args['tplname']: self.args['tplname'] = self.args['tplname'].split('/')[-1] 
        
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
        hdr = fits.Header()
        hdr.set('REC ID', 'viper RV estimation', 'Pipeline recipe')
        for i, (k, v) in enumerate(self.args.items()):
            hdr.set('REC PARAM'+str(i)+' NAME', k)
            hdr.set('REC PARAM'+str(i)+' VALUE', str(v))
        hdu0 = fits.PrimaryHDU(header=hdr)
        
        # create rvo table
        if self.data_rvo.size == 1:
            self.data_rvo = np.array([self.data_rvo])
                    
        table = []   
        for h, col in enumerate(self.header_rvo):
            unitpar = ''
            if col == 'filename':
                c = fits.Column(name=str(col), array=self.data_rvo[col], format='50A')
            else:
                if 'rv' in str(col).casefold() and str(col) != 'BERV':
                    unitpar = 'm/s'
                c = fits.Column(name=str(col), array=self.data_rvo[col], unit=unitpar, format='F')
            table.append(c)
        hdr = fits.Header()
        hdr.set('TYPE', 'rvo', 'RV data')
        table_rvo = fits.BinTableHDU.from_columns(table)
        
        # create par table
        table = []   
        for h, col in enumerate(self.header_par):
            unitpar = '' 
            if 'rv' in str(col).casefold() and str(col) != 'BERV':
                unitpar = 'm/s'
            elif 'wave' in str(col):
                unitpar = 'A/px^'+str(str(col)[-1])
            c = fits.Column(name=str(col), array=self.data_par[col], unit=unitpar, format='F')
            table.append(c)
        hdr = fits.Header()
        hdr.set('TYPE', 'par', 'parameter data')
        table_par = fits.BinTableHDU.from_columns(table)

        # combine and write to file
        hdul = fits.HDUList([hdu0, table_rvo, table_par])
        hdul.writeto(self.filename+'_rvo_par.fits', overwrite=True)
    
    def write_cpl(self):
        # convert data to fits files using ESO PyCPL libaries
        import cpl
        from cpl.core import Table
        from cpl.core import PropertyList, Property
        
        hdr0 = cpl.core.PropertyList()
        hdr0.append(Property('REC ID', 'viper RV estimation', 'Pipeline recipe'))
        for i, (k, v) in enumerate(self.args.items()):
            hdr0.append(Property('REC PARAM'+str(i)+' NAME', k))
            hdr0.append(Property('REC PARAM'+str(i)+' VALUE', str(v)))

        hdr = cpl.core.PropertyList()
        tbl = cpl.core.Table(self.data_rvo.size)

        if self.data_rvo.size == 1:
            self.data_rvo = np.array([self.data_rvo])
        
        hdr.append(Property('TYPE', 'rvo', 'RV data'))

        # save rvo data
        for h, col in enumerate(self.header_rvo):
           if str(col) == 'filename':
               tbl.new_column(str(col), cpl.core.Type.STRING)
           else:     
               tbl.new_column(str(col), cpl.core.Type.DOUBLE)
               if 'rv' in str(col).casefold() and str(col) != 'BERV':
                   tbl.set_column_unit(str(col), 'm/s')
           tbl[str(col)] = self.data_rvo[col]

        Table.save(tbl, hdr0, hdr, self.filename+'_rvo_par.fits', cpl.core.io.CREATE)
 
        # save par data       
        hdr = cpl.core.PropertyList()
        tbl = cpl.core.Table(self.data_par.size)
        
        hdr.append(Property('TYPE', 'par', 'parameter data'))

        for h, col in enumerate(self.header_par):
           tbl.new_column(str(col), cpl.core.Type.DOUBLE)
           tbl[str(col)] = self.data_par[col]

           if 'rv' in str(col).casefold() and str(col) != 'BERV':
               tbl.set_column_unit(str(col), 'm/s')
           elif 'wave' in str(col):
               tbl.set_column_unit(str(col), 'A/px^'+str(str(col)[-1]))

        Table.save(tbl, None, hdr, self.filename+'_rvo_par.fits', cpl.core.io.EXTEND) 

    def write_finalRV(self):
        # convert tmp.dat to fits
                
        data = np.genfromtxt(self.filename+'.dat', dtype=None, names=True,
            deletechars='', encoding=None).view(np.recarray)

        tab = Table(data, names=data.dtype.names)
        tab.write(self.filename+'.fits', format='fits', overwrite=True)
       
