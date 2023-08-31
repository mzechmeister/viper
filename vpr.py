#! /usr/bin/env python3

# ./vpr.py

import argparse
import os
import numpy as np
from astropy.io import fits

from wstat import wsem

from pause import pause
from gplot import gplot

gplot.bar(0.5)

def plot_RV(file):
    gplot.mxtics().mytics()
    gplot.xlabel("'BJD - 2 450 000'")
    gplot.ylabel("'RV [m/s]'")
    gplot.key("title '%s' noenhance" % file)
    gplot('"%s" us ($1-2450000):"RV":(sprintf("%%s\\nn: %%d\\nBJD: %%.6f\\nRV: %%f +/- %%f", strcol("filename"),$0+1,$1,$2,$3)) with labels hypertext point pt 0 t"", "" us ($1-2450000):"RV":"e_RV" w e pt 7 lc "#77000000"' % file)
    pause()

def plot_rvo(rv=None, e_rv=None, file=None):
    #gplot('"%s" matrix every 2::1 us (($1-1)/2+18):3 w p pt 7' % file)
    if file:
        mat = np.genfromtxt(file)
        bjd, RV, e_RV = mat[:3]
        rv, e_rv = mat[3::2], mat[4::2]

    ii = np.isfinite(e_rv)
    RV = np.mean(rv[ii])
    e_RV = np.std(rv[ii])/(ii.sum()-1)**0.5

    gplot.xlabel('"Order index"')
    gplot.ylabel('"RV [m/s]"')
    gplot(rv, e_rv, 'us 0:1:2 w e pt 7 t "%s", RV=%s, e_RV=%s, RV lc 3 t "RV = %.5f +/- %.5f m/s", RV+e_RV lc 3 dt 2 t "", RV-e_RV lc 3 dt 2 t ""' % (file, RV,e_RV,RV,e_RV))
    pause()

def arg2slice(arg):
    """Convert string argument to a slice."""
    # We want four cases for indexing: None, int, list of ints, slices.
    # Use [] as default, so 'in' can be used.
    if isinstance(arg, str):
        arg = eval('np.s_['+arg+']')
    return [arg] if isinstance(arg, int) else arg

def plot_cmp(vpr, vprcmp):
        gplot.mxtics().mytics()
        gplot.xlabel("'BJD - 2 450 000'")
        gplot.ylabel("'RV [m/s]'")
        gplot.key('invert')
        gplot.bind('''"$" "i=!i; set xlabel i? 'BJD - 2 450 000':'spectrum number n'; repl";  i=1''')

        gplot-(vprcmp.BJD, vprcmp.RV, vprcmp.e_RV, vprcmp.A.filename, ' us (i?$1-2450000:$0+1):2:(sprintf("%%s\\nn: %%d\\nBJD: %%.6f\\nRV: %%f +/- %%f", stringcolumn(4),$0+1,$1,$2,$3)) with labels  hypertext point pt 0 t "", "" us (i?$1-2450000:$0+1):2:3 w e pt 7 lc "#55FF0000" t "%s [o=%s] rms=%.2f m/s (med(σ)=%.2f m/s)" noenh' % (vprcmp.tag, str(vprcmp.oset).replace('\n',''), vprcmp.rms, vprcmp.medunc))
        gplot+(vpr.BJD, vpr.RV, vpr.e_RV, vpr.A.filename, ' us (i?$1-2450000:$0+1):2:(sprintf("%%s\\nn: %%d\\nBJD: %%.6f\\nRV: %%f +/- %%f", stringcolumn(4),$0+1,$1,$2,$3)) with labels  hypertext point pt 0 t "", "" us (i?$1-2450000:$0+1):2:3 w e pt 7 lc "#3300000" t "%s [o=%s] rms=%.2f m/s (med(σ)=%.2f m/s)" ' % (vpr.tag, str(vpr.oset).replace('\n',''), vpr.rms, vpr.medunc))
        pause('RV time serie')

def average(yi, e_yi=None, typ='wmean', **kwargs):
    # other futures types might be median, trimmed mean, maximum jitter likelihood mean
    if typ == 'mean':
        Y = np.mean(yi, **kwargs)
        e_Y = np.std(yi, **kwargs) / (yi.size//Y.size-1)**0.5
    else:
        e_yi[np.isnan(e_yi)] = np.inf
        Y, e_Y = wsem(yi, e=e_yi, **kwargs)
    return Y, e_Y

class VPR():
    def __init__(self, tag, oset=None, ocen=None, avg='wmean', gp='', sort='', cen=False, offset = 400):
        '''
        oset: slice,list
        '''
        if tag.endswith('.fits'):
            self.file = file = tag
            self.tag = tag.replace('.fits', '')    
        else:    
            self.tag = tag.replace('.rvo.dat', '')
            self.file = file = self.tag + '.rvo.dat'  

        self.oset = oset
        self.avgtyp = avg
        self.offset = offset
        print(self.tag)
        
        if not os.path.isfile(self.file):
            return        

        if gp:
           gplot.put(gp)

        if self.file.endswith('.fits'):
            hdu = fits.open(self.file, dtype=None)
            self.Afull = hdu[1].data.view(np.recarray)
            mat = hdu[1].data      
        else:
            try:
                self.Afull = np.genfromtxt(file, dtype=None, names=True,
                deletechars='',   # to keep the dash for chunks
                encoding=None).view(np.recarray)
            except:
                self.Afull = np.genfromtxt(file, dtype=None,
                deletechars='',   # to keep the dash for chunks
                names=True).view(np.recarray)    
            mat = np.genfromtxt(file, skip_header=1)  

            if np.size(self.Afull) == 1:
                mat = np.array([mat])
                self.Afull = (np.array([self.Afull])).view(np.recarray) 

        if sort:
            self.Afull.sort(order=sort)
        colnames = self.Afull.dtype.names

        onames_all = self.onames = np.array([col for col in colnames if col.startswith('rv')], dtype='O')
        orders_all, chunks_all = [*np.array([[*map(int, oname[2:].split('-'))] for oname in onames_all]).T, None][0:2]
        self.orders = orders_all
        self.chunks = chunks_all
        if oset is not None:
            self.oset = olist = np.r_[oset]
            ofilter = [int(o) in olist for o in orders_all]
            self.onames = self.onames[ofilter]
            self.orders = self.orders[ofilter]
            if self.chunks is not None:
                self.chunks = self.chunks[ofilter]

        # remove columns not in oset; keep BJD, etc. and other orders
        self.A = self.Afull[[col for col in colnames if col.endswith(tuple(self.onames)) or not col.startswith(('e_rv', 'rv'))]]

        self.full_rv, self.full_e_rv = mat.T[4:-1:2], mat.T[5::2]

        self.rv = np.array(self.A[self.onames].tolist()).T
        self.e_rv = np.array(self.A["e_"+self.onames].tolist()).T

        # recompute mean RV
        self.BJD = self.A.BJD
        self.RV, self.e_RV = average(self.rv, self.e_rv, axis=0, typ=self.avgtyp)

        if cen:
            off = np.nanmedian(self.RV)
            print('Subtracting', off, 'm/s from all RVs.')
            self.RV -= off
            self.rv -= off

        if self.file.endswith('.fits'):
            self.filenames = [name.decode('UTF-8') for name in self.A.filename]
        else:
            self.filenames = self.A.filename

        self.info()

        if ocen:
            RVo = np.mean(self.rv-self.RV, axis=1)
            print('Subtracting mean order offsets from all RVs.')
            # This should not change the RV mean values (if unweighted).
            # The idea is to reduce RV uncertainty overestimation coming from large order offsets.
            # Offsets are hints for systematics in wavelength solutions (template or fts-iod)
            self.rv -= RVo[:,None]
            self.RV, self.e_RV = average(self.rv, self.e_rv, axis=0, typ=self.avgtyp)
            self.info()

        self.stat_o = np.percentile(self.rv-self.RV, [17,50,83], axis=1)
        self.med_e_rvo = np.median(self.e_rv, axis=1)

    def plot_par(self, parcolx = None, parcoly = None):

        if self.file.endswith('.fits'):
            hdu = fits.open(self.file, dtype=None)
            par = hdu[2].data.view(np.recarray)
            parfile = self.file   
        else:
            parfile = self.file.replace('.rvo.dat', '.par.dat')
            # for now just to get the colnames; not handled: oset, ocen, cen, ...
            par = np.genfromtxt(parfile, dtype=None, names=True,
                            deletechars='',   # to keep the dash for chunks
                            encoding=None).view(np.recarray)

        gplot.mxtics().mytics()
        colnames = par.dtype.names[:-1]
        gplot.key("title '%s' noenhance" % (parfile))
 
        gplot.var(colx = 0)
        if parcolx:            
             gplot.xlabel("'%s'" % parcolx)
             gplot.var(colx = ['n', 'order', 'BJD'].index(parcolx))
        
        gplot.put(f'colnames(x) = word("{" ".join(colnames)}", x); Ncol={len(colnames)}')
        gplot.bind('''"$" "i=(i+1)%3; set xlabel i==0 ? 'Number n' : i==1? 'Order o' : 'BJD'; repl"; i=colx''')

        if parcoly:
            gplot.ylabel("'%s'" % parcoly)
            gplot.var(coly = colnames.index(parcoly)+1)
         #   gplot(par.BJD, par.n, par.order, par[parcoly], par['e_'+ parcoly], 'us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):4:4+1:(i==1?$2:$3) w e palett')
            gplot(par.BJD, par.n, par.order, par[parcoly], par['e_'+ parcoly], 'us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):4:(sprintf("BJD: %%.6f\\nn: %%d\\norder: %%d", $1,$2,$3))  with labels hypertext point pt 0 t"", "" us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):4:4+1:(i==1?$2:$3)  w e palett t "%s"' % str(parcoly))
        else:
            gplot.bind('''"@" "k = (k+2)%Ncol; set ylabel colnames(k); repl"; k=5''' )
        #gplot(par.BJD, par.n, par.order, par.ip0, par.e_ip0, 'us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):k:k+1:(i==1?$2:$3) w e palett')
        #    gplot(f'"{parfile}"' + 'us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):k:k+1:(i==1?$2:$3) w e palett')
            gplot(f'"{parfile}"' + 'us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):k:(sprintf("BJD: %%.6f\\nn: %%d\\norder: %%d", $1,$2,$3))  with labels hypertext point pt 0 t"", "" us (i==0? $2+$3/50 : i==1?  $3+$2/400 : $1-2450000):k:k+1:(i==1?$2:$3) w e palett t "%s"' % str(parfile))
            print('press "@" to browse the parameters')
        print('press "$" to toggle between n/o/BJD')
        pause('Parameters')

    def info(self):
        print('Number of chunks:', self.orders.size)
        self.rms = np.std(self.RV)
        self.medunc = np.median(self.e_RV)
        print('rms(RV) [m/s]:     ', self.rms)
        print('median(e_RV) [m/s]:', self.medunc)

    def plot_RV(self):
        gplot.mxtics().mytics()
        gplot.xlabel("'BJD - 2 450 000'")
        gplot.ylabel("'RV [m/s]'")
        gplot.key("title '%s' noenhance" % (self.tag))
        gplot.var(cx=1)
        gplot.bind('''"$" "cx = !cx; set xlabel cx? 'BJD - 2 450 000' : 'observation number';  repl"''')

        gplot(self.BJD, self.RV, self.e_RV, self.filenames, ' us (cx?$1-2450000:$0+1):2:(sprintf("%%s\\nn: %%d\\nBJD: %%.6f\\nRV: %%f +/- %%f", stringcolumn(4),$0+1,$1,$2,$3)) with labels  hypertext point pt 0 t"", "" us (cx?$1-2450000:$0+1):2:3 w e pt 7 lc "#77000000" t "orders = %s"' % str(self.oset).replace('\n',''))
        pause('RV time serie')

    def plot_rv(self, o=None, n=1):
        A = self
        gplot.var(n=n, N=len(self.rv.T))
        gplot.key('title "%s" noenhance'%self.tag)
        # print filename on terminal with carriage return and some trailing space
        gplot.macro('''; markspec = "repl; system('printf \\"\\r'.n.': '.spectrum.'       \\"')"''')
        gplot.bind('")" "n = n>=N? N : n+1; @markspec"')
        gplot.bind('"(" "n = n<=1? 1 : n-1;  @markspec"')
        gplot.bind('"]" "n = n+10>N? N : n+10; @markspec"')
        gplot.bind('"[" "n = n-10<1? 1 : n-10; @markspec"')
        gplot.bind('"^" "n = 1; @markspec"')
        gplot.bind('"$" "n = N; @markspec"')
        gplot.xlabel("'order o'")
        gplot.ylabel("'RV_{n,o} -- RV_{n}  [m/s]'")
        gplot.mxtics().mytics()

        Nch = 1
        gap = 0.4   # plot gaps between orders and chunks
        gapch = 0
        if self.chunks is not None:
            Nch = self.chunks.max() + 1
            gapch = 0.2 / (Nch-1)
        chksz = (1 - gap - (Nch-1) *gapch) / Nch
        xpos = self.orders - 0.5 + 0.5*gap + (0 if self.chunks is None else self.chunks * (chksz+gapch))

        gplot.put("replace(x, s1, s2) = (i=strstrt(x, s1), i ? x[1:i-1].s2.replace(x[i+1:], s1, s2) : x)")
        gplot('for [n=1:N]', xpos, (A.rv-A.RV).T, self.e_rv.T,
            f'us ($1+{chksz}*n/N):(column(1+n)):(column(1+n+N)) w e pt 6 lc "light-grey" t "", ' +
            f'"" us ($1+{chksz}*n/N):(column(1+n)):(column(1+n+N)) w e pt 6 lc "red" t "RV_{{".n.",o}} -- RV_{{".n."}}",',
            f'"" us ($1+{chksz}*n/N):(column(1+n)):'+'(sprintf("RV_{n=%d,o=%d} = %.2f ± %.2f m/s", n,$1, column(1+n), column(1+n+N))) w labels hypertext enh point pt 0 lc "red" t "",',
            f'"" us ($1+{chksz}*n/N):(column(1+n)-column(1+n+N)):'+'(sprintf("%.2f   ", column(1+n+N))) w labels noenh rotate right tc "red" t "",',  # red text
            A.BJD, A.RV+self.offset, A.e_RV, A.A.filename, ' us 1:2:(sprintf("%s\\nn: %d\\nBJD: %.6f\\nRV: %f ± %f",strcol(4),$0+1,$1,$2,$3)) w labels hypertext point pt 0 axis x2y1 t "",' +
            '"" us 1:2:3 w e lc "#77000000" pt 7 axis x2y1 t "RV_n",' +
            'spectrum="", "" us 1:(RVn=$2):(spectrum=strcol(4), e_RVn=$3) every ::n-1::n-1 w e lc "red" pt 7 axis x2y1 t replace(spectrum, "_", "\\\\_")." RV_{".n."}=".sprintf("%.2f±%.2f m/s", RVn-400, e_RVn),',
            xpos+chksz/2, A.stat_o, A.med_e_rvo, 'us 1:3:2:4 w e lc "blue" pt 4 t "order stat",' +
            ' "" us 1:3:(sprintf("o = %d\\noffset median: %.2f m/s\\nspread: %.2f m/s\\nmedian error: %.2f m/s", $1, $3, ($4-$2)/2, $5)) w labels hypertext rotate left point pt 0 lc 3 t "",' +
            '"" us 1:4:(sprintf("  %.2f",($4-$2)/2)) w labels rotate left tc "blue" t ""')
        print("Use '()[]^$' in gnuplot window to go through epochs n. Press Enter in terminal to quit.")
        pause('rv order dispersion\n')

    def plot_nrvo(self):
        A = self
        gplot.var(N=len(self.rv.T))
        gplot.key('title "%s" noenhance'%self.tag)
        # print filename on terminal with carriage return and some trailing space
        gplot.xlabel("'spectrum number n'")
        gplot.ylabel("'RV_{n,o} -- RV_{o}  [m/s]'")
        gplot.mxtics().mytics()

        gap = 40   # [m/s] separation between orders
        No = self.orders.size
        # parula http://www.gnuplotting.org/tag/palette + alpha 0x77
        colors = [0x0072bd, # blue
                  0xd95319, # orange
                  0xedb120, # yellow
                  0x7e2f8e, # purple
                  0x77ac30, # green
                  0x4dbeee, # light-blue
                  0xa2142f, # red
                 ]
        gplot.array(colors=colors, orders=self.orders)
        gplot.rmar(20)
        gplot.macro()

        do_lbl = ''
        for i, (o, rms, medunc) in enumerate(zip(self.orders, (self.stat_o[2]-self.stat_o[0])/2, self.med_e_rvo), 1):
            do_lbl += f'set label {i} "-- {o:3d}: {rms:.2f} ({medunc:.2f})" at graph 1, first ({i}+off)*{gap} tc rgb colors[({i}-1)%7+1]; '
        gplot.var(alpha=0x77000000, i=0, off=-No//2, do_lbl=f''' '{do_lbl}' ''')

        gplot.bind('''"$" "i = !i; set xlabel i? 'BJD - 2 450 000':'spectrum number n'; repl"; @do_lbl''')
        gplot.bind('''")" "off = off - 1; @do_lbl; repl"''')
        gplot.bind('''"(" "off = off + 1; @do_lbl; repl"''')

        gplot(f'for [o=1:{No}]', A.BJD, A.rv-A.RV, self.e_rv, f'us (i? $1-2450000:$0+1):(column(o+1)+(o+off)*{gap}):{No}+1+o:(alpha+colors[(o-1)%7+1]) w e lc rgb var t "",',
              A.BJD, A.RV, A.e_RV, A.A.filename, f' us (i? $1-2450000:$0+1):2:3 w e pt 7 lc "#77000000" t "{self.rms:.2f} m/s ({self.medunc:.2f})"')
        pause('nrvo\n')

    def save(self, filename):
        print('saving as', filename)

        if filename.endswith('.fits'):
            from astropy.table import Table
            tab = Table(np.asarray([*zip(self.BJD, self.RV, self.e_RV)]), names=["BJD", "RV", "e_RV"])
            tab.write(filename, format='fits', overwrite=True)
        else:
            np.savetxt(filename, [*zip(self.BJD, self.RV, self.e_RV)], header='BJD RV e_RV')

def plot_res(folder, o=[1], n=[1], sep=1.):
    '''
    Plot stacked residuals.

    Examples:
    vpr -res -oset [20] -nset [`seq -s, 1 50`]
    vpr -res -oset [`seq -s, 19 29`] -nset [17]
    '''
    gplot.var(No=len(o), o=1, obeg=1, oend=len(o))
    gplot.var(Nn=len(n), n=1, nbeg=1, nend=len(n))
    gplot.bind('''"("  'o=obeg=oend = o>1?o-1:1;   nbeg=1; nend=Nn; set key tit "(o=".A[o].")"; set y2label "spectrum number"; repl' ''')
    gplot.bind('''")"  'o=obeg=oend = o<No?o+1:No; nbeg=1; nend=Nn; set key tit "(o=".A[o].")"; set y2label "spectrum number"; repl' ''')
    gplot.bind('''"["  'n=nbeg=nend = n>1?n-1:1;   obeg=1; oend=No; set key tit "[n=".Sp[n]."]"; set y2label "order"; repl' ''')
    gplot.bind('''"]"  'n=nbeg=nend = n<Nn?n+1:Nn; obeg=1; oend=No; set key tit "[n=".Sp[n]."]"; set y2label "order"; repl' ''')
    gplot.key_invert()
    try:
        gplot.term_qt_size('600,1000')
    except:
        pass
    gplot.xlabel("'pixel x'").ylabel("'residuals'").y2label("'spectrum number'")
    gplot.array(A=o, Sp=n)
    gplot.array(d3color=[0x1F77B4, 0xFF7F0E, 0x2CA02C, 0xD62728, 0x9467BD, 0x8C564B])
    gplot.var(alpha=0x77<<24)
    print("type '(' and ')' to go through the orders o or '[' and ']' to go through epochs n")
    gplot.ytics('nomirr').y2tics().format('y2 "%g"').link('y2');
    gplot-(f'for [n=nbeg:nend] for [o=obeg:oend] sprintf("{folder}/%03d_%03d.dat", Sp[int(n)], A[int(o)]) us 1:($2+{sep}*(nbeg==nend?o:n)) lc rgb alpha+d3color[1+ int(nbeg==nend?o:n)%5] pt 7 ps 0.5 t "".Sp[n]."-".A[o]')
    # colored y2ticlabel are not really possible :(
    gplot+(f'[i=1:(nbeg==nend?No:Nn):1] "+" us (NaN):({sep}*$1):y2ticlabel(sprintf("%i", nbeg==nend?A[$1]:Sp[$1])) t ""')
    pause('residuals stacked: o=%s' % o)


def run(cmd=None):
    parser = argparse.ArgumentParser(description='Analyse viper RVs', add_help=False, formatter_class=argparse.RawTextHelpFormatter)
    argopt = parser.add_argument   # function short cut
    argopt('tag', nargs='?', help='tag', default='tmp', type=str)
    argopt('-avg', help='Average typ (mean, weighted mean).', default='wmean', choices=['mean', 'wmean'], type=str)
    argopt('-cen', help='center RVs to zero median', action='store_true')
    argopt('-cmp', help='compare two time series (default: cmp=None or, if cmposet is passed, cmp=tag)', type=str)
    argopt('-cmpocen', help='center orders (subtract order offset) of comparison', action='store_true')
    argopt('-cmposet', help='index for order subset of comparison', type=arg2slice)
    argopt('-gp', help='gnuplot commands', default='', type=str)
    argopt('-nset', help='index for spectrum subset (e.g. 1:10, ::5)', default=None, type=arg2slice)
    argopt('-ocen', help='center orders (subtract order offset)', action='store_true')
    argopt('-offset', help='RV plotting offset for rvo plot', default=400, type=float)
    argopt('-oset', help='index for order subset (e.g. 1:10, ::5)', default=None, type=arg2slice)
    argopt('-parcolx', help='Column name for plot par x-axis', default='', type=str)
    argopt('-parcoly', help='Column name for plot par y-axis', default='', type=str)
    argopt('-plot', help='List of plot tasks', nargs='+', default=['rv', 'rvo'], dest='tasks', choices=['rv', 'rvo', 'nrvo', 'par'])
    argopt('-save', nargs='?', help='Filename to save altered RVs.  Supported formats: .dat and .fits.', const='tmp.dat', metavar='FILENAME')
    argopt('-sort', nargs='?', help='sort by column name', const='BJD')
    argopt('-res', help='Plot residuals stacked (folder name)', nargs='?',  const='res', type=str)
    argopt('-ressep', help='Separation between traces (tip: ~5 rms).', default=1., type=float)

    args = vars(parser.parse_args(cmd))

    resopt = args.pop('res')
    ressep = args.pop('ressep')

    if resopt is not None:
        plot_res(resopt, o=args['oset'], n=args['nset'], sep=ressep)
        return

    args.pop('nset')
    tagcmp = args.pop('cmp')
    cmposet = args.pop('cmposet')
    cmpocen = args.pop('cmpocen')
    tasks = args.pop('tasks')
    savefile = args.pop('save')
    parcolx = args.pop('parcolx')
    parcoly = args.pop('parcoly')

    vpr = VPR(**args)

    if savefile:
        vpr.save(savefile)
        
    gplot.reset()    

    if tagcmp or cmposet or cmpocen:
        if tagcmp: args['tag'] = tagcmp
        if cmposet: args['oset'] = cmposet
        if cmpocen: args['ocen'] = cmpocen
        print()
        vprcmp = VPR(**args)
        if 'rv' in tasks:
            plot_cmp(vpr, vprcmp)
    elif args['oset'] is None and not args['cen'] and not args['ocen']:
        if 'rv' in tasks:
            #plot_RV(vpr.file)
            vpr.plot_RV()
    else:
        if 'rv' in tasks:
            vpr.plot_RV()

    if 'rvo' in tasks:
        vpr.plot_rv(n=1)

    if 'par' in tasks:
        vpr.plot_par(parcolx, parcoly)

    if 'nrvo' in tasks:
        vpr.plot_nrvo()

if __name__ == "__main__":
    run()
