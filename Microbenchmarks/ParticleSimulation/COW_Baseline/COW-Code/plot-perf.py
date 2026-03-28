#
# A quick and dirty plotting script for the LRZ Code Optimization Course.
#
# Requires Python3, Matplotlib 3, and pandas.
#
# On LRZ systems, run "module load anaconda3/2021.05"
#
#
# Usage:
#
#   python plot-perf.py <nbody-perf-files>
#
# Output is written to COW-perf-plot.pdf
#
# Data is plotted in the order it is given on the command line.
#

#{{{
# Python imports
import os, sys, pickle, tempfile, errno
import shutil
import scipy as sp
from scipy.interpolate import interp1d
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as pl
import itertools as it
from matplotlib import rc,rcParams,cm
from matplotlib import ticker as mticker
import subprocess
from collections import defaultdict, OrderedDict
from itertools import cycle,product
from datetime import datetime
from scipy.optimize import minimize
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
#}}}

mpl.style.use('seaborn-dark-palette')
display=False

#{{{ 
# Display properties and figure saving routines
class C:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ERROR = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

jls = {
     'loosely dotted'		: (0, (1, 10)),
     'dotted'			: (0, (1, 1)),
     'densely dotted'		: (0, (1, 1)),
     'loosely dashed'		: (0, (5, 10)),
     'dashed'			: (0, (5, 5)),
     'densely dashed'		: (0, (5, 1)),
     'loosely dashdotted'	: (0, (3, 10, 1, 10)),
     'dashdotted'		: (0, (3, 5, 1, 5)),
     'densely dashdotted'	: (0, (3, 1, 1, 1)),
     'dashdotdotted'		: (0, (3, 5, 1, 5, 1, 5)),
     'loosely dashdotdotted'	: (0, (3, 10, 1, 10, 1, 10)),
     'densely dashdotdotted'	: (0, (3, 1, 1, 1, 1, 1))
     }

if not display:
    rc('text', usetex=False)
    rc('font',**{'family':'serif','serif':['Palatino']})
    #rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('font',size=8)
    #rc('text',fontsize=14)
    rc('axes',labelsize=10)
    rc('axes',linewidth=0.1)
    rc('xtick',labelsize=8)
    rc('ytick',labelsize=8)
    rc('legend',fontsize=5)
    rc('lines',linewidth=0.1)
    figopts  = dict(figsize=np.array([4.0,1.5]) * 172.5/72.27,dpi=72.27)
    saveopts = dict(dpi=72.27, bbox_inches='tight')
    legopts  = dict(frameon=False, prop={'size':8})
else:
    #rc('text', usetex=True)
    #rc('font',**{'family':'serif','serif':['Palatino']})
    #rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('font',size=18)
    rc('text',fontsize=18)
    rc('axes',labelsize=18)
    rc('xtick',labelsize=14)
    rc('ytick',labelsize=14)
    rc('legend',fontsize=10)
    rc('lines',linewidth=0.1)
    rc('axes',linewidth=0.1)
    figopts  = dict() #figsize=np.array([2,1]) * 172.5/72.27,dpi=72.27)
    saveopts = dict(dpi=72.27, bbox_inches='tight')
    legopts  = dict(frameon=False, prop={'size':5})

#}}}

#{{{
def load_perf_data(fname):

    table = None

    if 1:
        try:
            table = OrderedDict()
            with open(fname, 'r') as fp:
                print('Reading', fname)
                npart     = 0
                compiler  = None
                flags     = None
                npart     = None
                nranks    = 1
                nthreads  = 1
                tile_size = 0
                keys      = None
                units     = None
                tag       = None
                for lineno,l in enumerate(fp):
                    l = l.strip()
                    if not l: continue
                    if l.startswith('#'): continue

                    tok = l.split()

                    kind = tok[0]

                    if kind == 'META':
                        if tok[1] == 'NEW':
                            pass

                    if kind == 'META':
                        if tok[1] == 'Compiler:':
                            compiler = tok[2]
                        elif tok[1] == 'Flags:':
                            flags = list(filter(None, tok[2].strip().split(',')))
                        elif tok[1] == 'Npart:':
                            npart = int(tok[2])
                        elif tok[1] == 'Nranks:':
                            nranks = int(tok[2])
                        elif tok[1] == 'Nthreads:':
                            nthreads = int(tok[2])
                        elif tok[1] == 'TileSize:':
                            tile_size = int(tok[2])
                        elif tok[1] == 'Tag:':
                            tag = tok[2]

                    if kind == 'HEAD':
                        keys = l[len(kind):].strip().split(',')
                    elif kind == 'UNIT':
                        units = l[len(kind):].strip().split(',')

                    if kind == 'DATA':
                        assert(npart > 0)
                        assert(compiler is not None)
                        assert(flags is not None)
                        assert(npart is not None)
                        assert(keys is not None)
                        assert(units is not None)

                        if nthreads > 2 and tile_size == 0:
                            print(C.WARNING + "nthreads > 2 and tile_size == 0 in %s"%fname + C.ENDC)

                        data = [tok[1]] + list(map(float,tok[2:]))

                        table_keys = ['Compiler', 'Npart', 'Nranks', 'Nthreads', 'TileSize', 'Tag'] + flags + keys
                        table_data = [compiler, npart, nranks, nthreads, tile_size, tag] + [True]*len(flags) + data

                        for k,d in zip(table_keys, table_data):
                            table[k] = table.get(k,[]) + [d]

                        good = True

        except Exception as e:
            print(C.ERROR + str(e) + C.ENDC)

    return pd.DataFrame(table)

def create_dataframe(fnames):

    table = pd.DataFrame()
    for fname in fnames:
        df = load_perf_data(fname) 
        table = table.append(df, ignore_index=True)

    table = table.fillna(False)
    return table

#}}}

def make_plots(df, pdf):

    if 1:
        f,a = pl.subplots(1,1, sharex=False, squeeze=False, gridspec_kw={'hspace': 0}, **figopts)
        ax = a[0,0]
        pl.sca(ax)

        xlabels = set()
        all_keys = ['soa','O3', 'arch', 'float', 'restrict','simd','align','malloc_align','Nthreads', 'openmp', 'first_touch', 'tiles', 'TileSize']

        offs = 0
        for f,marker in zip(['Kick', 'Drift', 'Accel'], ['x', 'o', 's']):
            D = df.query(f'Kind == @f')
            offs = 0
            for by0,Y0 in D.groupby(by=['Compiler', 'Npart'], sort=False):
                keys = [x for x in all_keys if x in Y0]
                #Y = Y0.sort_values(by=keys)
                Y = Y0
                pd.options.display.float_format = '{:,.4e}'.format
                print(Y[keys])

                ys = Y['FLOPS'].to_list()
                xs = np.arange(len(ys)) + offs
                offs += len(ys)

                for x,l in zip(xs,Y['Tag'].to_list()):
                    if l: xlabels.update(((x,l),))

                c = 'blue' if Y['Compiler'].iloc[0] == 'Intel' else 'green'
                ax.plot(xs,ys,c=c, marker=marker, mfc='none', markersize=8, markeredgewidth=0.5)

        ax.grid(which='both', axis='y', zorder=-1000, c='lightgrey')
        ax.set_yscale('log', base=10)
        ax.set_ylabel('FLOPS')

        xt,xl = list(zip(*xlabels))
        ax.set_xticks(xt)
        ax.set_xticklabels(xl, rotation=45, ha='right')

        legend_line_entry = lambda name, **args: [name, mpl.lines.Line2D((0,0.5), (1,0.5), mfc='none', markersize=4, markeredgewidth=0.5, visible=True, **args)]

        lh = [
            legend_line_entry('Kick',   marker='x',   c='k', lw=1, ls='-'),
            legend_line_entry('Drift',  marker='o',   c='k', lw=1, ls='-'),
            legend_line_entry('Accel',  marker='s',   c='k', lw=1, ls='-'),
        ]
        l,h = list(map(list, list(zip(*lh))))
        L = ax.legend(h,l, ncol=2, numpoints=1, loc='upper left', columnspacing=0.5, **dict(legopts, **dict(frameon=True, prop={'size':7})))
        L.get_frame().set_lw(0.5)
        L.get_frame().set_alpha(1.0)

        pdf.savefig(pl.gcf(), **saveopts)
        pl.close(pl.gcf())

    if 1:
        f,a = pl.subplots(1,1, sharex=False, squeeze=False, gridspec_kw={'hspace': 0}, **figopts)
        ax = a[0,0]
        pl.sca(ax)

        xlabels = set()

        offs = 0
        for f,marker in zip(['Kick', 'Drift', 'Accel'], ['x', 'o', 's']):
            D = df.query(f'Kind == @f')
            offs = 0
            for by0,Y0 in D.groupby(by=['Compiler', 'Npart'], sort=False):
                all_keys = ['soa','O3', 'arch', 'float', 'restrict','simd','align','malloc_align','Nthreads', 'openmp', 'first_touch', 'tiles', 'TileSize']
                #all_keys = ['soa']
                keys = [x for x in all_keys if x in Y0]
                #Y = Y0.sort_values(by=keys)
                Y = Y0
                pd.options.display.float_format = '{:,.4e}'.format
                print(Y[keys])
                #print(Y[['Npart','Compiler', 'Kind','Nthreads','FLOPS','O3', 'arch', 'float', 'restrict','simd','malloc_align','align', 'openmp', 'first_touch', 'tiles', 'TileSize']])
                #for by,Y in Y0.groupby(by=['Npart'], sort=True):

                #print(Y[['Npart','restrict','align','malloc_align','float']])
                #print(Y)
                ys = Y['Time'].to_list()
                #print(Y['Kind'], by, ys)
                xs = np.arange(len(ys)) + offs
                offs += len(ys)

                for x,l in zip(xs,Y['Tag'].to_list()):
                    if l: xlabels.update(((x,l),))

                c = 'blue' if Y['Compiler'].iloc[0] == 'Intel' else 'green'
                ax.plot(xs,ys,c=c, marker=marker, mfc='none', markersize=8, markeredgewidth=0.5)

        ax.grid(which='both', axis='y', zorder=-1000, c='lightgrey')
        ax.set_yscale('log', base=10)
        ax.set_ylabel('Time [s]')

        xt,xl = list(zip(*xlabels))
        ax.set_xticks(xt)
        ax.set_xticklabels(xl, rotation=45, ha='right')

        legend_line_entry = lambda name, **args: [name, mpl.lines.Line2D((0,0.5), (1,0.5), mfc='none', markersize=4, markeredgewidth=0.5, visible=True, **args)]

        lh = [
            legend_line_entry('Kick',   marker='x',   c='k', lw=1, ls='-'),
            legend_line_entry('Drift',  marker='o',   c='k', lw=1, ls='-'),
            legend_line_entry('Accel',  marker='s',   c='k', lw=1, ls='-'),
        ]
        l,h = list(map(list, list(zip(*lh))))
        L = ax.legend(h,l, ncol=2, numpoints=1, loc='best', columnspacing=0.5, **dict(legopts, **dict(frameon=True, prop={'size':7})))
        L.get_frame().set_lw(0.5)
        L.get_frame().set_alpha(1.0)

        pdf.savefig(pl.gcf(), **saveopts)
        pl.close(pl.gcf())


def main():

    df = create_dataframe(sys.argv[1:])
    with PdfPages('COW-perf-plot.pdf') as pdf:
    	make_plots(df, pdf)

if __name__ == '__main__':
    main()

