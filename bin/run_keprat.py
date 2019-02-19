#!/usr/bin/env python
import os
import glob
from argparse import ArgumentParser
from collections import OrderedDict

import pandas as pd
from matplotlib import pylab as plt

import keprat.io
def main():
    psr = ArgumentParser()
    subpsr = psr.add_subparsers(title="subcommands", dest='subcommand')
    psr_parent = ArgumentParser(add_help=False)

    psr2 = subpsr.add_parser('cksgaia-planets', parents=[psr_parent],)
    psr2.set_defaults(func=cksgaia_planets)

    psr2 = subpsr.add_parser('create-plot', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_plot)

    psr2 = subpsr.add_parser('create-table', parents=[psr_parent], )
    psr2.add_argument('name',type=str)
    psr2.set_defaults(func=create_table)

    psr2 = subpsr.add_parser('update-paper', parents=[psr_parent])
    psr2.set_defaults(func=update_paper)

    args = psr.parse_args()
    args.func(args)

def cksgaia_planets(args):
    import keprat.io
    df = keprat.io.load_table('cksgaia-planets')
    fn = 'data/cksgaia-planets.csv'
    df.to_csv(fn)
    print "created {}".format(fn)

def create_table(args):
    w = Workflow()
    w.create_file('table', args.name ) 

def create_plot(args):
    w = Workflow()
    w.create_file('plot', args.name ) 

def create_val(args):
    w = Workflow()
    w.create_file('val',args.name) 

def update_paper(args):
    w = Workflow()
    w.update_paper()

class Workflow(object):
    def __init__(self):
        d = OrderedDict()
        # register different plots here
        self.plot_dict = d

        d = OrderedDict()
        # register different tables here
        #d['star'] = lambda : ckscool.table.tab_star() 
        self.table_dict = d

        d = OrderedDict()
        #d['stat'] = ckscool.value.val_stat
        self.val_dict = d

        d = OrderedDict()
        d['table'] = self.table_dict
        d['plot'] = self.plot_dict
        d['val'] = self.val_dict
        self.all_dict = d

    def key2fn(self, key, kind):
        if kind=='plot':
            return 'fig_'+key+'.pdf'
        if kind=='table':
            return 'tab_'+key+'.tex'
        if kind=='val':
            return 'val_'+key+'.tex'
            
    def create_file(self, kind, name):
        i = 0
        for key, func in self.all_dict[kind].iteritems():
            if kind=='plot':
                if name=='all':
                    func()
                elif key.count(name)==1:
                    func()
                else:
                    continue
                    
                fn = self.key2fn(key, 'plot')
                plt.gcf().savefig(fn)

            elif kind=='table':
                if name=='all':
                    lines = func()
                elif key.count(name)==1:
                    lines = func()
                else:
                    continue
                    
                # Remove last \\
                fn = self.key2fn(key, 'table')
                with open(fn,'w') as f:
                    f.writelines("\n".join(lines))

            elif kind=='val':
                fn = self.key2fn(key, 'val')
                if name=='all':
                    lines = func()
                elif name==key:
                    lines = func()
                else:
                    continue

                lines1 = [
                    "\\newcommand{\%s}[1]{%%" % key,
                    "\IfEqCase{#1}{",
                ]

                lines2 = [
                    "}[XX]",
                    "}"
                ]
                lines = lines1 + lines + lines2

                with open(fn,'w') as f:
                    f.writelines("%\n".join(lines))

            i+=1

        if i==0:
            assert False, name + " not a valid key"

    def update_paper(self):
        for kind, d in self.all_dict.iteritems():
            for key, val in d.iteritems():
                fn = self.key2fn(key, kind)
                cmd = 'cp {} paper/'.format(fn)
                print cmd
                os.system(cmd)

if __name__=="__main__":
    main()

