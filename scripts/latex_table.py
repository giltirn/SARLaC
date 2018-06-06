import math
import numpy
import os
import sys
import re

def hdf5_print(filename, **kwargs):
    arg_str = ""
    for a in kwargs.keys():
        arg_str = arg_str+(" -%s %s" % (a,kwargs[a]))

    out = os.popen("hdf5_print %s %s" % (filename,arg_str)).read()
    out = out.rstrip('\r\n')
    return out


def latex_print(filename, elem, **kwargs):
    if('format' not in kwargs):
        kwargs['format'] = 'publication'

    kwargs['elem'] = elem
    kwargs['noindex'] = ""

    out = hdf5_print(filename,**kwargs)
    out = "$%s$" % out
    return out


def tabulate(cols,args,**kwargs):
    print r'\FloatBarrier'
    print r'\begin{table}[h]'
    print r'\centering'

    colstr=""
    if('colstr' in kwargs.keys()):
        colstr = kwargs['colstr']
    else:
        for i in range(cols):
            colstr=colstr+'c'
    
    print "\\begin{tabular}{%s}" % colstr
    print r'\hline\hline'

    np = len(args)
    ii = 0
    for i in range(np):
        if(args[i] == "\\hline"):
            if(ii > 0):
                sys.stdout.write("\\\\\n")
            ii=0
            sys.stdout.write("\\hline\n")
            continue

        if(ii > 0):
            sys.stdout.write(' & ')

        n = 1
        mcol = None
        if(isinstance(args[i],basestring)):
            mcol = re.search(r'multicolumn\{(\d+)\}',args[i])

        if(mcol):
             n = int(mcol.group(1))
             sys.stdout.write(args[i])
        elif(isinstance(args[i],basestring)):
            sys.stdout.write(args[i])
        elif(isinstance(args[i],int)):
            sys.stdout.write("%d" % args[i])
        elif(isinstance(args[i],list)):
            sys.stdout.write(latex_print(args[i][0],args[i][1],**args[i][2]))
        else:
            sys.exit("latex_table::tabulate Unknown format %s" % type(args[i]))

        ii = ii+n

        if(ii == cols):
            sys.stdout.write("\\\\\n")
            ii=0
    
    print "\\end{tabular}"
    print "\\end{table}"
    print "\\FloatBarrier"

