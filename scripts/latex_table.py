import math
import numpy
import os
import sys
import re

#--------------------------------------------------------------------------------------------------------------------------------------
# Automatically generate latex tables calling out to hdf5_print to print parameter entries
# Usage example: 
# 1) Import script:   from latex_table import *      (Make sure the directory of this script is in your PYTHONPATH)
# 2) Generate a 'format array' with table entries:
#     Basic options include strings such as latex commands, and integers
#     Special latex commands '\hline' for a horizontal line, '\multicolumn{ncols}{colfmt}{entry}' for multi-column entries are treated appropriately
#     Interface with hdf5_print is performed using an entry in the form of an array with the first entry the filename, the second the element index 
#     (for multi-dimensional arrays use a string with comma-separated indices), and finally a hash with any extra commands passed to the program (use {} for default)
#     Extra commands are the usual command line args but without the leading -, eg "pub_rounding"
# 3) Call 'tabulate' with the first argument the number of columns, the second the format array and finally a string with the column layout

# Example:

# from latex_table import * 

# fmt = [ '$M_0$',["params.hdf5", 6, {"pub_rounding":-2}],
#         '$M_1$',["params.hdf5", 7, {"pub_rounding":-2}],
#         '$\chi^2$',["chisq.hdf5",0,{}],
#         r'$\chi^2/{\rm dof}$',["chisq_per_dof.hdf5",0,{}]
# ]

# tabulate(2,fmt,colstr=r'c|c')
#--------------------------------------------------------------------------------------------------------------------------------------

def hdf5_print(filename, **kwargs):
    if(os.path.isfile(filename) == False):
        return "ERR"

    debug = False;
    if('debug' in kwargs.keys()):
        debug = True
        del kwargs['debug']

    arg_str = ""
    for a in kwargs.keys():
        arg_str = arg_str+(" -%s %s" % (a,kwargs[a]))

    if(debug == True):
        print >> sys.stderr, "DEBUG hdf5_print %s %s" % (filename, arg_str)

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

#Print the nth line stored in a file
def line_print(filename, line=0):
    if(os.path.isfile(filename) == False):
        return "ERR %s" % filename
    file = open(filename, 'r') 
    for i in range(line+1):
        out = file.readline()
    out = out.rstrip('\r\n')
    return out

def tabulate(cols,args,**kwargs):
    longtable = False
    if('longtable' in kwargs.keys() and kwargs['longtable'] == True):
        longtable = True

    colstr=""
    if('colstr' in kwargs.keys()):
        colstr = kwargs['colstr']
    else:
        for i in range(cols):
            colstr=colstr+'c'


    print r'\FloatBarrier'

    if longtable == True:
        print "\\begin{longtable}{%s}" % colstr

        if('shrink_to_fit' in kwargs.keys() and kwargs['shrink_to_fit'] == True):
            print "shrink_to_fit not supported for longtable"
            sys.exit(1)

    else:
        print r'\begin{table}[h]'
        print r'\centering'

        if('shrink_to_fit' in kwargs.keys() and kwargs['shrink_to_fit'] == True):
            print r'\resizebox{\textwidth}{!}{'
    
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
        elif(isinstance(args[i],float)):
            sys.stdout.write("%f" % args[i])
        elif(isinstance(args[i],list)):
            sys.stdout.write(latex_print(args[i][0],args[i][1],**args[i][2]))
        else:
            sys.exit("latex_table::tabulate Unknown format %s" % type(args[i]))

        ii = ii+n

        if(ii == cols):
            sys.stdout.write("\\\\\n")
            ii=0
    
    if longtable == False:
        print "\\end{tabular}"

        if('shrink_to_fit' in kwargs.keys() and kwargs['shrink_to_fit'] == True):
            print r'}'

    if('caption' in kwargs.keys() ):
        print "\\caption{%s}" % kwargs['caption']

    if longtable == True:
        print "\\end{longtable}"
    else:
        print "\\end{table}"

    print "\\FloatBarrier"



def hdf5_calc(outfile, expr, symbols, filenames, index_strings, **kwargs):
    for f in filenames:
        if(os.path.isfile(f) == False):
            return "ERR"

    nsymb = len(symbols)
    if(len(filenames) != nsymb or len(index_strings) != nsymb):
        print "hdf5_calc Input arrays must be same size"
        exit;

    arg_str = "\"%s\"" % expr
    for i in range(nsymb):
        arg_str = arg_str + " \"%s\" %s \"%s\"" % (symbols[i], filenames[i], index_strings[i])

    os.popen("hdf5_calc %s %s" % (outfile,arg_str)).read()
