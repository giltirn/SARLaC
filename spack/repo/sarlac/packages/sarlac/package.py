# Copyright 2013-2020 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack import *


class Sarlac(AutotoolsPackage):
    """The Statistical Analysis and Regression for Lattice Calculations (SARLaC) package"""

    homepage = "https://github.com/giltirn/SARLaC"
    git = "https://github.com/giltirn/SARLaC"

    version('main', branch='main', preferred=True)
    version('ckelly_develop', branch='ckelly_develop')

    variant('python', default=False, description='Build with online Python2 interface')
    variant('minuit', default=False, description='Build with Minuit fit wrappers')

    depends_on('python@2.7:2.9', when="+python", type=("build","link","run"))
    depends_on('minuit', when="+minuit",type=("build","link"))
    depends_on('gsl',type=("build","link"))
    depends_on('hdf5 +cxx',type=("build","link"))
    depends_on('boost @1.61: +serialization +timer +iostreams +filesystem +system', type=("build","link"))
    depends_on('boost @1.61: +serialization +timer +iostreams +filesystem +system +python', when="+python", type=("build","link"))

    depends_on('autoconf', type='build')
    depends_on('automake', type='build')
    depends_on('libtool',  type='build')
    depends_on('m4',       type='build')

    def configure_args(self):
        args = ["--with-gsl=%s" % self.spec['gsl'].prefix, "--with-hdf5=%s" % self.spec['hdf5'].prefix]
        if '+minuit' in self.spec:
            args.append("--with-minuit2=%s" % self.spec['minuit'].prefix)
        if '+python' in self.spec:
            args.append("--with-python=yes")
               
        return args
