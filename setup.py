import sys
from Cython.Build import cythonize

from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension
  
#heatsim2=numpy_Extension('heatsim',sources=['heatsim.c'],extra_compile_args=['-fopenmp'],extra_link_args=['-lgomp'])


ext_modules=cythonize("heatsim2/*.pyx")

emdict=dict([ (module.name,module) for module in ext_modules])

adi_c_pyx_ext=emdict['heatsim2.alternatingdirection_c_pyx']
adi_c_pyx_ext.sources.append("heatsim2/alternatingdirection_c.c")

#ext_modules.append(numpy_Extension("heatsim2.alternatingdirection_c",["heatsim2/alternatingdirection_c.c"]))

numpy_setup(name="heatsim2",
            description="heatsim2",
            author="Stephen D. Holland",
            url="http://thermal.cnde.iastate.edu",
            ext_modules=ext_modules,
            packages=["heatsim2"])
