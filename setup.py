import sys
from Cython.Build import cythonize

from numpy.distutils.core import setup as numpy_setup, Extension as numpy_Extension
  
#heatsim2=numpy_Extension('heatsim',sources=['heatsim.c'],extra_compile_args=['-fopenmp'],extra_link_args=['-lgomp'])


class install_lib_save_version(install_lib):
    """Save version information"""
    def run(self):
        install_lib.run(self)
        
        for package in self.distribution.command_obj["build_py"].packages:
            install_dir=os.path.join(*([self.install_dir] + package.split('.')))
            fh=open(os.path.join(install_dir,"version.txt"),"w")
            fh.write("%s\n" % (version))  # version global, as created below
            fh.close()
            pass
        pass
    pass



# Extract GIT version
if os.path.exists(".git") and distutils.spawn.find_executable("git") is not None:
    # Check if tree has been modified
    modified = subprocess.call(["git","diff-index","--quiet","HEAD","--"]) != 0
    
    gitrev = subprocess.check_output(["git","rev-parse","HEAD"]).strip()

    version = "git-%s" % (gitrev)

    # See if we can get a more meaningful description from "git describe"
    try:
        versionraw=subprocess.check_output(["git","describe","--tags","--match=v*"],stderr=subprocess.STDOUT).decode('utf-8').strip()
        # versionraw is like v0.1.0-50-g434343
        # for compatibility with PEP 440, change it to
        # something like 0.1.0+50.g434343
        matchobj=re.match(r"""v([^.]+[.][^.]+[.][^-.]+)(-.*)?""",versionraw)
        version=matchobj.group(1)
        if matchobj.group(2) is not None:
            version += '+'+matchobj.group(2)[1:].replace("-",".")
            pass
        pass
    except subprocess.CalledProcessError:
        # Ignore error, falling back to above version string
        pass

    if modified and version.find('+') >= 0:
        version += ".modified"
        pass
    elif modified:
        version += "+modified"
        pass
    pass
else:
    version = "UNKNOWN"
    pass

print("version = %s" % (version))


ext_modules=cythonize("heatsim2/*.pyx")

emdict=dict([ (module.name,module) for module in ext_modules])

adi_c_pyx_ext=emdict['heatsim2.alternatingdirection_c_pyx']
adi_c_pyx_ext.sources.append("heatsim2/alternatingdirection_c.c")

#ext_modules.append(numpy_Extension("heatsim2.alternatingdirection_c",["heatsim2/alternatingdirection_c.c"]))

numpy_setup(name="heatsim2",
            description="heatsim2",
            author="Stephen D. Holland",
            version=version,
            zip_safe=False,
            url="http://thermal.cnde.iastate.edu",
            ext_modules=ext_modules,
            packages=["heatsim2"],
            cmdclass={"install_lib": install_lib_save_version} )
