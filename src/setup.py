'''
Created on Jun 11, 2012

@author: Steven
'''

from numpy.distutils.core import setup,Extension
import sys

version_info = [0,4,0]
version = '.'.join([str(i) for i in version_info])
def generate_version_py():
    fid = open("__version.py",'w')
    try:
        fid.write("version_info = %s\n" % version_info)
        fid.write("version = %s\n" % version)
    finally:
        fid.close()
        
generate_version_py()

healpix_library_path = "/Users/Steven/Downloads/Healpix_2.20a"

shapelets = Extension('PyCosmology.shapelets.fort.shapelets',['PyCosmology/shapelets/fort/find_coeffs.f90'],
                      extra_f90_compile_args=['-Wtabs', '-O0'],
                      f2py_options=['--quiet'])
read_sim = Extension('PyCosmology.sims.fort.read_sim',['PyCosmology/sims/fort/ReadSim.f90'],
                     depends = ['PyCosmology/sims/fort/Auxiliary.f90'],
                     extra_f90_compile_args=['-Wtabs'],
                     f2py_options=['--quiet',"skip:","recentre","inertia",'eigorder','genrotate','logbins','jacobi',':'])
groupstructure = Extension('PyCosmology.sims.fort.groupstructure',['PyCosmology/sims/fort/GroupStructureOMP.f90'],
                     extra_f90_compile_args=['-Wtabs','-fopenmp'],
                     f2py_options=['--quiet'],
                     library_dirs = [healpix_library_path+'/lib'],
                     libraries = ["healpix","gomp"],
                     include_dirs = [healpix_library_path+'/include'])
if __name__=="__main__":
    
    setup(name = 'PyCosmology',
          version = version,
          author = 'Steven Murray',
          author_email = 'steven.jeanette.m@gmail.com',
          description = 'Bits-and-bobs cosmology project.',
          url = 'doesnt.have.one.yet.com',
          ext_modules = [shapelets,read_sim,groupstructure],
          packages = ['PyCosmology','PyCosmology.shapelets','PyCosmology.sims','PyCosmology.structure','PyCosmology.shapelets.fort','PyCosmology.sims.fort'],
          data_files = [("PyCosmology/data",["PyCosmology/data/healpix.dat"])]
    )