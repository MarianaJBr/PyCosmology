'''
Created on Jun 11, 2012

@author: Steven
'''

from numpy.distutils.core import setup,Extension


version_info = [0,2,5]
version = '.'.join([str(i) for i in version_info])
def generate_version_py():
    fid = open("__version.py",'w')
    try:
        fid.write("version_info = %s\n" % version_info)
        fid.write("version = %s\n" % version)
    finally:
        fid.close()
        
generate_version_py()

shapelets = Extension('PyCosmology.shapelets.fort.shapelets',['PyCosmology/shapelets/fort/find_coeffs.f90'],extra_f90_compile_args=['-Wtabs'],f2py_options=['--quiet'])
read_sim = Extension('PyCosmology.sims.fort.read_sim',['PyCosmology/sims/fort/ReadSim.f90'],extra_f90_compile_args=['-Wtabs'],f2py_options=['--quiet'])

if __name__=="__main__":
    
    setup(name = 'PyCosmology',
          version = version,
          author = 'Steven Murray',
          author_email = 'steven.jeanette.m@gmail.com',
          description = 'Bits-and-bobs cosmology project.',
          url = 'doesnt.have.one.yet.com',
          ext_modules = [shapelets,read_sim],
          packages = ['PyCosmology','PyCosmology.shapelets','PyCosmology.sims','PyCosmology.structure','PyCosmology.shapelets.fort','PyCosmology.sims.fort']
    )