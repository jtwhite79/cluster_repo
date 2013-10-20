from distutils.core import setup

#To use:
# python setup.py install


setup(
    name='cluster_repo',
    version='1.37',
    packages=['cluster_repo'],
    install_requires = ['numpy>=1.7','pandas>=0.12'],
    url='https://github.com/jtwhite79/cluster_repo.git',
    license='PSF',
    author='Jeremy T. White and  Joseph D. Hughes',
    author_email='jwhite@usgs.gov or jdhughes@usgs.gov',
    description='python utilities for working with MODFLOW-based input and output data'
)
