from distutils.core import setup

# To use:
#	   python setup.py bdist --format=wininst

setup(name='cluster_repo',
      version='1.37',
      author="Jeremy T. White and Joseph D. Hughes",
      author_email="jwhite@usgs.gov or jdhughes@usgs.gov",
      scripts=['MFArrayUtil.py','MFBinaryClass.py','MFData.py','MFInterpolators.py','pestUtil.py','shapefile.py'],
      install_requires = ['numpy>=1.7','pandas>=0.12'],
      url = 'https://github.com/jtwhite79/cluster_repo.git',
      description = 'python utilities for working with MODFLOW-based input and output data',
      license = 'PSF',
      keywords = ['MODFLOW','MODFLOW binary reader','MODFLOW array utility','shapefile library','PEST python utility'],
     )

