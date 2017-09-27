"""Setup RapidMoc."""

from setuptools import setup

setup(name='rapidmoc',
      description='Calculates observational-style decomposition of AMOC using '
                  'output from an ocean general circulation model.',
      packages=['rapidmoc'],
      package_dir={'rapidmoc': 'rapidmoc'},
      package_data={'rapidmoc': ['etc/*']},
      install_requires=['setuptools', 'netCDF4', 'numpy',  'matplotlib',
                        'scipy'],
      entry_points={
          'console_scripts':
          ['run_rapidmoc.py = rapidmoc.rapidmoc:main']},
      zip_safe=False)
