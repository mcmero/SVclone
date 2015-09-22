from setuptools import setup

setup(name='sv_process',
      version='1.0',
      description='Obtain read count data around breakpoints from BAM files',
      url='https://bitbucket.org/mcmero/sv_preprocess',
      author='Marek Cmero',
      author_email='cmerom@gmail.com',
      license='Unimelb',
      packages=['sv_process'],
      entry_points = {
        'console_scripts': ['build_sets=sv_process.cmd:main'],
      },
      zip_safe=False)

