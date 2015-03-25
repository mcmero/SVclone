from setuptools import setup

setup(name='proc_sv',
      version='0.1',
      description='Obtain post-processing data from BAM files based on SV calls',
      url='https://bitbucket.org/mcmero/sv_proc',
      author='Marek Cmero',
      author_email='cmerom@gmail.com',
      license='Unimelb',
      packages=['proc_svs'],
      test_suite='nose.collector',
      tests_require=['numpy','scipy','nose'],
      entry_points = {
        'console_scripts': ['build_sets=phyl_sv.cmd:main'],
      },
      zip_safe=False)

