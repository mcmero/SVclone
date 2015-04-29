from setuptools import setup

setup(name='phylo_sv',
      version='0.1',
      description='Obtain post-processing data from BAM files based on SV calls and build corresponding phylogenetic trees',
      url='https://bitbucket.org/mcmero/sv_proc',
      author='Marek Cmero',
      author_email='cmerom@gmail.com',
      license='Unimelb',
      packages=['phylo_sv'],
      test_suite='nose.collector',
      tests_require=['numpy','scipy','nose'],
      entry_points = {
        'console_scripts': ['build_sets=phylo_sv.cmd:main'],
      },
      zip_safe=False)

