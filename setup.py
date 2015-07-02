from setuptools import setup

setup(name='SVClone',
      version='0.1',
      description='Cluster structural variants of subclonal origin',
      url='https://bitbucket.org/mcmero/SVClone',
      author='Marek Cmero',
      author_email='cmerom@gmail.com',
      license='Unimelb',
      packages=['SVClone'],
      entry_points = {
        'console_scripts': ['build_sets=SVClone.cmd:main'],
      },
      zip_safe=False)


