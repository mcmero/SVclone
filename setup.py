from setuptools import setup

setup(name='SVClone',
      version='0.1.1',
      description='Cluster structural variants of subclonal origin',
      url='https://bitbucket.org/mcmero/svclone',
      author='Marek Cmero',
      author_email='cmerom@gmail.com',
      license='Unimelb',
      packages=['SVClone'],
      entry_points = {
        'console_scripts': ['build_sets=svclone.cmd:main'],
      },
      zip_safe=False)

