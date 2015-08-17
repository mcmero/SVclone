from setuptools import setup

setup(name='SVClone',
      version='0.1.1',
      description='Cluster structural variants of subclonal origin',
      url='https://bitbucket.org/mcmero/svclone',
      author='Marek Cmero',
      author_email='mcmero@student.unimelb.edu.au',
      license='Unimelb',
      packages=['SVClone'],
      entry_points = {
        'console_scripts': ['svclone = svclone.cmd:main'],
      },
      long_description=open('README.md').read(),
      zip_safe=False)

