from setuptools import setup

setup(name='SVClone',
      version='0.1.1',
      description='Cluster structural variants of subclonal origin',
      url='https://bitbucket.org/mcmero/svclone',
      author='Marek Cmero',
      author_email='mcmero@student.unimelb.edu.au',
      license='',
      packages=['SVClone'],
      scripts=['SVClone.py'],
      long_description=open('README.md').read(),
      zip_safe=False)

