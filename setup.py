from setuptools import setup, find_packages

setup(name='SVclone',
      version='0.1.1',
      description='Cluster structural variants of subclonal origin',
      url='https://bitbucket.org/mcmero/svclone',
      author='Marek Cmero',
      author_email='cmerom[at]gmail.com',
      license='',
      packages=['SVclone', 'SVclone.SVprocess'],
      scripts=['SVclone.py'],
      long_description=open('README.md').read(),
      zip_safe=False)

