from setuptools import setup, find_packages

setup(name='SVclone',
      version='v1.1.0',
      description='A computational method for inferring the cancer cell fraction of tumour structural variation from whole-genome sequencing data.',
      url='https://github.com/mcmero/SVclone',
      author='Marek Cmero',
      author_email='cmerom[at]gmail.com',
      license='BSD-3-Clause',
      packages=find_packages(),
      package_data={
          'SVclone': ['cluster_with_ccube.R', 'write_output.R', 'post_assign.R'],
      },
      entry_points={
          'console_scripts': [
              'svclone = SVclone.cli:main'
          ]
      },
      long_description=open('README.md').read(),
      zip_safe=False)

