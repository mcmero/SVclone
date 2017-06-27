# INSTALL #

SVclone runs on Linux and OS X. A Windows installation is currently not possible due to Windows incompatibility of the [PySam](http://pysam.readthedocs.org/en/latest/) package.

If you do not have numpy and SciPy installed, the easiest way to install the prerequisites is to install [Anaconda2](https://www.continuum.io/downloads) for Python 2.7. Then, install pymc:

    conda install pymc

(Note that pymc requires a Fortram compiler to be installed.)

Install [R](https://www.r-project.org/) to run the post-processing script. Running the script requires the following packages:

* [combinat](https://cran.r-project.org/web/packages/combinat/index.html)
* [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html)
* [gplots](https://cran.r-project.org/web/packages/gplots/index.html)
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
* [reshape](https://cran.r-project.org/web/packages/reshape/index.html)
* [gtools](https://cran.rstudio.com/web/packages/gtools/index.html)
* [plyr](https://cran.rstudio.com/web/packages/plyr/index.html)
* [circlize](https://cran.r-project.org/web/packages/circlize/index.html)

Run the following in R to install these dependencies:

    install.packages('combinat')
    install.packages('RColorBrewer')
    install.packages('gplots')
    install.packages('ggplot2')
    install.packages('gridExtra')
    install.packages('reshape')
    install.packages('gtools')
    install.packages('plyr')
    install.packages('circlize')

Then install SVclone and run the test example to check whether the install has been successful:

    curl -O https://github.com/mcmero/SVclone/archive/0.2.1.tar.gz
    tar -xvzf 0.2.1.tar.gz

    cd SVclone
    pip install -r requirements.txt
    python setup.py install

    chmod u+x run_example.sh
    ./run_example.sh

If you have admin rights and/or are running in a [python virtual environment](http://python-guide-pt-br.readthedocs.io/en/latest/dev/virtualenvs/) you should encounter no issues, otherwise you may need to add --user when installing:

    pip install -r requirements.txt --user
    python setup.py install --user

Below is a full list of Python dependencies:

Annotate/count steps:

* [Numpy](http://www.numpy.org/) - install for python 2
* [PySam](http://pysam.readthedocs.org/en/latest/)
* [PyVCF](https://pyvcf.readthedocs.org/en/latest/)

Cluster/post-assign steps:

* [Numpy](http://www.numpy.org/) - install for python 2
* [Pandas](http://pandas.pydata.org/)
* [scikit-learn](http://scikit-learn.org/stable/install.html)
* [PyMC](https://pymc-devs.github.io/pymc/INSTALL.html)
* [ipython](https://pypi.python.org/pypi/ipython)
* [ipdb](https://pypi.python.org/pypi/ipdb)
* [matplotlib](http://matplotlib.org/)
* [SciPy](https://http://www.scipy.org/)
