# modelTissueFlow

This is a program for analyzing time series movies for computing tissue velocity (PIV) and other signals (myosin)


## Installing

Use PyPI: <https://pypi.org/project/modelTissueFlow/0.1/>:
pip install modelTissueFlow==0.1

### To build from source

Download the package from the Github: https://github.com/HiBandan/modelTissueFlow/archive/refs/heads/main.zip
or clone using git

    git clone https://github.com/HiBandan/modelTissueFlow.git

Using distutils create a local (in the same directory) compilation of the Cython files:

    python setup.py build_ext --inplace

Or for the global installation, use:

    python setup.py install 
