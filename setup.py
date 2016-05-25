from setuptools import setup

setup(
    name='SNoBoL',
    version='0.2.1',
    description='Supernova Bolometric Lightcurves',
    author='Jeremy A. Lusk',
    author_email='jeremy.lusk@gmail.com',
    url='https://github.com/JALusk/SNoBoL',
    packages=['snobol'],
    package_data={'snobol' : ['data/sn_data.h5']},
    license='MIT License',
    install_requires=[
        "astropy",
        "scipy>=0.14",
        "numpy",
        "specutils"
    ],
    )   
