from setuptools import setup

setup(
    name='SuperBoL',
    version='0.3.0',
    description='Supernova Bolometric Lightcurves',
    author='Jeremy A. Lusk',
    author_email='jeremy.lusk@gmail.com',
    url='https://github.com/JALusk/SuperBoL',
    packages=['superbol'],
    package_data={'superbol' : ['data/sn_data.h5']},
    license='MIT License',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "Topic :: Scientific/Engineering :: Astronomy",
        ],
    install_requires=[
        "astropy",
        "scipy>=0.14",
        "numpy",
        "specutils"
    ],
    )   
