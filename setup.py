from setuptools import setup

with open('README.rst') as file:
        long_description = file.read()

setup(
    name='SuperBoL',
    version='0.4.0',
    description='Supernova Bolometric Lightcurves',
    long_description = long_description,
    author='Jeremy A. Lusk',
    author_email='jeremy.lusk@gmail.com',
    url='https://github.com/JALusk/SuperBoL',
    packages=['superbol'],
    #package_data={'superbol' : ['data/sn_data.h5']},
    license='MIT License',
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Astronomy",
        ],
    install_requires=[
        'astropy',
        'numpy',
        'pyerfa',
        'scipy'
    ],
    )   
