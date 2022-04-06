SuperBoL: Supernova Bolometric Lightcurves
==========================================

.. image:: https://github.com/JALusk/SuperBoL/actions/workflows/main.yml/badge.svg?branch=tdd_develop
   :target: https://github.com/JALusk/SuperBoL/actions/workflows/main.yml
   :alt: CI

To create a virtual environment: 
python -m venv env 

To activate env: 
.\env\Scripts\activate

To install requirements: 
pip install -r requirements.txt

To run tests in Windows: 
python -m unittest discover -s tests.unit -v
python -m unittest discover -s tests.integration -v
