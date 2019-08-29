.. Cocktail Shaker documentation master file, created by
   sphinx-quickstart on Wed Aug 28 23:41:51 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Cocktail Shaker's documentation!
===========================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. Cocktail Shaker documentation master file, created by sphinx-quickstart on Tue Mar 24 16:12:38 2015.

Cocktail Shaker
===============

.. sectionauthor:: Suliman Sharif <sharifsuliman1@gmail.com>

Cocktail Shaker is a **high-performance drug enumeration and expansion library**.
Cocktail Shaker leverages the computational power of  **RDKit** to create and enumerate large volumes of drug compounds.

    >>> cocktail = Cocktail('c1cc(CCCBr)ccc1')
    >>> new_compounds = cocktail.shake()
    >>> FileWriter('example', new_compounds, 'mol2')

Cocktail Shaker makes your drug enumeration and expansion life easy. It also generates your files for you in as many
formats needed for any cheminformatics software.

Features
--------


- File parsing of TXT, SDF, and Chemical Smiles.
- File writing in a variety of formats some of which include: cif, sdf, pdb, mol, mol2 and many others
- Ability to recognize and expand libraries of compounds some of which include halogens, acyl halides, aldehydes.
- Ability to enumerate in 1D, and 2D structures and produce those compounds.
- Supports Python versions 3.3+.
- Released under the `MIT license`_.

User guide
----------

A step-by-step guide to getting started with Cocktail Shaker.

.. toctree::
    :maxdepth: 2

    guide/installation
    guide/quickstart
    guide/functional_groups
    guide/contributing

API documentation
-----------------

Comprehensive API documentation with information on every function, class and method.

.. toctree::
    :maxdepth: 2
    guide/cocktail_shaker.rst

.. _`MIT license`: https://github.com/Sulstice/Cocktail-Shaker/blob/master/LICENSE