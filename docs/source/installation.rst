Installation
============

Requirements
------------

* Python >= 3.7
* protobuf>=3.20
* rdkit>=2022.9.5
* ord_schema>=0.3.99
* owlready2 >= 0.45
* rdflib>=6.2.0

Install from PyPI
-----------------

.. code-block:: bash

   pip install rxn_rdf_converter

Install from source
-------------------

.. code-block:: bash

   git clone https://github.com/quynhdtran17/rxn-rdf-converter 
   cd rxn_rdf_converter
   pip install -e .

Development Installation
------------------------

.. code-block:: bash

   git clone https://github.com/quynhdtran17/rxn-rdf-converter.git
   cd rxn_rdf_converter
   pip install -e ".[dev,docs]"