Installation
============

Requirements
------------

* Python >= 3.7
* protobuf
* rdkit
* ord_schema
* owlready2
* rdflib

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