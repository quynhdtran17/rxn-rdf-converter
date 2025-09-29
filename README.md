# Reaction Knowledge Graph Processor: rxn_rdf_converter

## Overview

The ``rxn_rdf_converter`` is a Python package which aims to process reaction data stored in Google Protocol Buffers with the Open Reaction Database (ORD) schema into **a knowledge graph representation** as **Resource Description Framework (RDF) triples** using **MDS-Onto** (a domain ontology for Materials Data Science) as the semantic model.

This package facilitates the transformation of raw experimental data into structured RDF triples (Turtle or JSON-LD format), making the data semantically searchable, linkable, and machine-readable for advanced data analysis and machine learning applications in chemistry and materials science.

---
## Motivation

We built this package to streamline the data integration of reaction/synthesis data into one centralized database with formulation, manufacturing, and degradation.

The bottleneck of mapping reaction data to an ontology can be tedious and error-prone. Thus, this package hopefully will provide an automated tool to reduce the time needed to integrate data. 

---
## Installation

### Prerequisites

* Python 3.11+
* The `ord-schema` library.
* The `owlready2` library for ontology handling.
* The `rdflib` library for RDF graph generation.
* **RDKit** for chemical identifier normalization (InChIKey, SMILES conversion).

### Setup
Install the required dependencies using pip 

---
## Package Usage

The package will convert a dataset (that contains hundreds to thousands of reactions) in ORD schema in Google Protocol Buffers format into Resource Description Framework (RDF) triples in JSON-LD or Turtle format using MDS-Onto as the semantic model. 

The core workflow involves initializing a DatasetProcessor to manage logging and file paths, and then iterating over its reactions using the ReactionKG class to build the individual knowledge graphs (Turtle or JSON-LD serialization format).

The package is capable of batch processing multiple datasets or only one dataset. In addition, using the ReactionKG class, a user can generate one single reaction at a time. 


---
## Limitation

This package currently only works for reaction data in ORD schema and MDS-Onto is passed as an argument. It will not work with data in other Google Protocol Buffer schemas or other ontologies. 

If a user has more than 50 datasets (each contains hundreds to thousands of reactions), running the package to process multiple datasets will cause a crash since there is not enough in-memory storage. The multiple-dataset batch processing was designed to run on distributed, parallel computing infrastructure with Hadoop ecosystem. 

---
## Affiliations: 
Materials Data Science for Stockpile Stewardship Center of Excellence (MDS3-COE),
Solar Durability and Lifetime Extension (SDLE) Research Center, 
Materials Science and Engineering,
Case Western Reserve University,
Cleveland, OH 44106, USA

---
## Python package documentation
https://rxn-rdf-converter.readthedocs.io/en/latest/

---

## Acknowledgements: 

This work was supported by the U.S. Department of Energy’s Office of Energy Efficiency and Renewable Energy (EERE) under Solar Energy Technologies Office (SETO) Agreement Numbers DE-EE0009353 and DE-EE0009347, Department of Energy (National Nuclear Security Administration) under Award Number DE-NA0004104 and Contract number B647887, and U.S. National Science Foundation Award under Award Number 2133576.