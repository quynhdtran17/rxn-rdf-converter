# Project Description for rxn-rdf-converter

## Introduction


``rxn_rdf_converter`` is a Python package which aims to process reaction data stored in Google Protocol Buffers with the Open Reaction Database (ORD) schema into Resource Description Framework (RDF) triples using MDS-Onto (a domain ontology for Materials Data Science) as the semantic model.

The package is meant to be a user-friendly tool for easily converting reaction data into RDF triples by providing a semantic model that is connected to vocabulary and terms beyond Chemistry. 

MDS-Onto enables reaction data to be linked and integrated with data from sources in Manufacturing, Formulation, Degradation, Biomedical, etc. 

The package is capable of batch processing datasets or only one dataset. 


## Motivation

We built this package to streamline the data integration of reaction/synthesis data into one centralized database with formulation, manufacturing, and degradation.

The bottleneck of mapping reaction data to an ontology can be tedious and error-prone. Thus, this package hopefully will provide an automated tool to reduce the time needed to integrate data. 

## Limitation

This package currently only works for reaction data in ORD schema and MDS-Onto is passed as an argument. It will not work with data in other Google Protocol Buffer schemas or other ontologies. 

## Affiliations: 
Materials Data Science for Stockpile Stewardship Center of Excellence (MDS3-COE),
Solar Durability and Lifetime Extension (SDLE) Research Center, 
Materials Science and Engineering,
Case Western Reserve University,
Cleveland, OH 44106, USA


## Package Usage: 

The package will convert a dataset (that contains hundreds to thousands of reactions) in ORD schema in Google Protocol Buffers format into Resource Description Framework (RDF) triples in JSON-LD or Turtle format using MDS-Onto as the semantic model. 

## Python package documentation
https://rxn-rdf-converter.readthedocs.io/en/latest/

## Acknowledgements: 

This work was supported by the U.S. Department of Energy’s Office of Energy Efficiency and Renewable Energy (EERE) under Solar Energy Technologies Office (SETO) Agreement Numbers DE-EE0009353 and DE-EE0009347, Department of Energy (National Nuclear Security Administration) under Award Number DE-NA0004104 and Contract number B647887, and U.S. National Science Foundation Award under Award Number 2133576.