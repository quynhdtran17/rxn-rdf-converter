Introduction
============

``rxn_rdf_converter`` is a Python package which aims to process reaction data stored in Google Protocol Buffers with the Open Reaction Database (ORD) schema into Resource Description Framework (RDF) triples using MDS-Onto (a domain ontology for Materials Data Science) as the semantic model.

The package is meant to be a user-friendly tool for easily converting reaction data into RDF triples by providing a semantic model that is connected to vocabulary and terms beyond Chemistry. 

MDS-Onto enables reaction data to be linked and integrated with data from sources in Manufacturing, Formulation, Degradation, Biomedical, etc. 

The package is capable of batch processing datasets or only one dataset. 

Motivation
==========

We built this package to streamline the data integration of reaction/synthesis data into one centralized database with formulation, manufacturing, and degradation.

The bottleneck of mapping reaction data to an ontology can be tedious and error-prone. Thus, this package hopefully will provide an automated tool to reduce the time needed to integrate data. 


Limitation
==========

This package currently only works for reaction data in ORD schema and MDS-Onto is passed as an argument. It will not work with data in other Google Protocol Buffer schemas or other ontologies. 