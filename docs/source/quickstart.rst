Quick Start Guide 
=================

This guide will get you up and running with rxn_rdf_converter in minutes 

Installation
------------

.. code-block:: bash

    pip install rxn_rdf_converter

Basic Usage
------------

Process Multiple Datasets CLI
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python rxn_rdf_converter all --dataset_root ../datasets/ --save_path ../save_path --onto_file_path ../onto_file_path --error_log_directory ../error_log_directory


Process Single Datasets CLI
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    python rxn_rdf_converter single-dataset dataset_file_path ../dataset_file_path  --save_path ../save_path --onto_file_path ../onto_file_path --error_log_directory ../error_log_directory


Process Multiple Datasets
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import ord_schema
    from ord_schema.message_helpers import load_message, write_message, message_to_row
    from ord_schema.proto import dataset_pb2, reaction_pb2
    import os
    from rdkit import Chem
    import re
    from owlready2 import get_namespace, get_ontology, Thing
    import rdflib
    from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
    from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
    import logging
    import csv
    import rxn_rdf_converter
    from rxn_rdf_converter import DatasetProcessor

    # Create logs directory if it doesn't exist
    error_log_directory = '../error_logs'
    os.makedirs(error_log_directory, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(error_log_directory, 'main.log')),
            logging.StreamHandler()
        ]
    )

    logger = logging.getLogger(__name__)

    def main(): 
        """ Main Executive Function """
        try: 
            logging.info("Starting Data Processing....")

            # set up path: 
            save_path = '../save_path'
            onto_file_path = '../MDS-Onto.owl'

            #logger.info(f"Found {len(file_list)} data files")

            dataset_reaction_list = []
            
            for dataset in input_df['file_list']:
                try: 
                    logger.info(f"Processing dataset {dataset}")

                    dataset_processor = rxn_rdf_converter.DatasetProcessor(
                        dataset_pb=dataset_pb2,
                        dataset_file_path=dataset,
                        owl_onto_file_path=onto_file_path,
                        output_directory=save_path,
                        error_log_directory=error_log_directory,
                        fmt='json-ld'
                    )

                    _, reaction_error, dataset_reaction_list = dataset_processor.extract_reaction(dataset_reaction_list)

                    logger.info(f"Successfully completed dataset {dataset}")
                except Exception as e: 
                    logger.error(f"Failed to process dataset {dataset} - Error: {e}")

                finally:
                    # Clean up logger resources
                    if 'dataset_process' in locals():
                        dataset_processor.cleanup_logger()

            csv_output_path = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/output_logs'
            os.makedirs(csv_output_path, exist_ok=True)
            
                # save the results
            with open(os.path.join(csv_output_path, 'dataset_reactions.csv'), 'w', newline='') as f: 
                writer = csv.writer(f)
                writer.writerow(['dataset_id', 'reaction_id'])
                writer.writerows(dataset_reaction_list)
            
            logger.info(f"Data processing completed successfully")

        except Exception as e: 
            logger.error(f"Error in main execution: {e}", exc_info=True)
            

    if __name__ == '__main__': 
        main()
    

Process Individual Dataset
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    import ord_schema
    from ord_schema.message_helpers import load_message, write_message, message_to_row
    from ord_schema.proto import dataset_pb2, reaction_pb2
    import os
    from rdkit import Chem
    import re
    from owlready2 import get_namespace, get_ontology, Thing
    import rdflib
    from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
    from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
    import logging
    import csv
    import rxn_rdf_converter
    from rxn_rdf_converter import DatasetProcessor

    # =================================================================
    #               SETUP LOGGING ERRORS
    # =================================================================
    logging.basicConfig(
        level=logging.INFO, 
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('error.log'),
            logging.StreamHandler()
        ]
    )

    logger = logging.getLogger(__name__)

    # =================================================================
    #               INITIATE FILE PATH
    # =================================================================
    
    def setup_file_path(dataset_path):
        """ Set up the file paths """

        file_list = []
        for root, dirs, files in os.walk(dataset_path):
            for name in files: 
                if name.startswith('ord_dataset'):
                    file_path = os.path.join(root, name)
                    file_list.append(file_path)
        
        return file_list

    # =================================================================
    #               MAIN FUNCTION
    # =================================================================

    dataset_path = '../datasets/'
    save_path = '../save_path'
    onto_file_path = '../MDS-Onto.owl''
    error_log_directory = '../error_log_directory'

    try: 
        logging.info("Starting Data Processing....")

        # set up path: 
        file_list,  = setup_file_path(dataset_path)
        logger.info(f"Found {len(file_list)} data files")

        dataset_file_path = file_list[] # add in here the index of the dataset of interest

        dataset_reaction_list = []

        dataset_1 = rxn_rdf_converter.DatasetProcessor(
            dataset_pb=dataset_pb2,
            dataset_file_path=dataset_file_path,
            owl_onto_file_path=mds_file_path,
            output_directory=save_path,
            fmt='json-ld'
        )

        _, reaction_error, dataset_reaction_list = dataset_1.extract_reaction(dataset_reaction_list)
            # save the results
        with open('dataset_reactions.csv', 'w', newline='') as f: 
            writer = csv.writer(f)
            writer.writerow(['dataset_id', 'reaction_id'])
            writer.writerows(dataset_reaction_list)

        print(f"Collected {len(dataset_reaction_list)} reaction mappings")
        print("Saved to dataset_reactions.csv")
    except Exception as e: 
        logger.error(f"Error in main execution: {e}", exc_info=True)


Process Individual Reaction
~~~~~~~~~~~~~~

.. code-block:: python

    import ord_schema
    from ord_schema.message_helpers import load_message, write_message, message_to_row
    from ord_schema.proto import dataset_pb2, reaction_pb2
    import os
    from rdkit import Chem
    import re
    from owlready2 import get_namespace, get_ontology, Thing
    import rdflib
    from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
    from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
    import logging
    import csv
    import rxn_rdf_converter
    from rxn_rdf_converter import ReactionKG

    def setup_file_path(dataset_path):
        """ Set up the file paths """

        file_list = []
        for root, dirs, files in os.walk(dataset_path):
            for name in files: 
                if name.startswith('ord_dataset'):
                    file_path = os.path.join(root, name)
                    file_list.append(file_path)
        
        return file_list
    dataset_path = '../datasets/'
    save_path = '../save_path'
    onto_file_path = '../MDS-Onto.owl''

    # Process all of the file paths of all of the datasets:
    file_list = setup_file_path(dataset_path) 

    # Load one dataset into a Python variable by calling the dataset index into file_list[], all reactions will be generated in a dataset from ORD 
    dataset = load_message(file_list[], dataset_pb2.Dataset,)

    # Process the reaction of interest by adding the index of the reaction in dataset.reactions[]
    reaction_1 = rxn_rdf_converter.ReactionKG(dataset.reactions[], fmt="json-ld").generate_reaction().generate_instances(onto_file_path).generate_data_graph(dataset.dataset_id, save_path)


Next Steps
----------

* Read the full :doc: 'rxn_rdf_converter' documentation
* Browse the :doc: 'modules' for detailed API reference
* Check out more examples in the main documentation