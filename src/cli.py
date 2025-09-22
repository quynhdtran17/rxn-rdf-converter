# =================================================================
#               IMPORT REQUIREMENTS
# =================================================================
from ord_schema.message_helpers import load_message, write_message, message_to_row
from ord_schema.proto import dataset_pb2, reaction_pb2
import os
from rdkit import Chem
import re
from owlready2 import get_namespace, get_ontology, Thing
import rdflib
from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
from datetime import datetime
import rxn_rdf_converter
import logging
import csv

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
dataset_path = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/OpenReactionDB/ord-data/data' # path of dataset 
savepath = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/Extracted_ORD_Data'
onto_file_path = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/chemOntologies/mdsChemRxn(v.0.3.0.7).owl'
error_log_directory = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/error_logs'

def setup_file_path(dataset_path, savepath, onto_file_path, error_log_directory):
    """ Set up the file paths """
    
    file_list = []
    for root, dirs, files in os.walk(dataset_path):
        for name in files: 
            if name.startswith('ord_dataset'):
                file_path = os.path.join(root, name)
                file_list.append(file_path)
    
    return dataset_path, savepath, onto_file_path, error_log_directory, file_list

# =================================================================
#               MAIN FUNCTION
# =================================================================

def main(): 
    """ Main Executive Function """
    try: 
        logging.info("Starting Data Processing....")

        # set up path: 
        dataset_path, savepath, onto_file_path, error_log_directory, file_list = setup_file_path(dataset_path, savepath, onto_file_path, error_log_directory)
        logger.info(f"Found {len(file_list)} data files")

        dataset_reaction_list = []

        dataset_1 = rxn_rdf_converter.DatasetProcessor(
            dataset_pb=dataset_pb2,
            dataset_file_path=file_list[2],
            owl_onto_file_path=onto_file_path,
            output_directory=savepath,
            error_log_directory=error_log_directory,
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
        
if __name__ == '__main__': 
    main()

# Extract one reaction into one JSON-LD 
#dataset = load_message(file_list[4], dataset_pb2.Dataset,)
#reaction_1 = rxn_rdf_converter.ReactionKG(dataset.reactions[5], fmt="json-ld").generate_reaction().generate_instances(onto_file_path).generate_data_graph(dataset.dataset_id, savepath)
