# =================================================================
#               IMPORT REQUIREMENTS
# =================================================================
from ord_schema.message_helpers import load_message, message_to_row
from ord_schema.proto import dataset_pb2, reaction_pb2
import os
import re
import rdflib
from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
from datetime import datetime
import rxn_rdf_converter 
import logging
import csv
import argparse
import sys
from rdkit import Chem
from owlready2 import get_namespace, get_ontology, Thing 

# =================================================================
#               SETUP LOGGING ERRORS
# =================================================================

logging.basicConfig(
    level=logging.INFO, 
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('error.log'),
        logging.StreamHandler(sys.stdout)
    ]
)

logger = logging.getLogger(__name__)

# =================================================================
#               UTILITY FUNCTIONS
# =================================================================

def setup_file_path(dataset_path, save_path, onto_file_path, error_log_directory):
    """ Searches the dataset_path recursively for all ORD dataset files. """
    
    file_list = []
    logger.info(f"Searching for datasets in: {dataset_path}")
    if not os.path.exists(dataset_path):
        logger.error(f"Dataset root path does not exist: {dataset_path}")
        return file_list, save_path, onto_file_path, error_log_directory

    for root, dirs, files in os.walk(dataset_path):
        for name in files: 
            # ORD datasets typically start with 'ord_dataset' and are .pb files
            if name.startswith('ord_dataset') and name.endswith('.pb'):
                file_path = os.path.join(root, name)
                file_list.append(file_path)
    
    logger.info(f"Found {len(file_list)} datasets.")
    return file_list, save_path, onto_file_path, error_log_directory

def save_results(dataset_reaction_list, error_log_directory):
    """Saves the list of (dataset_id, reaction_id) pairs to a CSV file."""
    if not dataset_reaction_list:
        logger.info("No reaction data to save to CSV.")
        return
        
    csv_output_path = os.path.join(error_log_directory, 'output_logs')
    os.makedirs(csv_output_path, exist_ok=True)
    
    output_file = os.path.join(csv_output_path, 'dataset_reactions.csv')
    
    with open(output_file, 'w', newline='') as f: 
        writer = csv.writer(f)
        writer.writerow(['dataset_id', 'reaction_id'])
        writer.writerows(dataset_reaction_list)
    
    logger.info(f"Successfully saved {len(dataset_reaction_list)} reaction mappings to {output_file}")


# =================================================================
#               PROCESSING LOGIC FUNCTIONS
# =================================================================

def process_all_datasets(file_list, save_path, onto_file_path, error_log_directory):
    """Processes every dataset file path in the provided list."""
    dataset_reaction_list = []
    
    for dataset_file_path in file_list:
        try: 
            logger.info(f"Processing dataset: {dataset_file_path}")
            
            dataset_processor = rxn_rdf_converter.DatasetProcessor(
                    dataset_pb=dataset_pb2,
                    dataset_file_path=dataset,
                    owl_onto_file_path=onto_file_path,
                    output_directory=save_path,
                    error_log_directory=error_log_directory,
                    fmt='json-ld'
                )
            
            _, reaction_error, dataset_reaction_list = dataset_processor.extract_reaction(dataset_reaction_list)

            logger.info(f"Successfully completed dataset {dataset_file_path}")
        except Exception as e: 
            logger.error(f"Failed to process dataset {dataset_file_path} - Error: {e}", exc_info=True)
        finally:
            # Cleanup logic here (e.g., dataset_processor.cleanup_logger())
            pass 
            
    return dataset_reaction_list

def process_single_dataset(dataset_file_path, save_path, onto_file_path, error_log_directory):
    """Processes a single, specified dataset file."""
    dataset_reaction_list = []
    
    try: 
        logger.info(f"Processing single dataset: {dataset_file_path}")

        dataset_1 = rxn_rdf_converter.DatasetProcessor(
            dataset_pb=dataset_pb2,
            dataset_file_path=dataset_file_path,
            owl_onto_file_path=onto_file_path,
            output_directory=save_path,
            error_log_directory=error_log_directory,
            fmt='json-ld'
        )
        _, reaction_error, dataset_reaction_list = dataset_1.extract_reaction(dataset_reaction_list)
        
        logger.info(f"Successfully completed single dataset: {dataset_file_path}")
    except Exception as e: 
        logger.error(f"Failed to process single dataset {dataset_file_path} - Error: {e}", exc_info=True)

    return dataset_reaction_list


# =================================================================
#               MAIN EXECUTION
# =================================================================

def main():
    try: 
        parser = argparse.ArgumentParser(
            description="Process Open Reaction Database (ORD) datasets into RDF."
        )
        
        # Define a Parent Parser for SHARED Arguments
        parent_parser = argparse.ArgumentParser(add_help=False)
        parent_parser.add_argument(
            '--save_path', type=str, required=True, 
            help="Output directory for processed RDF data."
        )
        parent_parser.add_argument(
            '--onto_file_path', type=str, required=True, 
            help="Path to the OWL ontology file (e.g., mdsChemRxn.owl)."
        )
        parent_parser.add_argument(
            '--error_log_directory', type=str, required=True, 
            help="Directory to store logs and CSV."
        )

        # Setup Subparsers - requires one of the following commands
        subparsers = parser.add_subparsers(
            dest='mode', required=True, help='Operational mode: all or single-dataset.'
        )

        # --- Sub-command 1: ALL ---
        parser_all = subparsers.add_parser(
            'all', parents=[parent_parser], 
            help='Process ALL datasets found in the root directory.'
        )
        parser_all.add_argument(
            '--dataset_root', type=str, required=True, 
            help="Root directory containing ORD datasets (to search for all files)."
        )
        
        # --- Sub-command 2: SINGLE DATASET ---
        parser_single_ds = subparsers.add_parser(
            'single-dataset', parents=[parent_parser], 
            help='Process a single, specific dataset file.'
        )
        parser_single_ds.add_argument(
            'dataset_file_path', type=str, 
            help="Full path to the SINGLE dataset file (e.g., /data/ds1.pb)."
        )

        # Parse Arguments
        args = parser.parse_args()

        # --- Execution Dispatch ---
        dataset_reaction_list = []
        
        save_path = args.save_path
        onto_file_path = args.onto_file_path
        error_log_directory = args.error_log_directory

        if args.mode == 'all':
            logger.info("Starting Data Processing in ALL mode....")
            
            # Find all dataset files
            file_list, _, _, _ = setup_file_path(
                args.dataset_root, save_path, onto_file_path, error_log_directory
            )
            # Process them
            dataset_reaction_list = process_all_datasets(
                file_list, save_path, onto_file_path, error_log_directory
            )
            
        elif args.mode == 'single-dataset':
            logger.info("Starting Data Processing in SINGLE-DATASET mode....")
            
            # Process the single specified dataset file
            dataset_reaction_list = process_single_dataset(
                args.dataset_file_path, save_path, onto_file_path, error_log_directory
            )

        # --- Save CSV Results (Applies to both modes) ---
        save_results(dataset_reaction_list, error_log_directory)
        logger.info("Data processing completed successfully.")

    except Exception as e: 
        logger.error(f"Error in main execution: {e}", exc_info=True)
        

if __name__ == '__main__': 
    main()