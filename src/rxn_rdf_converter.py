# Class and method constructions
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
import logging
import csv

# =================================================================
#               NAMESPACE DEFINITIONS
# =================================================================
global AFE, AFR, AFRL, AFQ, OBO, CCO, MDS, NCIT, QUDT, UNIT

AFE = Namespace('http://purl.allotrope.org/ontologies/equipment#')
AFR = Namespace('http://purl.allotrope.org/ontologies/result#')
AFRL = Namespace('http://purl.allotrope.org/ontologies/role#')
AFQ = Namespace('http://purl.allotrope.org/ontologies/quality#')
AFM = Namespace('http://purl.allotrope.org/ontologies/material#')
OBO = Namespace('http://purl.obolibrary.org/obo/')
CCO = Namespace('https://www.commoncoreontologies.org/')
MDS = Namespace('https://cwrusdle.bitbucket.io/mds/')
NCIT = Namespace('http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#')
QUDT = Namespace('http://qudt.org/schema/qudt/')
UNIT = Namespace('http://qudt.org/vocab/unit/')

# Set up logging
logger = logging.getLogger(__name__)

# =================================================================
#               CLASS & METHOD CONSTRUCTION
# =================================================================
        
class DatasetProcessor: 
    """ A class to process an Open Reaction Database (ORD) dataset """
    def __init__ (self, dataset_pb, dataset_file_path, owl_onto_file_path, output_directory, error_log_directory, fmt="turtle"): 
        self.dataset = load_message(dataset_file_path, dataset_pb.Dataset,)
        self.dataset_id = re.split('-', self.dataset.dataset_id)[1]
        dataset_folder = os.path.basename(os.path.dirname(dataset_file_path))
        os.makedirs(os.path.join(output_directory, dataset_folder), exist_ok=True)
        self.output_dir = os.path.join(output_directory, dataset_folder)
        self.owl_onto = owl_onto_file_path
        self.error_log_directory = error_log_directory
        self.fmt = fmt

        # Create dataset-specific logger
        self.dataset_logger, self.log_handler = self._create_dataset_logger()

        # Log initialization
        self.dataset_logger.info(f"DatasetProcessor initiated for dataset: {self.dataset_id}")
        self.dataset_logger.info(f"Dataset file: {dataset_file_path}")
        self.dataset_logger.info(f"Output Directory: {self.output_dir}")
        self.dataset_logger.info(f"Format: {fmt}")
    
    def _create_dataset_logger(self):
        """ Create simple logger for each dataset """

        log_file = f'{self.error_log_directory}/dataset_{self.dataset_id}.log'
        logger_name = f'dataset_{self.dataset_id}'

        # Create logger
        dataset_logger = logging.getLogger(logger_name)
        dataset_logger.setLevel(logging.INFO)
        dataset_logger.handlers.clear() # clear existing handlers

        # Create file handler
        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        dataset_logger.addHandler(handler)

        return dataset_logger, handler
    
    def extract_reaction(self, dataset_reaction_list=None):
        """ A method to extract reaction by instantiate the ReactionKG class """
        if dataset_reaction_list is None: 
            dataset_reaction_list = []

        reaction_error = []
        logger.info(f"Found {len(self.dataset.reactions)} reactions in {self.dataset_id}")

        for ind, reaction in enumerate(self.dataset.reactions): 
            try: 
                self.dataset_logger.info(f"Processing reaction {ind}/{len(self.dataset.reactions)}: {reaction.reaction_id}")

                reaction_data = ReactionKG(reaction, self.fmt).generate_reaction().generate_instances(self.owl_onto).generate_data_graph(self.dataset_id, self.output_dir)
                dataset_reaction_list.append([self.dataset_id, reaction.reaction_id])

                self.dataset_logger.info(f"Successfully processed reaction: {reaction.reaction_id}")
                
            except Exception as e: 
                reaction_error.append([reaction.reaction_id, str(e)])
                self.dataset_logger.error(f"Error processing reaction {reaction.reaction_id}: {e}")
                continue
        
        self.dataset_logger.info(f"Dataset {self.dataset_id} processing completed")
        
        return self, reaction_error, dataset_reaction_list

class ReactionKG: 
    """ A class to process reaction data based on the ord schema (in Google Protocol Buffers format) into Python lists of dictionaries and convert into RDF triples based on MDS-Onto to generate a RDF graph """
    
    def __init__ (self, reaction_pb, fmt="turtle"): 
        self.reaction_pb = reaction_pb
        self.fmt = fmt
        self.reaction_id = re.split('-', reaction_pb.reaction_id)[1] 
        
        # Initialize core lists to handle parts of a reaction data
        self.reaction_identifiers = []
        self.reaction_inputs = []
        self.reaction_setup = []
        self.reaction_conditions = []
        self.reaction_notes = []
        self.reaction_workups = []
        self.reaction_outcomes = []
        self.reaction_provenance = []

        # Initialize core instance variables for handling 
        self.chemical_reaction = None
        self.reaction_mixture = None
        self.reaction_environment = None
        self.crude_product = None

    
    def generate_reaction (self): 
        """ Main method that takes a reaction as an attribute of an instance object and generate lists of: 
            reaction identifiers dict
            reaction inputs dict
            reaction setup dict
            reaction conditions dict
            reaction outcomes dict

        """
        # reaction identifiers 
        for identifier in self.reaction_pb.identifiers:
            identifier_dict = {'reactionID':self.reaction_id}
            identifier_dict.update({item: None for item in list(identifier.DESCRIPTOR.fields_by_name)})
            identifier_dict.update(message_to_row(identifier))
            self.reaction_identifiers.append(identifier_dict)

        # reaction inputs, components, and compound identifiers
        for input in self.reaction_pb.inputs: 
            has_whitespace = any(char.isspace() for char in str(input.strip()))
            if has_whitespace == True: 
                input_string = re.sub(r'\s+', '', str(input.strip()))
            else: 
                input_string = input.strip()
            reaction_input_dict = {'reactionID':self.reaction_id, 'InputKey':input_string}
            reaction_input_dict.update({item: None for item in list(self.reaction_pb.inputs[input].DESCRIPTOR.fields_by_name)})
            reaction_input_dict.update(message_to_row(self.reaction_pb.inputs[input]))
            
            reaction_input_dict['components'] = True if self.reaction_pb.inputs[input].components else None
            reaction_input_dict['crude_components'] = True if self.reaction_pb.inputs[input].crude_components else None
            reaction_input_dict['addition_time'] = True if self.reaction_pb.inputs[input].addition_time else None
            reaction_input_dict['addition_speed'] = True if self.reaction_pb.inputs[input].addition_speed else None
            reaction_input_dict['addition_duration'] = True if self.reaction_pb.inputs[input].addition_duration else None
            reaction_input_dict['addition_device'] = True if self.reaction_pb.inputs[input].addition_device else None
            reaction_input_dict['addition_temperature'] = True if self.reaction_pb.inputs[input].addition_temperature else None
            reaction_input_dict['flow_rate'] = True if self.reaction_pb.inputs[input].flow_rate else None
            reaction_input_dict['texture'] = True if getattr(self.reaction_pb.inputs[input], 'texture') else None
            
            for ind, component in enumerate(self.reaction_pb.inputs[input].components): 
                identifier_list, inchi_key = self.generate_compound_identifiers(component.identifiers)
                if inchi_key:
                    reaction_input_dict[f"components[{ind}].INCHI_KEY"] = inchi_key
                if component.amount:
                    reaction_input_dict[f"components[{ind}].amount.type"] = component.amount.WhichOneof('kind')

            self.reaction_inputs.append(reaction_input_dict)

        # reaction setup
        setup_dict = {'reactionID':self.reaction_id}
        setup_dict.update({item:None for item in list(self.reaction_pb.setup.DESCRIPTOR.fields_by_name)})
        setup_dict.update(message_to_row(self.reaction_pb.setup))
        if self.reaction_pb.setup.vessel: 
            setup_dict['vessel'] = True
        if self.reaction_pb.setup.environment:
            setup_dict['environment'] = True

        self.reaction_setup.append(setup_dict)

        # reaction conditions
        conditions_dict = {'reactionID':self.reaction_id}
        conditions_dict.update({item:None for item in list(self.reaction_pb.conditions.DESCRIPTOR.fields_by_name)})
        conditions_dict.update(message_to_row(self.reaction_pb.conditions))
        if self.reaction_pb.conditions.temperature: 
            conditions_dict['temperature'] = True
        if self.reaction_pb.conditions.pressure: 
            conditions_dict['pressure'] = True
        if self.reaction_pb.conditions.stirring: 
            conditions_dict['stirring'] = True
        if self.reaction_pb.conditions.electrochemistry: 
            conditions_dict['electrochemistry'] = True
        if self.reaction_pb.conditions.flow: 
            conditions_dict['flow'] = True
        if self.reaction_pb.conditions.illumination: 
            conditions_dict['illumination'] = True
        self.reaction_conditions.append(conditions_dict)

        # reaction notes
        notes_dict = {'reactionID': self.reaction_id}
        notes_dict.update({item:None for item in list(self.reaction_pb.notes.DESCRIPTOR.fields_by_name)})
        notes_dict.update(message_to_row(self.reaction_pb.notes))
        self.reaction_notes.append(notes_dict)

        # reaction workups
        for ind, workup in enumerate(self.reaction_pb.workups):
            workups_dict = {'reactionID':self.reaction_id, 'Index':ind}
            workups_dict.update({item:None for item in list(workup.DESCRIPTOR.fields_by_name)})
            workups_dict.update(message_to_row(workup))
            if workup.input: 
                workups_dict['input'] = True
                for ind, component in enumerate(workup.input.components):
                    workup_input_list, inchi_key = self.generate_compound_identifiers(component.identifiers)
                    workups_dict[f"input.components[{ind}].INCHI_KEY"] = inchi_key
                    workups_dict['InputKey'] = f"workup_{0}"
            if workup.amount: 
                workups_dict['amount'] = True
                workups_dict[f"amount.type"] = workup.amount.WhichOneof('kind')

            if workup.duration: 
                workups_dict['duration'] = True
            if workup.temperature:
                workups_dict['temperature'] = True
            if workup.stirring:
                workups_dict['stirring'] = True

            self.reaction_workups.append(workups_dict)
        
        # reaction outcomes: 
        for ind, outcome in enumerate(self.reaction_pb.outcomes):
            outcome_dict = {'reactionID':self.reaction_id, 'Index': ind}
            outcome_dict.update({item:None for item in list(outcome.DESCRIPTOR.fields_by_name)})
            outcome_dict.update(message_to_row(outcome))
            
            if outcome.reaction_time: 
                outcome_dict['reaction_time'] = True
            if outcome.analyses: 
                outcome_dict['analyses'] = True
            if outcome.conversion:
                outcome_dict['conversion'] = True
            if outcome.products:
                outcome_dict['products'] = True
                for index, product in enumerate(outcome.products):
                    product_identifier, inchi_key = self.generate_compound_identifiers(product.identifiers)
                    outcome_dict[f"products[{index}].INCHI_KEY"] = inchi_key

            self.reaction_outcomes.append(outcome_dict)

        provenance_dict = {'reactionID': self.reaction_id}
        provenance_dict.update({item:None for item in list(self.reaction_pb.provenance.DESCRIPTOR.fields_by_name)})
        provenance_dict.update(message_to_row(self.reaction_pb.provenance))
        
        return self
        
    def generate_compound_identifiers (self, identifiers): 
        identifier_list = []
        desired_type = set()
        
        for identifier in identifiers:
            identifier_dict = {item:None for item in list(identifier.DESCRIPTOR.fields_by_name)}
            identifier_dict.update(message_to_row(identifier))
            desired_type.add(message_to_row(identifier)['type'])
            identifier_list.append(identifier_dict)
        
        inchi_key = None
        
        if 'INCHI_KEY' in desired_type: 
            inchi_key = next((item['value'] for item in identifier_list if item['type'] == 'INCHI_KEY'), None)
        elif 'INCHI' in desired_type: 
            inchi = next((item['value'] for item in identifier_list if item['type'] == 'INCHI'), None)
            rdkit_mol = Chem.MolFromInchi(inchi)
            if rdkit_mol:
                identifier_list.append({'type':'INCHI_KEY', 'details':None, 'value': Chem.MolToInchiKey(rdkit_mol)})
                inchi_key = Chem.MolToInchiKey(rdkit_mol)
            if 'SMILES' not in desired_type:
                identifier_list.append({'type': 'SMILES', 'details':None, 'value': Chem.MolToSmiles(rdkit_mol)})
        elif 'SMILES' in desired_type:
            smiles = next((item['value'] for item in identifier_list if item['type'] == 'SMILES'), None)
            rdkit_mol = Chem.MolFromSmiles(smiles)
            if rdkit_mol:
                identifier_list.append({'type':'INCHI', 'details':None, 'value':Chem.MolToInchi(rdkit_mol)})
                identifier_list.append({'type':'INCHI_KEY', 'details':None, 'value':Chem.MolToInchiKey(rdkit_mol)})
                inchi_key = Chem.MolToInchiKey(rdkit_mol)
            if 'CXSMILES' not in desired_type:
                identifier_list.append({'type': 'CXSMILES', 'details':None, 'value': Chem.MolToCXSmiles(rdkit_mol)})
        
        for ind, identifier in enumerate(identifier_list): 
            identifier_list[ind]['INCHI_KEY'] = inchi_key

        return identifier_list, inchi_key
    
    def generate_instances (self, mds_file_path): 
        """Main method orchestrates all other methods - coordinates all instance generation"""
        
        try: 
            self._initialize_instance_dict(mds_file_path)
        except Exception as e: 
            logger.error(f"Failed to initialize dictionaries of instances for reaction {self.reaction_id}: {e}")
            pass
        
        try:  
            self._process_reaction_identifiers() if self.reaction_identifiers else None
        except Exception as e: 
            logger.error(f"Failed to process reaction identifiers into RDF triples for reaction {self.reaction_id}: {e}")
            pass
        
        try: 
            if self.reaction_inputs: 
                self._process_reaction_inputs()
        except Exception as e: 
            logger.error(f"Failed to process reaction inputs into RDF triples for reaction {self.reaction_id}: {e}")
            pass   
        
        try:
            if self.reaction_conditions:
                self._process_reaction_conditions()
        except Exception as e: 
            logger.error(f"Failed to process reaction conditions into RDF triples for reaction {self.reaction_id}: {e}")
            pass
        
        try: 
            if self.reaction_notes:
                self._process_reaction_notes()
        except Exception as e: 
            logger.error(f"Failed to process reaction notes into RDF triples for reaction {self.reaction_id}: {e}")
            pass
        
        try: 
            if self.reaction_workups:
                self._process_reaction_workups()
        except Exception as e: 
            logger.error(f"Failed to process reaction workups into RDF triples for reaction {self.reaction_id}: {e}")
            pass
        
        try: 
            if self.reaction_outcomes: 
                self._process_reaction_outcomes()
        except Exception as e: 
            logger.error(f"Failed to process reaction outcomes into RDF triples for reaction {self.reaction_id}: {e}")
            pass
        
        try: 
            if self.reaction_setup:
                self._process_reaction_setup()
        except Exception as e: 
            logger.error(f"Failed to process reaction identifiers into RDF triples for reaction {self.reaction_id}: {e}")
            pass

        return self

    def _initialize_instance_dict(self, mds_file_path): 
        """ Initialization of: 
            1. the variables for OWL 2 ontology and the get the namespaces
            2. the instance dictionary with all property labels and 
            3. create core reaction instances
        
        Args: 
            
        Returns: 
            Object attributes: 
                1. mds, cco, and obo namespaces
                2. instance_dict
                3. chemical_reaction, reaction_mixture, reaction_environment, crude_product
        
        Generates a dictionary where the keys are human-readable labels of object/datatype properties, and the values are
        2-tuples that contain the URI of that property in the first entry and the type (object/datatype) in second entry.

        Parameters
        ----------
        mds_ontology_file_path : str
            Path to the RDF/OWL ontology file.

        Returns
        -------
        dict
            Dictionary of the form:
            {
                "has material": ("http://example.org/ontology#hasMaterial", "Object Property"),
                "has value": ("http://example.org/ontology#hasValue", "Datatype Property"),
                ...
            }

        Generate OWL 2 Onto attribute to store classes, properties, and instances
        Args: 
            file path to the ontology file in .owl (rdf/xml) format

        Returns: 
            onto attribute with classes, properties, and instances that can be called by other methods
        """

        self.onto = get_ontology(mds_file_path).load()

        self.onto.get_namespace('http://purl.allotrope.org/ontologies/equipment#')
        self.onto.get_namespace('http://purl.allotrope.org/ontologies/result#')
        self.onto.get_namespace('http://purl.allotrope.org/ontologies/role#')
        self.onto.get_namespace('http://purl.allotrope.org/ontologies/quality#')
        self.onto.get_namespace('http://purl.allotrope.org/ontologies/material#')
        self.onto.get_namespace('http://ncicb.nci.nih.gov/xml/owl/EVS/Thesaurus.owl#')
        cco = self.onto.get_namespace('https://www.commoncoreontologies.org/')
        mds = self.onto.get_namespace('https://cwrusdle.bitbucket.io/mds/')
        qudt = self.onto.get_namespace('http://qudt.org/schema/qudt/')
        rdf = self.onto.get_namespace('http://www.w3.org/1999/02/22-rdf-syntax-ns#type')
        obo = self.onto.get_namespace('http://purl.obolibrary.org/obo/')
        unit = self.onto.get_namespace('http://qudt.org/vocab/unit/')
        
        # Store namespaces as instance variables for use in other methods
        self.mds = mds
        self.cco = cco
        self.obo = obo
        self.qudt = qudt
        self.unit = unit
        
        # Setup and store units in unit_mapping dictionary to be used in other methods
        self.unit_mapping = {        
        'CELSIUS': self.unit.DEG_C,
        'FAHRENHEIT': self.unit.DEG_F,
        'KELVIN': self.unit.K,

        'KILOGRAM': self.unit.KiloGM,
        'GRAM': self.unit.GM,
        'MILLIGRAM': self.unit.MilliGM,
        'MICROGRAM': self.unit.MicroGM,

        'MOLE': self.unit.MOL,
        'MILLIMOLE' : self.unit.MilliMOL,
        'MICROMOLE': self.unit.MicroMOL,
        'NANOMOLE': self.unit.NanoMOL,

        'LITER': self.unit.L,
        'MILLILITER': self.unit.MilliL,
        'MICROLITER': self.unit.MicroL,
        'NANOLITER': self.unit.NanoL,

        'BAR': self.unit.BAR,
        'ATMOSPHERE': self.unit.ATM,
        'PSI': self.unit.PSI,
        'KPSI': f'{self.unit}KiloLB_F-PER-IN2',
        'PASCAL': self.unit.PA,
        'KILOPASCAL': self.unit.KiloPA,
        'TORR': self.unit.TORR,
        'MM_HG': self.unit.MilliM_HG,

        'DAY': self.unit.DAY,
        'HOUR': self.unit.HR,
        'MINUTE': self.unit.MIN,
        'SECOND': self.unit.SEC,

        'PERCENTAGE': self.unit.PERCENT,

        }

        self.prop_metadata_dict = {}

        for prop in self.onto.object_properties():
            if prop.label: 
                self.prop_metadata_dict[str(prop.label[0])] = (prop.iri, 'Object Property')
        for prop in self.onto.annotation_properties(): 
            if prop.label: 
                self.prop_metadata_dict[str(prop.label[0])] = (prop.iri, 'Annotation Property')
        for prop in self.onto.data_properties():
            if prop.label: 
                self.prop_metadata_dict[str(prop.label[0])] = (prop.iri, 'Datatype Property')
        self.prop_metadata_dict['type'] = ('http://www.w3.org/1999/02/22-rdf-syntax-ns#type', 'Object Property')
        
        self.instance_dict = {}
        
        for prop in self.onto.properties():
            self.instance_dict[str(prop.label[0])] = []
        self.instance_dict['type'] = []

        self.chemical_reaction = self.mds.ChemicalReaction('ChemicalReaction#' + self.reaction_id)
        self.reaction_mixture = self.mds.ReactionMixture('ReactionMixture#' + self.reaction_id)
        self.reaction_environment = self.mds.ReactionEnvironment('ReactionEnvironment#' + self.reaction_id)
        self.crude_product = self.mds.CrudeProduct('CrudeProduct#' + self.reaction_id)
        
        self.instance_dict['type'].append([self.chemical_reaction.iri, self.chemical_reaction.is_instance_of[0].iri])
        self.instance_dict['type'].append([self.reaction_mixture.iri, self.reaction_mixture.is_instance_of[0].iri])
        self.instance_dict['type'].append([self.reaction_environment.iri, self.reaction_environment.is_instance_of[0].iri])
        self.instance_dict['type'].append([self.crude_product.iri, self.crude_product.is_instance_of[0].iri])
        
        self.instance_dict['is input of'].append([self.reaction_mixture.iri, self.chemical_reaction.iri])
        self.instance_dict['environs'].append([self.reaction_environment.iri, self.chemical_reaction.iri])
        self.instance_dict['has output'].append([self.chemical_reaction.iri, self.crude_product.iri]) 
    
    def _extract_index_set(self, item_dict, pattern): 
        """ Extract indices from dictionary keys based on regex pattern
        
        Args: 
            item_dict: Dictionary containing keys to search
            pattern: Regex pattern to match keys and extract indices 

        Returns: 
            list: sorted list of integers for numeric indices or list of unique strings for text indices
        """
        
        numeric_indices = []
        string_indices = []
        compiled_pattern = re.compile(pattern)

        for key in item_dict.keys():
            match = compiled_pattern.match(key)
            if match:
                captured_value = match.group(1)
                if captured_value.isdigit():
                    numeric_indices.append(int(captured_value))
                else: 
                    string_indices.append(captured_value)
        
        if numeric_indices:
            result = sorted(set(numeric_indices))
        elif string_indices:
            result = list(set(string_indices)) # remove duplicates but preserve as strings 
        else: 
            result = []
        return result

    def _process_reaction_identifiers(self):
        """Process reaction identifiers"""
        identifier_mapping = {
            'UNSPECIFIED': None,
            'CUSTOM': None,
            'REACTION_TYPE': self.mds.ReactionType,
            'REACTION_SMILES': self.mds.ReactionSMILES,
            'REACTION_CXSMILES': self.mds.ReactionCXSMILES,
            'RDFILE': self.mds.ReactionDataFile,
            'RINCHI': self.mds.RInChI,
        }
            
        for item in self.reaction_identifiers: 
            identifier_class = identifier_mapping.get(item['type'])
            if identifier_class == None:
                identifier_type = 'Custom'
            else: 
                identifier_type = re.sub(r'\s+', '', str(identifier_class.label[0]))
            if identifier_class: 
                identifier = identifier_class(f"{identifier_type}{self.reaction_id}")
                self.instance_dict['type'].append([identifier.iri, identifier.is_instance_of[0].iri])
                self.instance_dict['designates'].append([identifier.iri, self.chemical_reaction.iri])
                self.instance_dict['has text value'].append([identifier.iri, item['value']])
                if item['details']: 
                    self.instance_dict['details'].append([identifier.iri, item['details']])
                if item['is_mapped']:
                    self.instance_dict['is mapped'].append([identifier.iri, item['is_mapped']])
            
    def _extract_components(self, component_list, process_node=None, context='input'):
        """ Private method to process the components in reaction inputs and reaction workups as well as products 
        
        Args: 
            Takes in the list of dictionaries for a component, the process node (input_addition or reaction_workup)

        Returns: 
            lists in instance_dict with component properties 
        """
        product_dict = {}

        if context == 'input': 
            pattern = rf'^components\[(\d+)\]'
            prefix = 'components'
        elif context == 'workup':
            pattern = rf'^input\.components\[(\d+)\]'
            prefix = 'input.components'
        elif context == 'product':
            pattern = rf'^products\[(\d+)\]'
            prefix = 'products'
        else: 
            raise ValueError(f"Invalid context: {context}. Must be 'input', 'workup', or 'product'")

        if (context == 'input' or context == 'workup') and process_node is None: 
            raise ValueError(f"A process node (input addition process or workup process) is required for components in 'input' and 'workup'")
        
        index_list = self._extract_index_set(component_list, pattern)
        
        for i in index_list:
            if (context == 'input' or context == 'workup'): 
                if f"{prefix}[{i}].INCHI_KEY" in component_list: 
                    reaction_component = self.mds.Component(f'Component#{self.reaction_id}_{component_list["InputKey"]}_{component_list[f"{prefix}[{i}].INCHI_KEY"]}')
                else:
                    reaction_component = self.mds.Component(f'Component#{self.reaction_id}_{component_list["InputKey"]}_component_{i}')
            elif context == 'product':
                if f"{prefix}[{i}].INCHI_KEY" in component_list: 
                    reaction_component = self.mds.Product(f'Product#{self.reaction_id}_{component_list[f"{prefix}[{i}].INCHI_KEY"]}')
                else: 
                    reaction_component = self.mds.Product(f'Product#{self.reaction_id}_Outcome_{component_list["Index"]}_Product_{i}')
                product_dict[i] = reaction_component
            
            self.instance_dict['type'].append([reaction_component.iri, reaction_component.is_instance_of[0].iri])

            if context == 'input' and self.reaction_mixture: 
                self.instance_dict['member part of'].append([reaction_component.iri, self.reaction_mixture.iri])
            if (context == 'input' or context == 'workup') and process_node:
                self.instance_dict['is input of'].append([reaction_component.iri, process_node.iri])
            elif context == 'product': 
                self.instance_dict['has output'].append([self.chemical_reaction.iri, reaction_component.iri])
            
            try:
                self._extract_compound_identifiers(component_list, i, context, prefix, reaction_component)
            except Exception as e:
                logger.error(f"Failed to process {context} {prefix} identifier {i} into RDF triples for reaction {self.reaction_id}: {e}")
                continue

            role_key = f"{prefix}[{i}].reaction_role"
            if role_key in component_list: 
                role_mapping = {
                    'UNSPECIFIED': None,
                    'REACTANT': self.cco.ont00000808,
                    'REAGENT': self.mds.ReagentArtifactFunction,
                    'SOLVENT': self.cco.ont00001160,
                    'CATALYST': self.cco.ont00000267,
                    'WORKUP': self.mds.WorkupArtifactFunction,
                    'INTERNAL_STANDARD': self.mds.InternalStandardArtifactFunction,
                    'AUTHENTIC_STANDARD': self.mds.AuthenticStandardArtifactFunction,
                    'PRODUCT': self.mds.Product,
                    'BYPRODUCT': self.mds.Byproduct,
                    'SIDE_PRODUCT': self.mds.SideProduct
                }
                role_class = role_mapping.get(component_list[role_key])
                if role_class == None:
                    role_type = 'Custom'
                else: 
                    role_type = re.sub(r'\s+', '', str(role_class.label[0]))
                if context == 'input' or context == 'workup':
                    role_id_base = f'{self.reaction_id}_{component_list["InputKey"]}_component_{i}'
                elif context == 'product':
                    role_id_base = f'{self.reaction_id}_Outcome_{component_list["Index"]}_Product_{i}'
                
                if role_class: 
                    component_role = role_class(f"{role_type}ArtifactFunction#{role_id_base}")
                    self.instance_dict['type'].append([component_role.iri, component_role.is_instance_of[0].iri])
                    self.instance_dict['inheres in'].append([component_role.iri, reaction_component.iri])

            if f"{prefix}[{i}].is_limiting" in component_list and component_list[f"{prefix}[{i}].is_limiting"]: 
                self.instance_dict['is limiting'].append([reaction_component.iri, component_list[f"{prefix}[{i}].is_limiting"]])
            
            if f"{prefix}[{i}].is_desired_product" in component_list: 
                self.instance_dict['is desired product'].append([reaction_component.iri, component_list[f"{prefix}[{i}].is_desired_product"]])
            
            if f"{prefix}[{i}].isolated_color" in component_list:
                self.instance_dict['isolated color'].append([reaction_component.iri, component_list[f"{prefix}[{i}].isolated_color"]])
            
            if f"{prefix}[{i}].amount.type" in component_list: 
                self._extract_component_amount(component_list, reaction_component, i, context)
            
            if 'addition_temperature.value' in component_list and 'addition_temperature.units' in component_list:
                self._process_temperature(component_list, reaction_component, context)
        
        return self, product_dict

    def _extract_component_amount (self, component_list, reaction_component, i=None, context='input'): 
        """ Process compound amounts in reaction inputs and reaction workups """

        if context == 'input': 
            prefix = f"components[{i}]."
        elif context == 'workup':
            prefix = f"input.components[{i}]."
        elif context == 'aliquot': 
            prefix = ''
        else: 
            raise ValueError(f"Invalid context: {context}. Must be 'input', 'workup', 'aliquot'")
        
        if f"{prefix}amount.type" in component_list:
            amount_type = component_list[f"{prefix}amount.type"] 
            if amount_type == 'moles':
                amount_added = self.cco.ont00000768(f'Moles#{self.reaction_id}_{component_list["InputKey"]}_component_{i}')
            elif amount_type == 'mass':
                amount_added = self.cco.ont00000768(f'Mass#{self.reaction_id}_{component_list["InputKey"]}_component_{i}')
            elif amount_type == 'volume':
                amount_added = self.cco.ont00000768(f'Volume#{self.reaction_id}_{component_list["InputKey"]}_component_{i}')
            else: amount_added = None
            if amount_added and f"{prefix}amount.{amount_type}.value" in component_list and f"{prefix}amount.{amount_type}.units" in component_list: 
                self.instance_dict['inheres in'].append([amount_added.iri, reaction_component.iri])
                self.instance_dict['has decimal value'].append([amount_added.iri, component_list[f"{prefix}amount.{amount_type}.value"]])
                self.instance_dict['uses measurement unit'].append([amount_added.iri, self.unit_mapping[component_list[f"{prefix}amount.{amount_type}.units"]].iri])

    def _extract_compound_identifiers(self, component_list, i, context, prefix, reaction_component):
        """ Extract compound identifiers for components in a reaction either as reaction inputs or reaction workups """
        
        pattern = rf'^{re.escape(prefix)}\[{i}\]\.identifiers\[(\d+)\]\.(.+)$'
        ind_list = self._extract_index_set(component_list, pattern)
        if not ind_list: 
            return

        if context == 'input' or context == 'workup': 
            identifier_id_base = f'{self.reaction_id}_{component_list["InputKey"]}_component_{i}'
        elif context == 'product': 
            identifier_id_base = f'{self.reaction_id}_Outcome_{component_list["Index"]}_Product_{i}'

        identifier_mapping = {
            'UNSPECIFIED': self.cco.ont00000649,
            'CUSTOM': self.cco.ont00000649,
            'SMILES': self.mds.SMILES,
            'INCHI': self.mds.InChI,
            'MOLBLOCK': self.mds.MolBlock,
            'IUPAC_NAME': self.mds.IUPACName,
            'NAME': self.mds.CompoundName, 
            'CAS_NUMBER': self.mds.CASNumber,
            'PUBCHEM_CID': self.mds.PubChemCompoundIdentifier,
            'CHEMSPIDER_ID': self.mds.ChemSpiderIdentifier,
            'CXSMILES': self.mds.CXSMILES,
            'INCHI_KEY': self.mds.InChIKey,
            'XYZ': self.mds.XYZ, 
            'UNIPROT_ID': self.mds.UniProtIdentifier,
            'PDB_ID': self.mds.ProteinDataBankIdentifier,
            'AMINO_ACID_SEQUENCE': self.mds.AminoAcidSequence, 
            'HELM': self.mds.HierarchicalEditingLanuageForMacromolecules, 
            'MDL': self.mds.MolecularDeisgnLimitedNumber,
        }

        for j in ind_list: 
            try:
                if f"{prefix}[{i}].identifiers[{j}].type" in component_list and f"{prefix}[{i}].identifiers[{j}].value" in component_list:
                    identifier_type_string = component_list[f"{prefix}[{i}].identifiers[{j}].type"]
                    identifier_class = identifier_mapping.get(component_list[f"{prefix}[{i}].identifiers[{j}].type"])
                    if identifier_type_string in['UNSPECIFIED', 'CUSTOM']:
                        identifier_type = 'Custom'
                    else:
                        identifier_type = re.sub(r'\s+', '', str(identifier_class.label[0]))
                    if identifier_class:
                        identifier = identifier_class(f"{identifier_type}#{identifier_id_base}_identifier_{j}")
                        self.instance_dict['designates'].append([identifier.iri, reaction_component.iri])
                        self.instance_dict['type'].append([identifier.iri, identifier.is_instance_of[0].iri])
                        self.instance_dict['has text value'].append([identifier.iri, component_list[f"{prefix}[{i}].identifiers[{j}].value"]])
                        if f"{prefix}[{i}].identifiers[{j}].details" in component_list: 
                            self.instance_dict['details'].append([identifier.iri, component_list[f"{prefix}[{i}].identifiers[{j}].details"]])
                        if f"{prefix}[{i}].identifiers[{j}].is_mapped" in component_list: 
                            self.instance_dict['is mapped'].append([identifier.iri, component_list[f"{prefix}[{i}].identifiers[{j}].is_mapped"]])
            except Exception as e: 
                logger.error(f"Failed to process identifier {j} for {prefix} {i} in reaction {self.reaction_id}: {e}")
                continue

    def _process_reaction_inputs(self):
        """ Process reaction input details such as addition speed, addition duration, addition time, addition speed, addition temperature """
        
        for item in self.reaction_inputs:
            input_addition = self.mds.InputAddition(f'InputAddition#{self.reaction_id}_{item["InputKey"]}')
            self.instance_dict['type'].append([input_addition.iri, input_addition.is_instance_of[0].iri])
            self.instance_dict['has process part'].append([self.chemical_reaction.iri, input_addition.iri])
            self.instance_dict['has output'].append([input_addition.iri, self.reaction_mixture.iri])
            
            self._extract_components(item, input_addition)
            
            if item['addition_order']: 
                self.instance_dict['addition order'].append([input_addition.iri, item['addition_order']])
            
            if item['addition_speed'] == True and 'addition_speed.type' in item: 
                addition_speed = self.mds.AdditionSpeed(f'AdditionSpeed#{self.reaction_id}')
                self.instance_dict['type'].append([addition_speed.iri, addition_speed.is_instance_of[0].iri])
                self.instance_dict['has occurent part'].append([input_addition.iri, addition_speed.iri])
                self.instance_dict['has text value'].append([addition_speed.iri, item['addition_speed.type']])
                if 'addition_speed.details' in item: 
                    self.instance_dict['details'].append([addition_speed.iri, item['addition_speed.details']])
            
            if item['addition_duration'] == True and 'addition_duration.value' in item and 'addition_duration.units' in item:
                addition_duration = self.obo.BFO_0000202(f'AdditionDuration#{self.reaction_id}')
                self.instance_dict['type'].append([addition_duration.iri, addition_duration.is_instance_of[0].iri])
                self.instance_dict['occupies temporal region'].append([input_addition.iri, addition_duration.iri])
                self.instance_dict['has datetime value'].append([addition_duration.iri, item['addition_duration.value']])
                self.instance_dict['uses measurement unit'].append([addition_duration.iri, self.unit_mapping[item['addition_duration.units']].iri])
            
            if item['addition_time'] == True and 'addition_time.value' in item and 'addition_time.units' in item:
                addition_time = self.obo.BFO_0000202(f'AdditionTime#{self.reaction_id}')
                self.instance_dict['type'].append([addition_time.iri, addition_time.is_instance_of[0].iri])
                self.instance_dict['occupies temporal region'].append([input_addition.iri, addition_time.iri])
                self.instance_dict['has datetime value'].append([addition_time.iri, item['addition_time.value']])
                self.instance_dict['uses measurement unit'].append([addition_time.iri, self.unit_mapping[item['addition_time.units']].iri])

            if item['flow_rate'] == True and 'flow_rate.value' in item and 'flow_rate.units' in item:
                flow_rate = self.mds.FlowRate(f'FlowRate#{self.reaction_id}')
                self.instance_dict['type'].append([flow_rate.iri, flow_rate.is_instance_of[0].iri])
                self.instance_dict['has occurent part'].append([input_addition.iri, flow_rate.iri])
                self.instance_dict['has decimal value'].append([flow_rate.iri, item['flow_rate.value']])
                self.instance_dict['uses measurement unit'].append([flow_rate.iri, item['flow_rate.units']])

            if item['addition_device'] == True and 'addition_device.type' in item:
                addition_device = self.cco.ont00000581(f"AdditionDevice#{self.reaction_id}_{item['addition_device.type']}")
                self.instance_dict['type'].append([addition_device.iri, addition_device.is_instance_of[0].iri])
                self.instance_dict['participates in'].append([addition_device.iri, input_addition.iri])
                if 'addition_device.details' in item: 
                    self.instance_dict['details'].append([addition_device.iri, item['addition_device.details']])
    
    def _process_temperature(self, component_list, reaction_component=None, context='condition'):
        """ Process and extract temperature information about a reaction input, the reaction environment, or a reaction workup"""

        if context == 'condition' or context == 'workup': 
            prefix = 'temperature.setpoint'
            pattern = rf'^temperature\.measurements\[(\d+)\]'
            
            if context == 'workup': 
                reaction_temperature = self.cco.ont00000441(f'WorkupTemperature#{self.reaction_id}')
            elif context == 'condition':
                reaction_temperature = self.mds.ReactionTemperature(f'ReactionTemperature#{self.reaction_id}')
        
        elif context == 'input':
            reaction_temperature = self.cco.ont00000441(f'AdditionTemperature#{self.reaction_id}')
            prefix = 'addition_temperature'
        else: 
            raise ValueError(f"Invalid context: {context}. Must be 'condition', 'input', or 'workup'")

        if context == 'input' and reaction_component is None: 
            raise ValueError(f"If _process_temperature method is used with {context}, reaction_component must be provided.")

        self.instance_dict['type'].append([reaction_temperature.iri, reaction_temperature.is_instance_of[0].iri])
        temperature_measurement = self.mds.AnalyticalResult(f'TemperatureMeasurement#{self.reaction_id}_{context}')
        self.instance_dict['type'].append([temperature_measurement.iri, temperature_measurement.is_instance_of[0].iri])
        self.instance_dict['is about'].append([temperature_measurement.iri, reaction_temperature.iri])

        self.instance_dict['has decimal value'].append([temperature_measurement.iri, component_list[f'{prefix}.value']])
        self.instance_dict['uses measurement unit'].append([temperature_measurement.iri, self.unit_mapping[component_list[f'{prefix}.units']].iri])
        
        if f'{prefix}.control.type' in component_list and (context == 'input' or context == 'workup'): 
            control_tool = self.cco.ont00000581(f'TemperatureControl#{self.reaction_id}_{component_list[f"{prefix}.control.type"]}')
            self.instance_dict['type'].append([control_tool.iri, control_tool.is_instance_of[0].iri])
            self.instance_dict['affects'].append([control_tool, reaction_temperature.iri])
        
        if context == 'condition':
            self.instance_dict['bearer of'].append([self.reaction_environment.iri, reaction_temperature.iri])
        elif context == 'workup': 
            self.instance_dict['bearer of'].append([self.reaction_mixture.iri, reaction_temperature.iri])
        elif context == 'input':
            self.instance_dict['bearer of'].append([reaction_component.iri, reaction_temperature.iri])

        if context == 'condition' or context == 'workup': 
            index_list = self._extract_index_set(component_list, pattern)
            
            if context == 'condition' or context == 'workup': 
                prefix = 'temperature.measurements'

            for i in index_list: 
                type_key = f"{prefix}[{i}].type"
                details_key = f"{prefix}[{i}].details"
                time_key = f"{prefix}[{i}].time.value"
                temp_key = f"{prefix}[{i}].temperature.value"
                if type_key in component_list: 
                    temperature_measurement = self.mds.AnalyticalResult(f'TemperatureMeasurement#{self.reaction_id}_{context}_{i}')
                if type_key in component_list and temperature_measurement: 
                    self.instance_dict['details'].append([temperature_measurement.iri, component_list[type_key]])
                if details_key in component_list and temperature_measurement:
                    self.instance_dict['details'].append([temperature_measurement.iri, component_list[details_key]])
                if temp_key in component_list and f"{prefix}[{i}].temperature.units" in component_list: 
                    self.instance_dict['has decimal value'].append([temperature_measurement.iri, component_list[temp_key]])
                    self.instance_dict['uses measurement unit'].append([temperature_measurement.iri, self.unit_mapping[component_list[f"{prefix}[{i}].temperature.units"]].iri])
                if time_key in component_list and f"{prefix}[{i}].time.units" in component_list: 
                    time_measurement = self.obo.BFO_0000148(f"TemperatureMeasurementTime#{self.reaction_id}_{context}_{i}")
                    self.instance_dict['has datetime value'].append([time_measurement.iri, component_list[time_key]])
                    self.instance_dict['uses measurement unit'].append([temperature_measurement.iri, self.unit_mapping[component_list[f"{prefix}[{i}].time.units"]].iri])
                
    def _process_reaction_conditions(self):
        """ Process reaction conditions such as reaction temperature, reaction pressure, and other processes such as stirring, illumination, electrochemistry, and flow """
        atmosphere_mapping = {      
            'UNSPECIFIED': None,
            'CUSTOM': None,
            'AIR': 'AmbientAir',
            'NITROGEN': 'Nitrogen',
            'ARGON': 'Argon',
            'OXYGEN': 'Oxygen',
            'HYDROGEN': 'Hydrogen',
            'CARBON_MONOXIDE': 'CarbonMonoxide',
            'CARBON_DIOXIDE': 'CarbonDioxide',
            'METHANE': 'Methane',
            'AMMONIA': 'Ammonia', 
            'OZONE': 'Ozone',
            'ETHYLENE': 'Ethylene',
            'ACETYLENE': 'Acetyle'
        }
        for item in self.reaction_conditions: 
            if item['temperature'] == True and 'temperature.setpoint.value' in item and 'temperature.setpoint.units' in item: 
                self._process_temperature(item, context='condition')
                
            if item['pressure'] == True and 'pressure.setpoint.value' in item and 'pressure.setpoint.units' in item: 
                reaction_pressure = self.mds.ReactionPressure(f'ReactionPressure#{self.reaction_id}')
                self.instance_dict['type'].append([reaction_pressure.iri, reaction_pressure.is_instance_of[0].iri])
                self.instance_dict['has decimal value'].append([reaction_pressure.iri, item['pressure.setpoint.value']])
                self.instance_dict['bearer of'].append([self.reaction_environment.iri, reaction_pressure.iri])
                self.instance_dict['uses measurement unit'].append([reaction_pressure.iri, self.unit_mapping[item['pressure.setpoint.units']].iri])
                self.instance_dict['bearer of'].append([self.reaction_environment.iri, reaction_pressure.iri])
                if 'pressure.control.type' in item: 
                    self.instance_dict['details'].append([reaction_pressure.iri, item['pressure.control.type']])
                if 'pressure.atmosphere.type' in item: 
                    atmosphere_type = atmosphere_mapping.get(item['pressure.atmosphere.type'])
                    reaction_atmosphere = self.mds.ReactionAtmosphere(f"ReactionAtmosphere#{self.reaction_id}_{atmosphere_type}")
                    self.instance_dict['type'].append([reaction_atmosphere.iri, reaction_atmosphere.is_instance_of[0].iri])
                    self.instance_dict['is made of'].append([self.reaction_environment.iri, reaction_atmosphere.iri])

            if item['stirring'] == True and 'stirring.type' in item :
                stirring_condition = self.obo.CHMO_0002774(f'StirringProcess#{self.reaction_id}')
                self.instance_dict['type'].append([stirring_condition.iri, stirring_condition.is_instance_of[0].iri])
                self.instance_dict['has process part'].append([self.chemical_reaction.iri, stirring_condition.iri])
                if 'stirring.rate.type' in item:
                    stirring_rate = self.mds.StirringRate(f'StirringRate#{self.reaction_id}')
                    self.instance_dict['type'].append([stirring_rate.iri, stirring_rate.is_instance_of[0].iri])
                    self.instance_dict['has occurent part'].append([stirring_condition.iri, stirring_rate.iri])
                    self.instance_dict['details'].append([stirring_rate.iri, item['stirring.rate.type']])
            
            if item['illumination'] == True and 'illumination.type' in item:
                illumination_conditions = self.mds.Illumination(f'IlluminationProcess#{self.reaction_id}')
                self.instance_dict['type'].append([illumination_conditions.iri, illumination_conditions.is_instance_of[0].iri])
                self.instance_dict['has process part'].append([self.chemical_reaction.iri, illumination_conditions.iri])
                if 'illumination.peak_wavelength.value' in item and 'illumination.peak_wavelength.units' in item:
                    peak_wavelength = self.mds.PeakWavelength(f'PeakWavelength#{self.reaction_id}')
                    self.instance_dict['type'].append([peak_wavelength.iri, peak_wavelength.is_instance_of[0].iri])
                    self.instance_dict['has occurent part'].append([illumination_conditions.iri, peak_wavelength.iri])
                    self.instance_dict['has decimal value'].append([peak_wavelength.iri, item['illumination.peak_wavelength.value']])
                    self.instance_dict['uses measurement unit'].append([peak_wavelength.iri, item['illumination.peak_wavelength.units']])
                #if 'illumination.distance_to_vessel.value' in item and 'illumination.distance_to_vessel.units' in item: 

            if item['electrochemistry'] == True and 'electrochemistry.type' in item: 
                electrochemistry_condition = self.mds.ElectrochemicalReaction(f'ElectrochemicalReaction#{self.reaction_id}')
                self.instance_dict['type'].append([electrochemistry_condition.iri, electrochemistry_condition.is_instance_of[0].iri])
                self.instance_dict['has process part'].append([self.chemical_reaction.iri, electrochemistry_condition.iri])
                if 'electrochemistry.details' in item: 
                    self.instance_dict['details'].append([electrochemistry_condition.iri, item['electrochemistry.details']])
                if 'electrochemistry.current.value' in item and 'electrochemistry.current.units' in item:
                    self.instance_dict['has decimal value'].append([electrochemistry_condition.iri, item['electrochemistry.current.value']])
                    self.instance_dict['uses measurement unit'].append([electrochemistry_condition.iri, item['electrochemistry.current.units']])
                if 'electrochemistry.voltage.value' in item and 'electrochemistry.voltage.units' in item:
                    self.instance_dict['has decimal value'].append([electrochemistry_condition.iri, item['electrochemistry.voltage.value']])
                    self.instance_dict['uses measurement unit'].append([electrochemistry_condition.iri, item['electrochemistry.voltage.units']])
                if 'electrochemistry.anode_material' in item: 
                    self.instance_dict['has text value'].append([electrochemistry_condition.iri, item['electrochemistry.anode_material']])
                if 'electrochemistry.cathode_material' in item: 
                    self.instance_dict['has text value'].append([electrochemistry_condition.iri, item['electrochemistry.cathode_material']])
                if 'electrochemistry.electrode_separation.value' in item and 'electrochemistry.electrode_separation.units' in item: 
                    electrode_separation = cco.ont00000738(f"ElectrodeSeparation#{self.reaction_id}")
                    self.instance_dict['type'].append([electrode_separation.iri, electrode_separation.is_instance_of[0].iri])
                    #self.instance_dict['inheres in'].append([tube_diameter.iri, flow_tube.iri])
                    #self.instance_dict['has decimal value'].append([tube_diameter.iri, item['flow.tubing.diameter.value']])
                    #self.instance_dict['uses measurement unit'].append([tube_diameter.iri, item['flow.tubing.diamter.units']])
                
                    self.instance_dict['has decimal value'].append([electrochemistry_condition.iri, item['electrochemistry.electrode_separation.value']])
                    self.instance_dict['uses measurement unit'].append([electrochemistry_condition.iri, item['electrochemistry.electrode_separation.units']])
                #if 'electrochemistry.cell.type' in item:
                #    self.instance_dict[]
            if item['flow'] == True and 'flow.type' in item: 
                flow_condition = self.mds.ContinuousFlow(f'FlowConditions#{self.reaction_id}')
                self.instance_dict['type'].append([flow_condition.iri, flow_condition.is_instance_of[0].iri])
                self.instance_dict['has process part'].append([self.chemical_reaction.iri, flow_condition.iri])
                if 'flow.details' in item: 
                    self.instance_dict['details'].append([flow_condition.iri, item['flow.details']])
                if 'flow.pump_type' in item: 
                    self.instance_dict['details'].append([flow_condition.iri, item['flow.pump_type']])
                if 'flow.tubing.type' in item: 
                    tube_mapping = {
                        'UNSPECIFIED': None,
                        'CUSTOM': None,
                        'STEEL': 'SteelTube',
                        'COPPER': 'CopperTube',
                        'PFA': 'PerfluoroalkoxyTube',
                        'FEP': 'FluorinatedEthylenePropyleneTube',
                        'TEFLONAF': 'TeflonAmorphousFluoropolymersTube', 
                        'PTFE': 'PolytetrafluoroethyleneTube',
                        'GLASS': 'GlassTube',
                        'QUARTZ': 'QuartzTube',
                        'SILICON': 'SiliconTube',
                        'PDMS': 'PolydimethylsiloxaneTube',
                    }
                    flow_tube = cco.ont00000581(f"{tube_mapping[item['flow.tubing.type']]}#{self.reaction_id}")
                    self.instance_dict['type'].append([flow_tube.iri, flow_tube.is_instance_of[0].iri])
                    self.instance_dict['details'].append([flow_tube.iri, item['flow.tubing.details']]) if 'flow.tubing.details' in item else None
                    if 'flow.tubing.diameter' in item and 'flow.tubing.diameter.value' in item and 'flow.tubing.diamter.units' in item:
                        tube_diameter = cco.ont00000738(f"TubeDiameter#{self.reaction_id}")
                        self.instance_dict['type'].append([tube_diameter.iri, tube_diameter.is_instance_of[0].iri])
                        self.instance_dict['inheres in'].append([tube_diameter.iri, flow_tube.iri])
                        self.instance_dict['has decimal value'].append([tube_diameter.iri, item['flow.tubing.diameter.value']])
                        self.instance_dict['uses measurement unit'].append([tube_diameter.iri, item['flow.tubing.diamter.units']])
    
    def _process_reaction_notes(self):
        for item in self.reaction_notes:
            if 'is_heterogeneous' in item: 
                self.instance_dict['is heterogeneous'].append([self.chemical_reaction.iri, item['is_heterogeneous']])
            if 'forms_precipitate' in item: 
                self.instance_dict['forms precipitate'].append([self.chemical_reaction.iri, item['forms_precipitate']])
            if 'is_exothermic' in item: 
                self.instance_dict['is exothermic'].append([self.chemical_reaction.iri, item['is_exothermic']])
            if 'offgasses' in item: 
                self.instance_dict['off gasses'].append([self.chemical_reaction.iri, item['offgasses']])
            if 'is_sensitive_to_moisture' in item: 
                self.instance_dict['is sensitive to moisture'].append([self.chemical_reaction.iri, item['is_sensitive_to_moisture']])
            if 'is_sensitive_to_oxygen' in item: 
                self.instance_dict['is sensitive to oxygen'].append([self.chemical_reaction.iri, item['is_sensitive_to_oxygen']])
            if 'is_sensitive_to_light' in item: 
                self.instance_dict['is sensitive to light'].append([self.chemical_reaction.iri, item['is_sensitive_to_light']])
            if 'safety_notes' in item: 
                self.instance_dict['safety notes'].append([self.chemical_reaction.iri, item['safety_notes']])
            if 'procedure_details' in item: 
                self.instance_dict['procedure details'].append([self.chemical_reaction.iri, item['procedure_details']])

    def _process_reaction_workups(self):
        """ Process reaction workups that take in the crude product from the chemical reaction to produce purified products """ 
        workup_mapping = {
            'UNSPECIFIED': None,
            'CUSTOM': None,
            'ADDITION': 'Addition',
            'ALIQUOT': 'Aliquot',
            'TEMPERATURE': 'Temperature',
            'CONCENTRATION': 'Concentration',
            'EXTRACTION': 'Extraction',
            'FILTRATION': 'Filtration',
            'WASH': 'Wash',
            'DRY_IN_VACUUM': 'DryInVacuum',
            'DRY_WITH_MATERIAL': 'DryWithMaterial',
            'FLASH_CHROMATOGRAPHY': 'FlashChromatography',
            'OTHER_CHROMATOGRAPHY': 'OtherChromatography',
            'SCAVENGING': 'Scavenging',
            'WAIT': 'Wait',
            'STIRRING': 'Stirring',
            'PH_ADJUST': 'pHAdjust',
            'DISSOLUTION': 'Dissolution',
            'DISTILLATION': 'Distillation'
        }
        for item in self.reaction_workups: 
            workup_type = workup_mapping.get(item['type'])
            workup_instance = self.mds.ReactionWorkup(f'ReactionWorkup#{self.reaction_id}_{workup_type}_{item["Index"]}')
            if item['details'] is not None: 
                self.instance_dict['details'].append([workup_instance.iri, item['details']])
            if item['type'] == 'WAIT' and 'duration.value' in item and 'duration.units' in item: 
                duration_instance = self.obo.BFO_0000202(f'WorkupDuration#{self.reaction_id}_{item["Index"]}')
                self.instance_dict['type'].append([duration_instance.iri, duration_instance.is_instance_of[0].iri])
                self.instance_dict['occupies temporal region'].append([workup_instance.iri, duration_instance.iri])
                self.instance_dict['has datetime value'].append([duration_instance.iri, item['duration.value']])
                self.instance_dict['uses measurement unit'].append([duration_instance.iri, self.unit_mapping[item['duration.units']].iri])
            self.instance_dict['is input of'].append([self.crude_product.iri, workup_instance.iri])
            self.instance_dict['precedes'].append([self.chemical_reaction.iri, workup_instance.iri])
            if (item['type'] == 'ADDITION' or item['type'] ==  'DRY_WITH_MATERIAL' or item['type'] ==  'WASH' or item['type'] == 'SCAVENGING') and item['input'] == True:
                self._extract_components(item, workup_instance, context='workup')
            if (item['type'] == 'ALIQUOT') and 'amount.value' in item and 'amount.units' in item and 'workup.amount.type' in item: 
                self._extract_component_amount(item, workup_instance, context='aliquot')
            if (item['type'] == 'TEMPERATURE' or item['type'] == 'DRY_IN_VACUUM' or item['type'] == 'DISTILLATION') and 'workup.temperature.setpoint.value' in item and 'workup.temperature.setpoint.units' in item : 
                self._process_temperature(item, context='workup')
            if (item['type'] == 'EXTRACTION' or item['type'] == 'FILTRATION') and 'keep_phase' in item:
                self.instance_dict['keep phase'].append([workup_instance.iri, item['keep_phase']])
            if item['type'] == 'PH_ADJUST' and 'target_ph' in item: 
                self.instance_dict['details'].append([workup_instance.iri, item['target_ph']])
            if item['type'] == 'STIRRING' and item['stirring'] == True: 
                stirring_condition = self.obo.CHMO_0002774(f'StirringProcess#{self.reaction_id}_Workup_{item["Index"]}')
                self.instance_dict['type'].append([stirring_condition.iri, stirring_condition.is_instance_of[0].iri])
                self.instance_dict['has process part'].append([workup_instance.iri, stirring_condition.iri])
                if 'stirring.rate.type' in item:
                    stirring_rate = self.mds.StirringRate(f'StirringRate#{self.reaction_id}')
                    self.instance_dict['type'].append([stirring_rate.iri, stirring_rate.is_instance_of[0].iri])
                    self.instance_dict['has occurent part'].append([stirring_condition.iri, stirring_rate.iri])
                    self.instance_dict['details'].append([stirring_rate.iri, item['stirring.rate.type']])
            if (item['type'] == 'FLASH_CHROMATOGRAPHY' or item['type'] == 'OTHER_CHROMATOGRAPHY') and 'is_automated' in item:
                self.instance_dict['is automated'].append([workup_instance.iri, item['is_automated']])
                

    def _extract_product_measurement(self, outcome_list, product_dict):
        """ Process measurements of the products of the chemical reaction """
        pattern = rf'^products\[(\d+)\]'
        index_list = self._extract_index_set(outcome_list, pattern)
        measurement_dict = {}

        if not index_list:
            return self, measurement_dict
        
        for i in index_list: 
            product_prefix = f"products[{i}].measurements"
            pattern = rf'^products\[{i}\]\.measurements\[(\d+)\]\.(.+)$'
            ind_list = self._extract_index_set(outcome_list, pattern)
            
            for j in ind_list: 
                if i not in product_dict: 
                    continue
                if f"{product_prefix}[{j}].type" in outcome_list and f"{product_prefix}[{j}].analysis_key" in outcome_list: 
                    product_measurement = self.mds.AnalyticalResult(f'ProductMeasurement#{self.reaction_id}_{outcome_list[f"{product_prefix}[{j}].analysis_key"]}_Product_{i}_Measurement_{j}')
                    measurement_dict[outcome_list[f"{product_prefix}[{j}].analysis_key"]] = product_measurement
                
                elif f"{product_prefix}[{j}].type" in outcome_list: 
                    product_measurement = self.mds.AnalyticalResult(f'ProductMeasurement#{self.reaction_id}_Product_{i}_Measurement_{j}')
                else: 
                    product_measurement = None
                if product_measurement: 
                    self.instance_dict['type'].append([product_measurement.iri, product_measurement.is_instance_of[0].iri])
                    self.instance_dict['is about'].append([product_measurement.iri, product_dict[i].iri])

                    if f"{product_prefix}[{j}].percentage.value" in outcome_list: 
                        self.instance_dict['has decimal value'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].percentage.value"]])
                        self.instance_dict['uses measurement unit'].append([product_measurement.iri, self.unit_mapping['PERCENTAGE'].iri])
                    
                    if f"{product_prefix}[{j}].float_value.value" in outcome_list: 
                        self.instance_dict['has decimal value'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].float_value.value"]])

                    if f"{product_prefix}[{j}].string_value" in outcome_list: 
                        self.instance_dict['has text value'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].string_value"]])

                    if outcome_list[f"{product_prefix}[{j}].type"] == 'AMOUNT' and f"{product_prefix}[{j}].amount.value" in outcome_list and f"{product_prefix}[{j}].amount.units" in outcome_list: 
                        amount_added = self.cco.ont00000768('Mass#' + outcome_list['reactionID'] + '_' + f"_Product_{str(i)}_Measurement_{str(j)}")
                        self.instance_dict['type'].append([amount_added.iri, amount_added.is_instance_of[0].iri])
                        self.instance_dict['inheres in'].append([amount_added.iri, product_dict[i].iri])
                        self.instance_dict['is about'].append([product_measurement.iri, amount_added.iri])
                        self.instance_dict['has decimal value'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].amount.value"]])
                        self.instance_dict['uses measurement unit'].append([product_measurement.iri, self.unit_mapping[outcome_list[f"{product_prefix}[{j}].amount.units"]].iri])

                    if f"{product_prefix}[{j}].retention_time.value" in outcome_list and f"{product_prefix}[{j}].retention_time.units" in outcome_list: 
                        retention_time = self.mds.RetentionTime('RetentionTime#' + outcome_list['reactionID'] + f"_Product_{str(i)}_Measurement_{str(j)}")
                        self.instance_dict['type'].append([retention_time.iri, retention_time.is_instance_of[0].iri])
                        self.instance_dict['is about'].append([retention_time.iri, product_dict[i].iri])
                        self.instance_dict['is about'].append([product_measurement.iri, retention_time.iri])
                        self.instance_dict['has datetime value'].append([retention_time.iri, outcome_list[f"{product_prefix}[{j}].retention_time.value"]])
                        self.instance_dict['uses measurement unit'].append([retention_time.iri, self.unit_mapping[outcome_list[f"{product_prefix}[{j}].retention_time.units"]].iri])
                    
                    if f"{product_prefix}[{j}].uses_internal_standard" in outcome_list: 
                        self.instance_dict['uses internal standard'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].uses_internal_standard"]])
                    if f"{product_prefix}[{j}].uses_authentic_standard" in outcome_list:
                        self.instance_dict['uses authentic standard'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].uses_authentic_standard"]])
                    if f"{product_prefix}[{j}].is_normalized" in outcome_list:
                        self.instance_dict['is normalized'].append([product_measurement.iri, outcome_list[f"{product_prefix}[{j}].is_normalized"]])
        return self, measurement_dict
    
    def _process_reaction_outcomes(self):
        """ Method to process reaction outcomes of a chemical reaction: products, analyses, and measurements of products (i.e. result from analyses) """

        for item in self.reaction_outcomes: 
            if item.get('reaction_time') == True and 'reaction_time.value'in item and 'reaction_time.units' in item:
                reaction_time = self.mds.ReactionTime(f'ReactionTime#{self.reaction_id}')
                self.instance_dict['type'].append([reaction_time.iri, reaction_time.is_instance_of[0].iri])
                self.instance_dict['occupies temporal region'].append([self.chemical_reaction.iri, reaction_time.iri])
                self.instance_dict['has datetime value'].append([reaction_time.iri, item['reaction_time.value']])
                self.instance_dict['uses measurement unit'].append([reaction_time.iri, self.unit_mapping[item['reaction_time.units']].iri])
            
            try: 
                _, product_dict = self._extract_components(item, context='product')
            except Exception as e: 
                logger.error(f"Failed to process product components into RDF triples for reaction {self.reaction_id}: {e}")
                pass
            
            try:
                _, measurement_dict = self._extract_product_measurement(item, product_dict)
            except Exception as e:
                logger.error(f"Failed to process product measurements into RDF triples for reaction {self.reaction_id}: {e}")
                pass

            if item.get('analyses') == True: 
                pattern = rf'^analyses\["([^"]*)"\]'  
                analyses_set = self._extract_index_set(item, pattern)

                for analysis in analyses_set:
                    if f'analyses["{analysis}"].type' in item: 
                        analysis_type = self.mds.AnalyticalTechnique(f'AnalyticalTechnique#{self.reaction_id}_{analysis}')
                        self.instance_dict['type'].append([analysis_type.iri, analysis_type.is_instance_of[0].iri])
                        self.instance_dict['details'].append([analysis_type.iri, item[f'analyses["{analysis}"].type']])
                        #self.instance_dict['has output'].append([analysis_type.iri, measurement_dict[analysis]])
                    
                    if f'analyses["{analysis}"].details' in item: 
                        self.instance_dict['details'].append([analysis_type.iri, item[f'analyses["{analysis}"].details']])
                    if f'analyses["{analysis}"].is_of_isolated_species' in item: 
                        self.instance_dict['is of isolated species'].append([analysis_type.iri, item[f'analyses["{analysis}"].is_of_isolated_species']])
                    #if f'analyses["{analysis}"].chmo_id' in item:
                    #    self.instance_dict['CHMO ID'].append([analysis_type.iri, item[f'analyses["{analysis}"].chmo_id']])
                    if f'analyses["{analysis}"].instrument_manufacturer' in item:
                        self.instance_dict['has text value'].append([analysis_type.iri, item[f'analyses["{analysis}"].instrument_manufacturer']])
                    
                    pattern = rf'^analyses\["{re.escape(analysis)}"\]\.data\["([^"]*)"\]'
                    data_set = self._extract_index_set(item, pattern)
                    
                    if data_set: 
                        for data in data_set: 
                            if f'analyses["{analysis}"].data["{data}"].float_value.value' in item: 
                                self.instance_dict['has decimal value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].float_value.value']])
                            if f'analyses["{analysis}"].data["{data}"].integer_value' in item: 
                                self.instance_dict['has integer value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].integer_value']])
                            if f'analyses["{analysis}"].data["{data}"].bytes_value' in item: 
                                self.instance_dict['has text value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].bytes_value']])
                            if f'analyses["{analysis}"].data["{data}"].string_value' in item: 
                                self.instance_dict['has text value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].string_value']])
                            if f'analyses["{analysis}"].data["{data}"].description' in item: 
                                self.instance_dict['has text value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].description']])
                            if f'analyses["{analysis}"].data["{data}"].format' in item: 
                                self.instance_dict['has text value'].append([analysis_type.iri, item[f'analyses["{analysis}"].data["{data}"].format']])
    
    def _process_reaction_setup(self):
        for item in self.reaction_setup: 
            reaction_setup = self.mds.ReactionSetup(f'ReactionSetup#{self.reaction_id}')
            self.instance_dict['precedes'].append([reaction_setup.iri, self.chemical_reaction.iri])

    def generate_data_graph(self, dataset_id, save_file_path): 
        
        # Initialize graph, binding to the corresponding namespaces
        self.graph = Graph()
        self.graph.bind('afe', AFE)
        self.graph.bind('afr', AFR)
        self.graph.bind('afrl', AFRL)
        self.graph.bind('afq', AFQ)
        self.graph.bind('afm', AFM)
        self.graph.bind('obo', OBO)
        self.graph.bind('cco', CCO)
        self.graph.bind('mds', MDS)
        self.graph.bind('ncit', NCIT)
        self.graph.bind('qudt', QUDT)
        self.graph.bind('unit', UNIT)
        
        try: 
            for key in self.prop_metadata_dict.keys(): 
                prop_metadata = self.prop_metadata_dict.get(key)
                for item in self.instance_dict[key]: 
                    if prop_metadata[1] == "Object Property": 
                        subj_uri = URIRef(item[0])
                        obj_uri = URIRef(item[1])
                        pred_uri = URIRef(prop_metadata[0])
                        self.graph.add((subj_uri, pred_uri, obj_uri))
                    elif prop_metadata[1] == "Datatype Property":
                        subj_uri = URIRef(item[0])
                        pred_uri = URIRef(prop_metadata[0])
                        self.graph.add((subj_uri, pred_uri, Literal(item[1])))
            
            if self.fmt == "turtle": 
                new_path = os.path.join(save_file_path, f"mds_dataset-{dataset_id}_reaction-{self.reaction_id}.ttl")
                self.graph.serialize(new_path, format=self.fmt)
                logger.info(f"Saved turtle file: {new_path}")
            elif self.fmt == "json-ld": 
                new_path = os.path.join(save_file_path, f"mds_dataset-{dataset_id}_reaction-{self.reaction_id}.jsonld")
                self.graph.serialize(new_path, format=self.fmt, auto_compact=True)
                logger.info(f"Saved JSON-LD file: {new_path}")
            else: 
                raise ValueError(f"unsupported format: {self.fmt}")
        
        except Exception as e: 
            logger.error(f"Failed to save graph for reaction {self.reaction_id}: {e}")
            return None


        return self


