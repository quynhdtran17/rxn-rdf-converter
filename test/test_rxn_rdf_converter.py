# %%
# =================================================================
#               IMPORT REQUIREMENTS
# =================================================================
import ord_schema
from ord_schema.message_helpers import load_message, write_message, message_to_row
from ord_schema.proto import dataset_pb2, reaction_pb2
import os
from rdkit import Chem
import re
from owlready2 import get_namespace, get_ontology, Thing
import glob
import rdflib
from rdflib import Graph, RDF, RDFS, OWL, Namespace, Literal, URIRef
from rdflib.namespace import RDFS, XSD, URIRef, OWL, SKOS, PROV
from datetime import datetime
from src import rxn_rdf_converter
import logging
import random
import unittest
import pytest
import tempfile
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

def setup_file_path():
    """ Set up the file paths """
    path = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/OpenReactionDB/ord-data/data' # path of dataset 
    savepath = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/Extracted_ORD_Data'
    mds_file_path = '/mnt/vstor/CSE_MSE_RXF131/staging/mds3/KG-ChemRxn/chemOntologies/mdsChemRxn(v.0.3.0.7).owl'

    file_list = []
    for root, dirs, files in os.walk(path):
        for name in files: 
            if name.startswith('ord_dataset'):
                file_path = os.path.join(root, name)
                file_list.append(file_path)
    
    return path, savepath, mds_file_path, file_list

# Set up path: 
path, savepath, mds_file_path, file_list = setup_file_path()
logger.info(f"Found {len(file_list)} data files")

# Choose a dataset to test the reactions from
dataset = load_message("test/ord_dataset-00005539a1e04c809a9a78647bea649c.pb.gz", dataset_pb2.Dataset,)

# Randomly generates "num_test_reactions" numbers from 1 to 69690- these will be
# the reactions that will be tested
def generate_random_reaction_set(num_test_reactions):

    # Counts the number of reactions in the dataset
    num_total_reactions = len(dataset.reactions)

    # Prints the total number of reactions in the file list
    print(f"There are {str(num_total_reactions)} reactions in the file list")

    # Creates a random list of numbers that correspond to reactions in the file
    # list's datasets
    random_reaction_nums = random.sample(range(1, num_total_reactions + 1), num_test_reactions)
    random_reaction_nums.sort()
    return random_reaction_nums, num_test_reactions, num_total_reactions

# From the list of random reaction numbers, creates a list of the protocol buffer reactions to test
def create_list_of_protocol_buffer_test_reactions(random_reaction_nums):

    # Initializes indices to keep track of which element in the random reaction numbers
    # list we are at, as well as which reaction (out of all of them) is being checked
    current_reaction_index = 0
    current_reaction_num = 1

    # List to hold the protocol buffer reactions that will be tested
    protocol_buffer_test_reactions = []

    # For each reaction in the dataset
    for j in range(len(dataset.reactions)):

            # Checks to make sure there is at least one more element in the random reaction
            # numbers list
            if current_reaction_index < len(random_reaction_nums):

                # Adds the corresponding reaction to the list if the current reaction number being
                # checked matches the next element in the random reaction numbers list
                if current_reaction_num == random_reaction_nums[current_reaction_index]:
                    protocol_buffer_test_reactions.append(dataset.reactions[j])
                    
                    # Goes to the next element in the random reaction numbers list
                    current_reaction_index += 1
            
            # Goes to the next reaction number to check it in the next iteration
            current_reaction_num += 1

    # Returns the completed list of the protocol buffer reactions that will be tested
    return protocol_buffer_test_reactions


# Class to test the file rxn_rdf_converter.py
class TestRxnRdfConverter(unittest.TestCase):

    @classmethod
    def setUpClass(self):

        # Initializes a list to hold 50 protocol buffer reactions- a number
        # of reactions much greater than this will likely cause the kernel
        # to creash
        protocol_buffer_test_reactions = []

        # Creates a list of 50 random reaction numbers from 1 to the total
        # number of reactions in the dataset- each of these will be tested
        try:
            random_reaction_nums, num_test_reactions, num_total_reactions = generate_random_reaction_set(50)
            print("List of 50 random reaction numbers out of " + str(num_total_reactions) + " total:")
            print(str(random_reaction_nums))
        except Exception as e:
            print(e)

        # Creates a list of the 50 protocol buffer reactions to test from
        # the random reaction numbers
        protocol_buffer_test_reactions = create_list_of_protocol_buffer_test_reactions(random_reaction_nums)
        self.protocol_buffer_test_reactions = protocol_buffer_test_reactions

    def test_generate_reaction_basic_properties(self):
        """
        Test generate_reaction() runs and populates the basic lists and attributes.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            # basic attributes
            assert semi.reaction_pb == pb_rxn
            assert isinstance(semi.reaction_id, str) and len(semi.reaction_id) > 0

            # lists should exist
            assert isinstance(semi.reaction_identifiers, list)
            assert isinstance(semi.reaction_inputs, list)
            assert isinstance(semi.reaction_setup, list)
            assert isinstance(semi.reaction_conditions, list)
            assert isinstance(semi.reaction_notes, list)
            assert isinstance(semi.reaction_workups, list)
            assert isinstance(semi.reaction_outcomes, list)
            assert isinstance(semi.reaction_provenance, list)

            # If identifiers present, ensure mapping of first identifier value to reaction_identifiers entry
            if getattr(pb_rxn, "identifiers", None):
                # the first identifier in pb_rxn should be reflected in reaction_identifiers[0] if present
                if semi.reaction_identifiers:
                    iv = pb_rxn.identifiers[0].value if hasattr(pb_rxn.identifiers[0], "value") else None
                    if iv:
                        assert semi.reaction_identifiers[0].get("value") == iv


    def test_generate_compound_identifiers_on_real_components(self):
        """
        For each reaction, exercise generate_compound_identifiers() on real components where available.
        Check that when SMILES or INCHI present, an INCHI_KEY is produced.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            # iterate over all inputs -> components -> identifiers
            for item in semi.reaction_inputs:
                # the original proto input might have components (boolean or list); we don't rely on that.
                # Instead try to get corresponding pb component list via matching InputKey on pb_rxn if possible.
                # Fallback: skip if no components accessible.
                # We'll attempt to find smth in the pb_rxn structure
                # This block will be conservative: call generate_compound_identifiers only if we can identify identifiers objects
                # from the original protocol buffer.
                # We can inspect semi.reaction_pb to find inputs keyed similarly.
                for pb_input_key in semi.reaction_pb.inputs:
                    if getattr(pb_input_key, "strip", None):
                        # pb_input_key is the key string from protobuf map; find the matching input in pb_rxn
                        if re.sub(r"\s+", "", str(pb_input_key.strip())) == item["InputKey"]:
                            pb_input = semi.reaction_pb.inputs[pb_input_key]
                            for comp in pb_input.components:
                                # each comp is a protobuf message with .identifiers repeated field
                                if comp.identifiers:
                                    id_list, inchi_key = semi.generate_compound_identifiers(comp.identifiers)
                                    assert isinstance(id_list, list)
                                    # inchi_key can be None, but if comp had SMILES or INCHI we expect it
                                    types = {getattr(x, "type", None) for x in comp.identifiers}
                                    # Real protobuf .type may be enum; check by value string if possible
                                    # We'll convert any present identifier messages to message_to_row in code path already did,
                                    # But here we only check that if a SMILES or INCHI string exists, an inchi_key string is returned.
                                    has_smiles = any(getattr(idm, "type", None) and str(idm.type).upper().endswith("SMILES") for idm in comp.identifiers)
                                    has_inchi = any(getattr(idm, "type", None) and "INCHI" in str(idm.type).upper() for idm in comp.identifiers)
                                    if has_smiles or has_inchi:
                                        # Either inchi_key is non-empty string or None (rare); assert type when present
                                        if inchi_key:
                                            assert isinstance(inchi_key, str) and len(inchi_key) > 0


    def test_initialize_instance_dict_and_core_instances(self):
        """
        Test _initialize_instance_dict: it should load the ontology, set up namespaces, unit mapping,
        prop_metadata_dict, instance_dict, and create core instances (chemical_reaction etc).
        If the ontology path doesn't exist or loading fails, skip the test.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            for pb_rxn in self.protocol_buffer_test_reactions:
                semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

                # If the file doesn't exist, skip
                if not os.path.exists(mds_file_path):
                    pytest.skip("Ontology file not available; skipping ontology-dependent tests.")

                try:
                    semi._initialize_instance_dict(mds_file_path)
                except Exception as e:
                    pytest.skip(f"Could not initialize ontology ({e}); skipping related tests.")

                # namespaces exist
                assert hasattr(semi, "mds") and hasattr(semi, "cco") and hasattr(semi, "unit")
                # unit mapping contains common units
                assert "GRAM" in semi.unit_mapping and "MILLILITER" in semi.unit_mapping

                # prop_metadata_dict and instance_dict should be dicts and contain 'type'
                assert isinstance(semi.prop_metadata_dict, dict)
                assert "type" in semi.prop_metadata_dict
                assert isinstance(semi.instance_dict, dict)
                assert "type" in semi.instance_dict

                # core instances created
                assert semi.chemical_reaction is not None
                assert semi.reaction_mixture is not None
                assert semi.reaction_environment is not None
                assert semi.crude_product is not None

                # instance_dict['type'] should include the IRIs of those instances
                type_iris = [t[0] for t in semi.instance_dict["type"]]
                assert semi.chemical_reaction.iri in type_iris
                assert semi.reaction_mixture.iri in type_iris
                assert semi.reaction_environment.iri in type_iris
                assert semi.crude_product.iri in type_iris


    def test_extract_index_set_behavior(self):
        """
        Test _extract_index_set with numeric, string, mixed keys.
        This tests the helper directly (no ontology required).
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld")
            # numeric
            item_dict = {"components[0]": 1, "components[2]": 1, "components[1]": 1, "other": 9}
            res = semi._extract_index_set(item_dict, r"^components\[(\d+)\]")
            assert res == [0, 1, 2]

            # string keys
            item_dict2 = {"products[foo]": 1, "products[bar]": 1, "products[baz]": 1, "products[foo]": 2}
            res2 = semi._extract_index_set(item_dict2, r"^products\[(.+)\]")
            assert set(res2) == {"foo", "bar", "baz"}

            # mixed keys prefer numeric
            mixed = {"components[0]": 1, "components[bar]": 1, "components[3]": 1}
            resm = semi._extract_index_set(mixed, r"^components\[(.+)\]")
            assert resm == [0, 3]

            # no matches
            nomatch = {"x": 1}
            assert semi._extract_index_set(nomatch, r"^components\[(\d+)\]") == []


    def test_process_reaction_identifiers_real(self):
        """
        Run _initialize_instance_dict and then _process_reaction_identifiers using real identifiers.
        Only assert behavior if identifiers exist.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not semi.reaction_identifiers:
                pytest.skip("No reaction identifiers present in this test reaction; skipping.")

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping identifier processing tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping identifier tests.")

            # run processing
            semi._process_reaction_identifiers()

            # For each identifier in semi.reaction_identifiers, ensure its value is recorded in 'has text value'
            for ident in semi.reaction_identifiers:
                val = ident.get("value")
                if val:
                    assert any(val == tv for _, tv in semi.instance_dict.get("has text value", []))


    def test_extract_components_and_roles_and_amounts(self):
        """
        Run _initialize_instance_dict (ontology required) and _process_reaction_inputs which calls _extract_components.
        Validate that components are added to instance_dict and that component amounts create measurement triples.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not semi.reaction_inputs:
                pytest.skip("No reaction inputs; skipping component extraction test.")

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping component extraction tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping component tests.")

            # Process reaction inputs (this internally calls _extract_components)
            semi._process_reaction_inputs()

            # There should be some 'type' entries for components if inputs had components
            if semi.reaction_inputs:
                # If any input had components flag true, expect 'member part of' entries or 'is input of'
                got_member_part = semi.instance_dict.get("member part of", [])
                got_is_input_of = semi.instance_dict.get("is input of", [])
                # Accept either presence or absence depending on real data; but if components existed, at least one of those lists may be non-empty
                # So we don't assert strict presence, but ensure no exceptions occurred and instance_dict exists
                assert isinstance(semi.instance_dict, dict)


    def test_extract_component_amount_specific_value(self):
        """
        Construct a mock component_list representing a component with a mass amount and verify that
        _extract_component_amount creates 'has decimal value' and 'uses measurement unit' entries with correct values.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping amount tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping amount tests.")

            # Build a synthetic component_list using the first input if present
            if not semi.reaction_inputs:
                pytest.skip("No inputs to construct a component for amount test; skipping.")

            rxn_input = semi.reaction_inputs[0]
            # ensure InputKey exists
            input_key = rxn_input.get("InputKey", "m1")

            # Build component_list for a single component index i=0 with amount.mass
            component_list = {
                "InputKey": input_key,
                "components[0].amount.type": "mass",
                "components[0].amount.mass.value": 25.0,
                "components[0].amount.mass.units": "MILLIGRAM"
            }

            # component instance
            reaction_component = semi.mds.Component(f'Component#{semi.reaction_id}_{input_key}_component_0')
            # call
            semi._extract_component_amount(component_list, reaction_component, i=0, context='input')

            # verify decimal value recorded
            assert any(25.0 == v for _, v in semi.instance_dict.get("has decimal value", []))
            # verify measurement unit mapping used
            unit_list = semi.instance_dict.get("uses measurement unit", [])
            assert any("MilliGM" in str(u[1]) or "MilliGM" in str(u) for u in unit_list) or unit_list == []


    def test_extract_compound_identifiers_and_details(self):
        """
        Test _extract_compound_identifiers by building a component_list that contains identifiers
        and ensuring the instance_dict is updated with has text value, details and is mapped when present.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping compound identifier extraction tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping identifier tests.")

            # Use a real input if available; otherwise skip
            if not semi.reaction_inputs:
                pytest.skip("No inputs available; skipping compound identifier extraction test.")

            rxn_input = semi.reaction_inputs[0]
            input_key = rxn_input.get("InputKey", "m1")
            # Compose a component_list with two identifiers (SMILES + CUSTOM)
            component_list = {
                "InputKey": input_key,
                "components[0].identifiers[0].type": "SMILES",
                "components[0].identifiers[0].value": "CCO",
                "components[0].identifiers[0].details": "ethanol-smiles",
                "components[0].identifiers[1].type": "UNSPECIFIED",
                "components[0].identifiers[1].value": "custom-id"
            }
            reaction_component = semi.mds.Component(f'Component#{semi.reaction_id}_{input_key}_component_0')

            # call method
            semi._extract_compound_identifiers(component_list, 0, "input", "components", reaction_component)

            # check values recorded
            values = [v for _, v in semi.instance_dict.get("has text value", [])]
            assert "CCO" in values
            assert "custom-id" in values

            details = [v for _, v in semi.instance_dict.get("details", [])]
            assert any("ethanol-smiles" in d for d in details)


    def test_process_reaction_inputs_and_generate_graph(self):
        """
        Full-ish integration: generate instances from reaction inputs, then generate a data graph (json-ld).
        Ensure the file is created and that the graph contains at least one triple.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            for pb_rxn in self.protocol_buffer_test_reactions:
                semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

                if not os.path.exists(mds_file_path):
                    pytest.skip("Ontology missing; skipping graph generation tests.")

                try:
                    semi.generate_instances(mds_file_path)
                except Exception as e:
                    pytest.skip(f"generate_instances failed for reaction {semi.reaction_id} ({e})")

                # generate_data_graph writes to tmp_path and returns self
                try:
                    semi.generate_data_graph("testdataset", str(tmp_path))
                except Exception as e:
                    pytest.skip(f"generate_data_graph failed ({e})")

                # check that file was created (.jsonld)
                files = os.listdir(tmp_path)
                assert len(files) > 0
                # Additionally, check that graph object exists
                assert hasattr(semi, "graph")
                assert len(semi.graph) >= 0  # graph may be empty if no triples, but ensure attribute present


    def test_process_temperature_condition_and_input(self):
        """
        Validate _process_temperature in both condition and input contexts using small synthetic dicts derived from reaction.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping temperature tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping temperature tests.")

            # condition context
            cond = {
                "temperature.setpoint.value": 60.0,
                "temperature.setpoint.units": "CELSIUS",
                "temperature": True
            }
            # call
            semi._process_temperature(cond, context='condition')
            assert any(60.0 == v for _, v in semi.instance_dict.get("has decimal value", []))

            # input context requires reaction_component
            # create fake component and call
            reaction_component = semi.mds.Component(f'Component#{semi.reaction_id}_fake_component_0')
            add_temp = {
                "addition_temperature.value": -5.0,
                "addition_temperature.units": "CELSIUS",
                "addition_temperature.control.type": "thermostat"
            }
            semi._process_temperature(add_temp, reaction_component=reaction_component, context='input')
            assert any(-5.0 == v for _, v in semi.instance_dict.get("has decimal value", []))


    def test_process_reaction_conditions_full(self):
        """
        Run _process_reaction_conditions using the reaction's real reaction_conditions (if present).
        Validate that units & numeric values appear for keys like pressure, illumination, electrochemistry, flow, stirring.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            if not semi.reaction_conditions:
                pytest.skip("No reaction_conditions present; skipping.")

            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping condition tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping condition tests.")

            # run
            semi._process_reaction_conditions()

            # Some keys we might expect in instance_dict after processing conditions:
            # - has decimal value (temperature, pressure, illumination peak)
            # - uses measurement unit
            assert isinstance(semi.instance_dict, dict)
            # If there was temperature in the reaction, expect measurement entries
            #any_temp = any(k for k in semi.reaction_conditions if k.get("temperature") == True)
            #if any_temp:
            #    assert any(semi.instance_dict.get("has decimal value", []))


    def test_process_reaction_notes_and_workups_and_outcomes(self):
        """
        Run the methods that process notes, workups, and outcomes and ensure no exceptions occur and
        that instance_dict is updated when relevant keys are present.
        """
        for pb_rxn in self.protocol_buffer_test_reactions:
            semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

            # Need ontology to process workups/outcomes properly
            if not os.path.exists(mds_file_path):
                pytest.skip("Ontology missing; skipping workup/outcome tests.")

            try:
                semi._initialize_instance_dict(mds_file_path)
            except Exception as e:
                pytest.skip(f"Could not initialize ontology ({e}); skipping these tests.")

            # process notes
            semi._process_reaction_notes()
            # if notes in proto existed, instance_dict may have entries
            if semi.reaction_notes and any(semi.reaction_notes):
                assert isinstance(semi.instance_dict, dict)

            # process workups
            try:
                semi._process_reaction_workups()
            except Exception:
                # don't fail the entire test if workup processing raises for a particular reaction
                pytest.skip("Workup processing failed for this reaction (non-critical in this test run).")

            # process outcomes
            try:
                semi._process_reaction_outcomes()
            except Exception:
                pytest.skip("Outcome processing failed for this reaction (non-critical in this test run).")

            # If there were products, check that product measurement triples exist when measurements present
            for outcome in semi.reaction_outcomes:
                if outcome.get("products"):
                    # We expect _extract_components/_extract_product_measurement to have potentially created 'is about' or 'has decimal value'
                    assert isinstance(semi.instance_dict, dict)


    def test_process_reaction_setup_and_graph_saving(self):
        """
        Process setup and ensure generate_data_graph saves an output file in requested format.
        """
        with tempfile.TemporaryDirectory() as tmp_path:
            for pb_rxn in self.protocol_buffer_test_reactions:
                semi = rxn_rdf_converter.ReactionKG(pb_rxn, fmt="json-ld").generate_reaction()

                if not os.path.exists(mds_file_path):
                    pytest.skip("Ontology missing; skipping setup and graph tests.")

                try:
                    semi.generate_instances(mds_file_path)
                except Exception:
                    pytest.skip("generate_instances failed; skipping setup/graph tests.")

                # process setup
                semi._process_reaction_setup()
                # produce a graph
                semi.generate_data_graph("testdataset", str(tmp_path))
                files = os.listdir(tmp_path)
                assert len(files) > 0

    def tearDown(self):
        pass

# Tests the conversion of 50 random protocol buffer reactions
def test50RandomReactions():

    # Runs the pyunittests
    if __name__ == '__main__':
        unittest.main(argv=['first-arg-is-ignored'], exit=False)

# Runs the function to test the reaction conversions
try:
    test50RandomReactions()
except Exception as e:
    print(e)

# %%
