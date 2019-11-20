#!/usr/bin/env python
#
# Runs R Group Converter for Library Creation
#
# ----------------------------------------------------------

# Imports
# ---------
from os.path import join
from rdkit import Chem
import ruamel.yaml as yaml

# Relative Imports
# ----------------
from cocktail_shaker.validation import RaiseMoleculeError

# Load datasources
# ----------------
def load_datasources():

    """

    Load all the datasources for running this package in local context.

    This might slow down performance later -- we can opt in to load sources of data dependent on the functional group.

    """
    from pathlib import Path
    datasource_location = Path(__file__).absolute().parent
    with open(join(str(datasource_location), "datasources", "R_Groups.yaml")) as stream:
        try:
            global R_GROUP_DATASOURCE
            R_GROUP_DATASOURCE = yaml.safe_load(stream)

            global R_GROUPS
            R_GROUPS = R_GROUP_DATASOURCE['R_Groups']
        except yaml.YAMLError as exc:
            print ("Datasources not loading correctly, Please contact lead developer")
            print(exc)


class LoadCustomLibrary(object):

    """

    Load a custom datasource per the user.

    Convert their SMILES to SMARTS

    """

    __version__ = "1.0.2"


    def __init__(self, molecules):

        # imports
        # -------

        from molvs import Validator
        from validation import MoleculeValidator

        self.molecules = molecules
        self.custom_datasources = False

        validator = MoleculeValidator(self.molecules, smiles=True)


class Cocktail(object):
    """

    This class is used to take in a molecule and replace any R groups with a one of the groups from the R-Group.

    """

    __version_parser__ = 1.0
    __allow_update__ = False

    def __init__(self, molecules):

        # imports
        # -------
        import re
        from molvs import Validator

        # I will allow the user to pass a string but for easier sake down the road
        # I will reimplement it as a list.

        load_datasources()
        rdkit_rendered_molecules = []
        self.markush_structure = False
        for molecule in molecules:
            if 'R1' in molecule:
                molecule = molecule.replace('R1', '[*:1]', 1)
                self.markush_structure = True
            elif 'R' in molecule:
                molecule = molecule.replace('R', '[*:1]', 1)
                self.markush_structure = True

        rdkit_rendered_molecules.append(Chem.MolFromSmiles(molecule))

        self.molecules = rdkit_rendered_molecules
        
        for molecule in self.molecules:
            self.original_smiles = Chem.MolToSmiles(molecule)
            # Validation
            validator_format = '%(asctime)s - %(levelname)s - %(validation)s - %(message)s'
            self.validate = Validator(log_format=validator_format)

    def detect_functional_groups(self):

        """

        Find functional groups that ligand library loader supports

        :return:

        """

        pattern_payload = {}
        load_datasources()

        print ("Detecting Functional Groups...")

        for molecule in self.molecules:
            for functional_group, pattern in R_GROUPS.items():
                for i in range(0, len(pattern)):
                    for key, value in pattern[i].items():
                        smart_pattern = Chem.MolFromSmarts(value[1])
                        if molecule.GetSubstructMatches(smart_pattern,uniquify=False):
                            print ("Found Functional Group: %s | Pattern Count: %s" % (key,
                                                                                       len(molecule.GetSubstructMatches(
                                                                                           smart_pattern,uniquify=False))))
                        pattern_payload[key] = [value[0], value[1]]

        return pattern_payload

    def shake(self, functional_groups=["all"], shape=None):

        """

        Used to swap out molecules based on the patterns found from the detection.

        Arguments:
            self (Object): Cocktail object of the list of molecules

        Return:
            modified_molecules (List): List of the RDKit molecule objects that have had their structures replaced.

        TODO: Do this faster than O(n)^3 as this algorithm is not the most efficient.

        """

        # Run detection first to see and validate what functional groups have been found.

        patterns_found = self.detect_functional_groups()
        print ("Shaking Compound....")
        modified_molecules = []
        if functional_groups[0] == 'all':
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                        smarts_mol = Chem.MolFromSmarts(value[1])
                        for functional_group, pattern in R_GROUPS.items():
                            if functional_group == 'R-Groups':
                                continue
                            else:
                                for i in range(0, len(pattern)):
                                    for r_functional_group, r_data in pattern[0].items():
                                        if r_data[1] == value[1]:
                                            continue
                                        try:
                                            if self.markush_structure:
                                                modified_molecule = Chem.ReplaceSubstructs(molecule, Chem.MolFromSmiles('[*:1]'),
                                                                                           Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                            else:
                                                modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                                          Chem.MolFromSmiles(r_data[0]), replaceAll=True)
                                            modified_molecules.append(modified_molecule[0])
                                        except RaiseMoleculeError:
                                            print ("Molecule Formed is not possible")
        else:
            for molecule in self.molecules:
                for key, value in patterns_found.items():
                    smarts_mol = Chem.MolFromSmarts(value[1])
                    for functional_group, pattern in R_GROUPS.items():
                        if functional_group == 'R-Groups':
                            continue
                        elif functional_group not in functional_groups:
                            continue
                        else:
                            for r_functional_group, r_data in pattern[0].items():
                                # Skip redundacies if the r group is already matched.
                                if r_data[1] == value[1]:
                                    continue
                                try:
                                    if self.markush_structure:
                                        modified_molecule = Chem.ReplaceSubstructs(molecule, Chem.MolFromSmiles('[*:1]'),
                                                                                   Chem.MolFromSmiles(r_data[0]))
                                    else:
                                        modified_molecule = Chem.ReplaceSubstructs(molecule, smarts_mol,
                                                                                   Chem.MolFromSmiles(r_data[0]))
                                    modified_molecules.append(modified_molecule[0])
                                except RaiseMoleculeError:
                                    print ("Molecule Formed is not possible")
        self.modified_molecules = modified_molecules

        print ("Molecules Generated: {}".format(len(modified_molecules)))

        return modified_molecules
    
    def validate(self, sanitizeFlags=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL):
        
        """

        Remove duplicated molecules and sanitize them.

        Arguments:
            self (Object): Cocktail object of the list of molecules
            
            sanitizeFlags (rdkit.Chem.rdmolops.SanitizeFlags): RDKit sanitization processes to go through

        Return:
            modified_molecules (List): List of the RDKit molecule objects that have had their structures replaced.
        """

        return [mol for mol in self.modified_molecules \
                if Chem.SanitizeMol(x, sanitizeOps=sanitizeFlags catchErrors=True) == Chem.rdmolops.SanitizeFlags.SANITIZE_NONE]
    
    def enumerate(self, enumeration_complexity='1D', dimensionality=None):

        """

        Enumerate the drug library based on dimension.

        Arguments:
            molecules (List): a list of molecules that the user would like enumerated.
            enumeration_complexity (String): Declares how many times we will want to discover another molecule
                                             configuration
            dimensionality (String): Enumerate based on dimensionality (1D, 2D, 3D)

        Returns:
            enumerated_molecules (List): Dependent on the dimensionality of the user it can be -> a list of smiles, or
                                         a list of RDKit Molecule Objects.

        """

        # Enumeration comes from the user iwatobipen
        # https://iwatobipen.wordpress.com/2018/11/15/generate-possible-list-of-smlies-with-rdkit-rdkit/

        print ("Enumerating Compunds....")

        if enumeration_complexity.lower() == 'low':
            complexity = 10
        elif enumeration_complexity.lower() == 'medium':
            complexity = 100
        elif enumeration_complexity.lower() == 'high':
            complexity = 1000
        else:
            complexity = 10

        enumerated_molecules = []
        for molecule in self.modified_molecules:
            if dimensionality == '1D' and smiles_enumerated not in enumerated_molecules:
                for i in range(complexity):
                    enumerated_molecules.append(smiles_enumerated)
            elif dimensionality == '2D' and smiles_enumerated not in enumerated_molecules:
                    enumerated_molecules.append(Chem.MolFromSmiles(smiles_enumerated))
            elif dimensionality == '3D' and smiles_enumerated not in enumerated_molecules:
                    Chem.rdDistGeom.EmbedMolecule(Chem.MolFromSmiles(smiles_enumerated), Chem.rdDistGeom.ETKDGv2())
                    enumerated_molecules.append(Chem.MolFromSmiles(smiles_enumerated))
        return enumerated_molecules



