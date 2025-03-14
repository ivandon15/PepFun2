"""
PepFun 2.0: improved protocols for the analysis of natural and modified peptides

From publication "PepFun 2.0: improved protocols for the analysis of natural and modified peptides"
Author: Rodrigo Ochoa
Year: 2023
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa"]
__license__ = "MIT"
__version__ = "2.0"
__email__ = "rodrigo.ochoa@udea.edu.co, rodrigo.ochoa@boehringer-ingelheim.com"

########################################################################################
# Modules to import
########################################################################################

# System dependent
import os
import math
import warnings
import numpy as np

# Modeller
import modeller
from modeller.automodel import *

# BioPython
from Bio.Align import PairwiseAligner
from Bio.PDB.Polypeptide import is_aa
from Bio.Data.SCOPData import protein_letters_3to1 as aa3to1
from Bio.PDB import PDBIO, PDBParser, NeighborSearch
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.vectors import Vector, rotaxis, calc_dihedral, calc_angle

import collections
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

BB_SMRT = "[N:1][C:2][C:3](=[O:4])[O:5]"
BB_ATOMS = ['N', 'CA', 'C', 'O', 'OXT']
BBS_SMRT = "[N:1][C:2][C:3](=[S:4])[O:5]"
BBS_ATOMS = ['N', 'CA', 'C', 'S', 'OXT']
########################################################################################
# Classes and functions
########################################################################################

class Filling:

    @staticmethod
    def fill_peptide_complex(code, template, model, path="", output_name=""):
        """
        Function to model new monomers of a bound peptide. It is assumed the peptide is the second smaller chain in the PDB

        Arguments:
        :param code: code of the structure used to model the peptide
        :param template: sequence of the peptide in the complex structure
        :param model: new sequence we want to model
        :param output_name: name of the final modelled structure

        :return: PDB structure of the protein-peptide complex with new amino acids
        """
        print("Fill a peptide in complex with a protein target")

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1]=="/":
                os.system("cp {}{}.pdb .".format(path,code))
            else:
                os.system("cp {}/{}.pdb .".format(path,code))

        structure = parser.get_structure('PEP', code + ".pdb")
        chains=[]
        for i,chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 1:
                if chSeq != template:
                    print("The second chain is not the peptide sequence. Please verify the input")
                    exit(1)
        if len(chains) != 2:
            print("The complex does not have exactly 2 chains (protein and peptide). Please verify the input")
            exit(1)

        aligner = PairwiseAligner()
        alignments = aligner.align(template, model)
        alignment = alignments[0]
        pepT=str(alignment).split()[0]
        pepM=str(alignment).split()[2]
        print("The template sequence is: {}, and the sequence to model is: {}".format(pepT, pepM))

        # Start Modeller environment and write the sequence file
        e = modeller.Environ()
        m = modeller.Model(e, file=code)
        aln = modeller.Alignment(e)
        aln.append_model(m, align_codes=code)
        aln.write(file=code +'.seq')

        # Obtain the protein sequence with the peptide
        infoSeq =[x.strip() for x in open(code +'.seq')]
        header =[]
        sequenceLine =''
        for info in infoSeq:
            if ">" not in info and ":" not in info:
                if info:
                    sequenceLine += info
            else:
                if info: header.append(info)

        # Split the sequence information as required by Modeller
        last_line =sequenceLine.split("/")
        sequenceTemp =last_line[0 ] +"/ " +pepT +"*"
        sequenceMod =last_line[0 ] +"/ " +pepM +"*"

        seqTempList = [sequenceTemp[i: i +75] for i in range(0, len(sequenceTemp), 75)]
        seqModList = [sequenceMod[i: i +75] for i in range(0, len(sequenceMod), 75)]

        # Create the alignment file
        alignmentFile =open("alignment.ali" ,"w")
        for h in header: alignmentFile.write( h +"\n")
        for s in seqTempList: alignmentFile.write( s +"\n")
        alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
        for s in seqModList: alignmentFile.write( s +"\n")
        alignmentFile.close()

        # Directories for input atom files
        e.io.atom_files_directory = ['.', '../atom_files']
        a = AutoModel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
        a.starting_model= 1
        a.ending_model  = 1
        a.make()

        # Export the generated model
        if output_name:
            if ".pdb" in output_name:
                os.system('mv {}_fill.B99990001.pdb {}'.format(code, output_name))
            else:
                os.system('mv {}_fill.B99990001.pdb {}.pdb'.format(code, output_name))
        else:
            os.system('mv {}_fill.B99990001.pdb filled_{}.pdb'.format(code, code))

        os.system("rm {c}.seq {c}_fill.D00000001 {c}_fill.ini {c}_fill.rsr {c}_fill.sch {c}_fill.V99990001 alignment.ali".format(c=code))
        if path:
            os.system("rm {}.pdb".format(code))


    ########################################################################################
    @staticmethod
    def fill_peptide_alone(code, template, model, path="", output_name=""):
        """
        Function to model new monomers of a single peptide.

        Arguments:
        :param code: code of the structure used to model the peptide
        :param template: sequence of the peptide
        :param model: new sequence we want to model
        :param output_name: name of the final modelled structure

        :return: PDB structure of the peptide with new amino acids
        """
        print("Fill a single peptide")

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1]=="/":
                os.system("cp {}{}.pdb .".format(path,code))
            else:
                os.system("cp {}/{}.pdb .".format(path,code))

        structure = parser.get_structure('PEP', code + ".pdb")
        chains = []
        for i, chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 0:
                if chSeq != template:
                    print("The template peptide is not alone in the structure. Please verify the input")
                    exit(1)
        if len(chains) > 1:
            print("The PDB must have only 1 chain with the peptide")
            exit(1)

        aligner = PairwiseAligner()
        alignments = aligner.align(template, model)
        alignment = alignments[0]
        pepT = str(alignment).split()[0]
        pepM = str(alignment).split()[2]
        print("The template sequence is: {}, and the sequence to model is: {}".format(pepT, pepM))

        # Start the Modeller environment
        e = modeller.Environ()
        m = modeller.Model(e, file=code)
        aln = modeller.Alignment(e)
        aln.append_model(m, align_codes=code)
        aln.write(file=code+'.seq')

        # Edit the information of the sequence to store the sequences in the Modeller format
        infoSeq=[x.strip() for x in open(code+'.seq')]
        header=[]
        sequenceLine=''
        for info in infoSeq:
            if ">" not in info and ":" not in info:
                if info:
                    sequenceLine+=info
            else:
                if info: header.append(info)

        # Store the sequences in variables according to Modeller format
        sequenceTemp=pepT+"*"
        sequenceMod=pepM+"*"

        seqTempList = [sequenceTemp]
        seqModList = [sequenceMod]

        # Create the alignment file
        alignmentFile=open("alignment.ali","w")
        for h in header: alignmentFile.write(h+"\n")
        for s in seqTempList: alignmentFile.write(s+"\n")
        alignmentFile.write("\n>P1;{}_fill\nsequence:::::::::\n".format(code))
        for s in seqModList: alignmentFile.write(s+"\n")
        alignmentFile.close()

        # Directories for input atom files
        e.io.atom_files_directory = ['.', '../atom_files']
        a = AutoModel(e, alnfile='alignment.ali', knowns='{}'.format(code), sequence='{}_fill'.format(code))
        a.starting_model= 1
        a.ending_model  = 1
        a.make()

        if output_name:
            if ".pdb" in output_name:
                os.system('mv {}_fill.B99990001.pdb {}'.format(code, output_name))
            else:
                os.system('mv {}_fill.B99990001.pdb {}.pdb'.format(code, output_name))
        else:
            os.system('mv {}_fill.B99990001.pdb filled_{}.pdb'.format(code, code))

        os.system("rm {c}.seq {c}_fill.D00000001 {c}_fill.ini {c}_fill.rsr {c}_fill.sch {c}_fill.V99990001 alignment.ali".format(c=code))
        if path:
            os.system("rm {}.pdb".format(code))

    # End of Filling class
    ############################################################

########################################################################################
class Capping:

    ########################################################################################
    @staticmethod
    def capping_peptide(pdb_file, peptide, path="", output_name = "", mode='both'):
        '''
        Function to cap an existing PDB file

        :param pdb_file: PDB of the peptide
        :param peptide: sequence of the original peptide
        :param output_name: name of the output capped structure
        :param mode: Mode to select which positions to cap. Choose between, N-term, C-term, or both

        :return: a peptide capped structure
        '''

        print("Capping a peptide with ACE and/or NME groups")

        template=peptide
        if mode in ['N-term','C-term','both']:
            if mode=='N-term':
                model='G' + peptide
            elif mode=='C-term':
                model = peptide + 'G'
            elif mode=='both':
                model = 'G' + peptide + 'G'
        else:
            print("The mode option should be N-term, C-term or both. Please check the option")
            exit(1)
        code=pdb_file.split('.')[0]

        def _get_pdb_sequence(structure):
            """
            Retrieves the AA sequence from a PDB structure.
            """
            _aainfo = lambda r: aa3to1.get(r.resname, "X")
            seq = [_aainfo(r) for r in structure.get_residues() if is_aa(r)]
            final_seq = ''.join(seq)

            return final_seq

        # Verifying the chains and sequence
        parser = PDBParser()
        if path:
            if path[-1] == "/":
                os.system("cp {}{} .".format(path, pdb_file))
            else:
                os.system("cp {}/{} .".format(path, pdb_file))
        structure = parser.get_structure('PEP', pdb_file)
        chains = []
        for i, chain in enumerate(structure[0]):
            chains.append(chain)
            chSeq = _get_pdb_sequence(chain)

            if i == 0:
                if chSeq != peptide:
                    print("The peptide sequence does not coincide with the PDB file. Please verify input")
                    exit(1)
        try:
            Filling.fill_peptide_alone(code, template, model, output_name='modeller_peptide')
        except Exception as e:
            print("Failed to fill the peptide with cap groups, Error: {}".format(e))
            exit(1)

        parser = PDBParser()

        reference = parser.get_structure('REF', "modeller_peptide.pdb")
        model = reference[0]
        chain = model['A']

        for i, res in enumerate(chain):
            if i == 0:
                if mode == 'both' or mode == 'N-term':
                    # Delete the other atoms leaving only the atoms of the backbone
                    ids = []
                    for a in res:
                        atomId = a.id
                        if atomId not in ("CA", "O", "C"): ids.append(atomId)
                    for i in ids: res.detach_child(i)

                    res.resname = 'ACE'

                    for a in res:
                        atomId = a.id
                        if atomId == 'CA':
                            a.fullname = 'CH3'
                            a.name = 'CH3'
                            a.element = 'C'

            if i == len(chain)-1:
                if mode == 'both' or mode == 'C-term':
                    # Delete the other atoms leaving only the atoms of the backbone
                    ids = []
                    for a in res:
                        atomId = a.id
                        if atomId not in ("CA", "N"): ids.append(atomId)
                    for i in ids: res.detach_child(i)

                    res.resname = 'NME'

                    for a in res:
                        atomId = a.id
                        if atomId == 'CA':
                            a.fullname = ' C'
                            a.name = ' C'
                            a.element = 'C'

        # Saving the new structure
        io = PDBIO()
        io.set_structure(reference)
        if output_name:
            if ".pdb" in output_name:
                io.save("{}".format(output_name))
            else:
                io.save("{}.pdb".format(output_name))
        else:
            io.save("capped_{}".format(pdb_file))

        os.system("rm modeller_peptide.pdb")
        if path:
            os.system("rm {}".format(pdb_file))

    # End of Capping class
    ############################################################

########################################################################################

class Mutation:

    def __init__(self, pdb_file, chainID, pepPosition, nnaa_name, pdbNNAA, smilesNNAA=None, mutated_file=None, ignore_pairs=None):
        '''
        Class to mutate a NNAA in a position of a natural AA

        :param pdb_file: PDB file of the peptide that will be mutated
        :param chainID: chain ID of the peptide
        :param pepPosition: position on the peptide that will be mutated
        :param pdbNNAA: PDB file of the NNAA used in the mutation
        :param smilesNNAA: if NNAA's smiles is provided, then generate and save in pdbNNAA
        :param chainNNAA: chain of the NNAA PDB file
        :param nnaa_name: abbreviated name of the NNAA
        :param ignore_pairs: a list of tuple pair(in order), should contains the position of cyclization, [(1,13), (3,13)]
        '''

        self.pdb_file = pdb_file
        self.chainID = chainID
        self.pepPosition = pepPosition
        self.pdbNNAA = pdbNNAA
        self.smilesNNAA = smilesNNAA
        self.chainNNAA = "A"
        self.nnaa_name = nnaa_name
        self.mutated_file = "mutated_{}.pdb".format(self.nnaa_name) if not mutated_file else mutated_file
        self.ignore_pairs = ignore_pairs
        self.BBS_FLG = False

    ########################################################################################

    def get_bb_dihedral(self):
        '''
        Function to get the backbone dihedral from the structure that will be mutated

        :return: The backbone dihedrals of the input PDB file
        '''

        # Read the reference PDB file
        parser = PDBParser()
        reference = parser.get_structure('REF',self.pdb_file)
        model = reference[0]
        chain = model[self.chainID]

        # Inititalize variables
        rad = 180.0 / math.pi
        flagCap = False
        
        # deal with HETAM
        resid_map = {}
        for i,res in enumerate(chain):
            resid_map[i] = res.get_id()

        for i,res in enumerate(chain):
            # First residue
            if i == 0:
                try:
                    phi = 0
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()

                    N_next = chain[resid_map[i + 1]]["N"]
                    n_next = N_next.get_vector()

                    psi = calc_dihedral(n, ca, c, n_next)

                    # Store the data for the position of interest
                    if self.pepPosition == i + 1:
                        print(c,n,ca)
                        print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(),self.pepPosition, phi * rad, psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad

                    # Update the previous coordinates
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                except:
                    flagCap=True

            # Last residue
            elif i == len(chain)-1:
                try:
                    n1 = N_prev.get_vector()
                    ca1 = CA_prev.get_vector()
                    c1 = C_prev.get_vector()

                    # Get current AA coordinates
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()

                    # Calculate dihedral
                    psi = 0
                    phi = calc_dihedral(c1, n, ca, c)

                    # Store the data for the position of interest
                    if self.pepPosition == i + 1:
                        print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(),self.pepPosition, phi * rad,psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad
                except:
                    flagCap=True
            else:
                if flagCap is False:
                    # Get previous AA coordinates
                    n1 = N_prev.get_vector()
                    ca1 = CA_prev.get_vector()
                    c1 = C_prev.get_vector()

                    # Get current AA coordinates
                    C_curr = res["C"]
                    N_curr = res["N"]
                    CA_curr = res["CA"]
                    c = C_curr.get_vector()
                    n = N_curr.get_vector()
                    ca = CA_curr.get_vector()
                    # print("resid_map[i + 2], raw:", resid_map[i + 1], i + 2)
                    N_next = chain[resid_map[i + 1]]["N"]
                    n_next = N_next.get_vector()

                    # Calculate dihedral
                    psi = calc_dihedral(n, ca, c, n_next)
                    phi = calc_dihedral(c1, n, ca, c)

                    # Store the data for the position of interest
                    if self.pepPosition == i+1:
                        # print('The dihedrals for residue {}{} are: Phi ({}) - Psi ({})'.format(res.get_resname(), self.pepPosition, phi * rad, psi * rad))
                        self.ref_psi = psi * rad
                        self.ref_phi = phi * rad

                    # Update the previous coordinates
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                else:
                    N_prev = res["N"]
                    CA_prev = res["CA"]
                    C_prev = res["C"]
                    flagCap=False

        try:
            self.near_phi = int(self.ref_phi / 10) * 10
            self.near_psi = int(self.ref_psi / 10) * 10
        except:
            print("Failed to map the backbone dihedrals. Check the input variables")

    ########################################################################################
    # HACK
    def check_clashes(self, pdb_file, chain_id, distance_threshold=2.1, ignore_pairs=None):
        '''
        Check for clashes between atoms in a specified chain and other atoms in the structure.

        :param pdb_file: Path to the PDB file of the structure.
        :param chain_id: ID of the chain to check for clashes.
        :param distance_threshold: Distance threshold to define a clash (default: 2.1 Å).
        :return: List of clash pairs, where each pair contains the residue number and atom name of the clashing atoms.
        '''
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        model = structure[0]

        atom_list = list(model.get_atoms())
        key_atoms = [atom for atom in atom_list if atom.get_id()[0] != 'H']

        ns = NeighborSearch(key_atoms)

        def is_backbone(atom_name):
            return atom_name in ["N", "CA", "C", "O"]
        
        clash_pairs = set()
        for chain in model:
            if chain.id != chain_id:
                continue  
            for residue in chain:
                for atom in residue:
                    if is_backbone(atom.get_name()):
                        continue
                    if atom.get_id()[0] == 'H':
                        continue
                    neighbors = ns.search(atom.coord, distance_threshold)

                    for neighbor in neighbors:
                        if neighbor.get_parent() == residue:
                            continue
                        pair = tuple(sorted((residue.id[1], neighbor.get_parent().id[1])))
                        if not ignore_pairs or pair not in ignore_pairs:
                            clash_pairs.add(pair)
        return clash_pairs


    def get_connect(self, pdbfile):
        '''
        Extract CONECT lines from the NNAA PDB file.

        :return: List of CONECT lines from the NNAA PDB file.
        '''
        conect_lines = []
        with open(pdbfile, "r") as f:
            for line in f:
                if line.startswith("CONECT"):
                    conect_lines.append(line)
        return conect_lines


    def add_connect(self, mutated_pdb_file):
        '''
        Add CONECT lines from the NNAA PDB file to the mutated PDB file.

        :param mutated_pdb_file: Path to the mutated PDB file where CONECT lines will be added.
        :return: None. The CONECT lines are appended to the mutated PDB file.
        '''
        parser = PDBParser()
        reference = parser.get_structure('REF', mutated_pdb_file)
        model = reference[0]
        chain = model[self.chainID]
        conect_lines = self.get_connect(self.pdbNNAA)
        atom_idx_map = {}        
        for i, res in enumerate(chain):
            if self.pepPosition == i + 1:
                for atom in res.get_atoms():
                    atom_idx_map[self.name2wtidx[atom.get_id()]] = atom.serial_number
        
        updated_conect_lines = []
        for line in conect_lines:
            parts = line.split()
            new_parts = [parts[0]]
            
            for atom_id in parts[1:]:
                new_id = atom_idx_map.get(int(atom_id), ' ')
                new_parts.append(str(new_id))
                
            updated_conect_lines.append(" ".join(new_parts))

        with open(mutated_pdb_file, 'a') as f:
            for conect in updated_conect_lines:
                f.write(conect + '\n')


    def generate_conformer(self, smiles):
        """
        Generate a 3D conformer for a molecule from its SMILES string.

        Arguments:
        :param smiles: SMILES string of the molecule.
        :return: RDKit Mol object with embedded 3D coordinates.
        """
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)

        ps = AllChem.ETKDGv2()
        id = AllChem.EmbedMolecule(mol, ps)
        if id == -1:
            # use random coords
            ps.useRandomCoords = True
            AllChem.EmbedMolecule(mol, ps)
            AllChem.MMFFOptimizeMolecule(mol, confId=0)
        return mol
    

    def _atomidx2atomname(self, mol, keep_existing_names=False):
        """
        Map atom indices to atom names for a molecule.

        Arguments:
        :param mol: RDKit Mol object representing the molecule.
        :param keep_existing_names: If True, preserve existing atom names; otherwise, generate new names.
        :return: A dictionary mapping atom indices to atom names.
        """
        mol = Chem.Mol(mol)
        backbone =  Chem.MolFromSmarts(BB_SMRT) 
        backbone_idxes = mol.GetSubstructMatches(backbone)
        bb_atoms = BB_ATOMS
        if not backbone_idxes:
            backbone =  Chem.MolFromSmarts(BBS_SMRT) 
            backbone_idxes = mol.GetSubstructMatches(backbone)
            bb_atoms = BBS_ATOMS
            self.BBS_FLG = True

        backbone_idxes = backbone_idxes[0]
        backbone_map = dict(zip(backbone_idxes, bb_atoms))
        origin_new_name_map = {}
        element_counts = collections.Counter()
        for atom in mol.GetAtoms():
            if not atom.HasProp('atom_name') or not keep_existing_names:
                atom_idx = atom.GetIdx()
                if atom_idx in backbone_map.keys():
                    new_name = backbone_map[atom_idx]
                else:
                    element = atom.GetSymbol()
                    element_counts[element] += 1
                    new_name = f'{element}{element_counts[element]}'
                origin_new_name_map[atom_idx+1] = new_name
        return origin_new_name_map
    

    def rename_pdb(self, input_pdb, output_pdb, mol, aa_name):
        """
        Rename atoms in a PDB file based on a provided index-to-name mapping.

        Arguments:
        :param input_pdb: Path to the input PDB file.
        :param output_pdb: Path to the output PDB file.
        :param mol: Mol object need to rename
        :param aa_name: Amino acid name to be added to the PDB file.
        :return: None
        """
        idx_map = self._atomidx2atomname(mol)

        with open(input_pdb, "r") as f_in, open(output_pdb, "w") as f_out:
            for line in f_in:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    atom_idx = int(line[6:11].strip())
                    if atom_idx in idx_map:
                        new_atom_name = idx_map[atom_idx].ljust(4)[:4]
                        line = line[:13] + new_atom_name + aa_name.ljust(4)[:4] + "A" + line[22:]
                f_out.write(line)


    def smiles2pdb(self, smiles, pdbNNAA, nnaa_name):
        """
        Convert a SMILES string to a PDB file with renamed atoms.

        Arguments:
        :param smiles: SMILES string of the molecule.
        :param aa_name: Amino acid name to be added to the PDB file.
        :param pdb_file: Path to the output PDB file.
        :return: None
        """
        mol = self.generate_conformer(smiles)
        tmp_pdb_file = pdbNNAA.replace(".","tmp.")
        Chem.MolToPDBFile(mol, tmp_pdb_file)
        self.rename_pdb(tmp_pdb_file, pdbNNAA, mol, nnaa_name)
        os.remove(tmp_pdb_file)


    def find_connected_atoms(self, res, start_atom_name, stop_atom_names, ns):
        """
        Find all atoms connected to start_atom_name using BFS, stopping at stop_atom_names.
        Record the order in which atoms are traversed and the edges (pairs of connected atoms).

        :param res: The residue object (Bio.PDB.Residue).
        :param start_atom_name: The name of the starting atom.
        :param stop_atom_names: A set of atom names at which to stop traversal.
        :param ns: NeighborSearch object for finding connected atoms.
        :return: A tuple containing:
                - A list of atom names in the order they were traversed.
                - A list of tuples representing edges (pairs of connected atoms).
        """
        # Initialize BFS queue with the starting atom
        queue = deque()
        queue.append(start_atom_name)

        # Set to track visited atoms
        visited = set()

        # List to store the order of traversal
        traversal_order = []

        # List to store the edges (pairs of connected atoms)
        edges = []

        while queue:
            current_atom_name = queue.popleft()

            # Skip if already visited or if it's a stop atom
            if current_atom_name in visited or current_atom_name in stop_atom_names:
                continue

            # Mark the current atom as visited
            visited.add(current_atom_name)
            traversal_order.append(current_atom_name)

            # Get the current atom object
            if current_atom_name not in res:
                continue
            current_atom = res[current_atom_name]

            # Find connected atoms
            connected_atoms_list = ns.search(current_atom.coord, 2.0)
            connected_names = {atom.get_id() for atom in connected_atoms_list if atom.get_id() != current_atom_name}

            # Add connected atoms to the queue and record edges
            for atom_name in connected_names:
                if atom_name not in visited and atom_name not in stop_atom_names:
                    queue.append(atom_name)
                    edges.append((current_atom_name, atom_name))  # Record the edge

        return traversal_order, edges


    def has_multiple_paths_N_to_CA(self, residue):
        """
        Check if there are multiple paths between the N atom and the CA atom in a residue.

        :param residue: A Bio.PDB.Residue object representing the residue.
        :return: True if there are multiple paths, False otherwise.
        """
        # Get the N and CA atoms
        n_atom = residue["N"]
        ca_atom = residue["CA"]

        # Initialize a NeighborSearch object to find bonded atoms
        ns = NeighborSearch(list(residue.get_atoms()))

        # Set to store visited atoms to avoid cycles
        visited = set()

        # Counter to track the number of unique paths
        path_count = 0

        def dfs(current_atom, target_atom, path):
            nonlocal path_count
            if current_atom == target_atom:
                path_count += 1
                if path_count > 1:
                    return True  # Early exit if multiple paths are found
                return False

            visited.add(current_atom)
            for neighbor in ns.search(current_atom.coord, 2.0):  # 2.0 Å is a typical bond cutoff
                if neighbor not in visited and neighbor != current_atom:
                    if dfs(neighbor, target_atom, path + [neighbor]):
                        return True
            visited.remove(current_atom)
            return False

        # Start DFS from N atom
        dfs(n_atom, ca_atom, [n_atom])

        # Return True if more than one path exists
        return path_count > 1

        
    def get_NNAA_basic_info(self):
        # Input PDB information
        rad = 180.0 / math.pi
        parser = PDBParser()
        if self.smilesNNAA:
            self.smiles2pdb(self.smilesNNAA, self.pdbNNAA, self.nnaa_name)

        reference = parser.get_structure('REF', self.pdbNNAA)
        model = reference[0]
        chain = model[self.chainNNAA]
        
        # Get all the atoms in the NNAA and store those that are not hydrogens
        atom_list = [x for x in model.get_atoms()]
        self.bonds = []
        self.assoc_bonds = {}
        self.name2wtidx = {atom.get_id(): i+1 for i, atom in enumerate(model.get_atoms())}
        self.sidechain_atom_order = [] # sidechain
        self.Nmod_atom_order = [] # N modification chain
        key_atoms = [atom for atom in atom_list if atom.get_id()[0] != 'H']
        ns = NeighborSearch(key_atoms)
        for res in chain:
            is_N2CA_cyclic = self.has_multiple_paths_N_to_CA(res)
            if not is_N2CA_cyclic:
                # eluminate cyclic situations
                self.Nmod_atom_order, bonds = self.find_connected_atoms(res, "N", {"CA"}, ns)
                if len(self.Nmod_atom_order) == 1:
                    # no N modification
                    self.Nmod_atom_order = []
                else:
                    self.bonds.extend(bonds)
            self.sidechain_atom_order, bonds = self.find_connected_atoms(res, "CA", {"N", "C"}, ns)
            if len(self.sidechain_atom_order) == 1:
                # no sidechain
                self.sidechain_atom_order = []
            else:
                self.bonds.extend(bonds)
            break
        for ele1, ele2 in self.bonds:
            self.assoc_bonds[ele2] = ele1

        self.dict_NNAA = {}
        for bond in self.bonds:
            self.dict_NNAA[bond]={}
            ele1 = bond[0]
            ele2 = bond[1]
            if sorted(list(bond)) == sorted(["N", "CA"]):
                # no need to consider
                continue
            if ele1 in self.sidechain_atom_order:
                ind1 = self.sidechain_atom_order.index(ele1)
                # Exception
                if ele1 != "CA" and self.assoc_bonds[ele1] == "CA":
                    ind1 = 1
                for i,res in enumerate(chain):

                    if ind1 == 0:
                        atom_1 = res["N"]
                        atom_2 = res["C"]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]
                    elif ind1 == 1:
                        atom_1 = res["N"]
                        atom_2 = res[self.assoc_bonds[ele1]]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]
                    else:
                        atom_1 = res[self.assoc_bonds[self.assoc_bonds[ele1]]]
                        atom_2 = res[self.assoc_bonds[ele1]]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]

                    atom1_vec = atom_1.get_vector()
                    atom2_vec = atom_2.get_vector()
                    atom3_vec = atom_3.get_vector()
                    atom4_vec = atom_4.get_vector()

                    length = atom_4 - atom_3
                    angle = calc_angle(atom2_vec, atom3_vec, atom4_vec) * rad
                    diangle = calc_dihedral(atom1_vec, atom2_vec, atom3_vec, atom4_vec) * rad

                    self.dict_NNAA[bond]={'length': length, 'angle': angle, 'diangle': diangle}
            else:
                ind1 = self.Nmod_atom_order.index(ele1)
                # Exception
                if ele1 != "N" and self.assoc_bonds[ele1] == "N":
                    ind1 = 1
                for i,res in enumerate(chain):

                    if ind1 == 0:
                        atom_1 = res["CA"]
                        atom_2 = res["C"]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]
                    elif ind1 == 1:
                        atom_1 = res["CA"]
                        atom_2 = res[self.assoc_bonds[ele1]]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]
                    else:
                        atom_1 = res[self.assoc_bonds[self.assoc_bonds[ele1]]]
                        atom_2 = res[self.assoc_bonds[ele1]]
                        atom_3 = res[ele1]
                        atom_4 = res[ele2]

                    atom1_vec = atom_1.get_vector()
                    atom2_vec = atom_2.get_vector()
                    atom3_vec = atom_3.get_vector()
                    atom4_vec = atom_4.get_vector()

                    length = atom_4 - atom_3
                    angle = calc_angle(atom2_vec, atom3_vec, atom4_vec) * rad
                    diangle = calc_dihedral(atom1_vec, atom2_vec, atom3_vec, atom4_vec) * rad

                    self.dict_NNAA[bond]={'length': length, 'angle': angle, 'diangle': diangle}
    ########################################################################################
    def assign_mutation(self):
        '''
        Function to assign the mutation based on the stored parameters

        :return: A mutated PDB file that can be minimized using external tools
        '''

        parser = PDBParser()
        reference = parser.get_structure('REF', self.pdb_file)
        model = reference[0]
        chain = model[self.chainID]

        for i, res in enumerate(chain):
            # 遍历sequence
            if self.pepPosition == i + 1:
                # 找到突变位置
                C = res["C"]
                N = res["N"]
                CA = res["CA"]
                # used for following connection information
                res.id = ('H_', res.id[1], res.id[2])
                # Delete the other atoms leaving only the atoms of the backbone
                ids = []
                for a in res:
                    atomId = a.id
                    if atomId not in ("N", "CA", "O", "C"): ids.append(atomId)
                for i in ids: res.detach_child(i)

                # deal with thio-amid
                if self.BBS_FLG:
                    O_atom = res["O"]
                    S_atom = Atom(
                        name="S",
                        coord=O_atom.coord,
                        bfactor=O_atom.bfactor,
                        occupancy=O_atom.occupancy,
                        altloc=O_atom.altloc,
                        fullname=" S ",
                        serial_number=O_atom.serial_number,
                        element="S"  
                    )
                    res.detach_child("O")
                    res.add(S_atom)

                res.resname = self.nnaa_name
                atom_to_select = {'CA': CA, 'N': N}
                for bond in self.bonds:
                    ele1 = bond[0]
                    ele2 = bond[1]
                    if ele1 in self.sidechain_atom_order:
                        ind1 = self.sidechain_atom_order.index(ele1)
                        if ele1 != "CA" and self.assoc_bonds[ele1] == "CA":
                            ind1 = 1
                        if ind1 == 0:
                            atom_coord = calculateCoordinates(
                                    N, C, atom_to_select[ele1], self.dict_NNAA[bond]['length'], self.dict_NNAA[bond]['angle'],
                                    self.dict_NNAA[bond]['diangle']
                            )
                            atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
                        elif ind1 == 1:
                            atom_coord = calculateCoordinates(
                                    N, atom_to_select[self.assoc_bonds[ele1]], atom_to_select[ele1],
                                    self.dict_NNAA[bond]['length'],
                                    self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                            )
                            atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
                        else:
                            atom_coord = calculateCoordinates(
                                    atom_to_select[self.assoc_bonds[self.assoc_bonds[ele1]]], atom_to_select[self.assoc_bonds[ele1]],
                                    atom_to_select[ele1], self.dict_NNAA[bond]['length'],
                                    self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                            )
                            if ele2[0] not in ('C', 'N', 'O', 'S'):
                                if ele2 == 'BRZ':
                                    ele2mod = 'BrZ'
                                else:
                                    ele2mod = ele2
                                atom_object = Atom(ele2mod, atom_coord, 0.0, 1.0, " ", " {}".format(ele2mod), 0,
                                                ele2[0:len(ele2) - 1])
                            else:
                                atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
                    else:
                        ind1 = self.Nmod_atom_order.index(ele1)
                        if ele1 != "N" and self.assoc_bonds[ele1] == "N":
                            ind1 = 1
                        if ind1 == 0:
                            atom_coord = calculateCoordinates(
                                    CA, C, atom_to_select[ele1], self.dict_NNAA[bond]['length'], self.dict_NNAA[bond]['angle'],
                                    self.dict_NNAA[bond]['diangle']
                            )
                            atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
                        elif ind1 == 1:
                            atom_coord = calculateCoordinates(
                                    CA, atom_to_select[self.assoc_bonds[ele1]], atom_to_select[ele1],
                                    self.dict_NNAA[bond]['length'],
                                    self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                            )
                            atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
                        else:
                            atom_coord = calculateCoordinates(
                                    atom_to_select[self.assoc_bonds[self.assoc_bonds[ele1]]], atom_to_select[self.assoc_bonds[ele1]],
                                    atom_to_select[ele1], self.dict_NNAA[bond]['length'],
                                    self.dict_NNAA[bond]['angle'], self.dict_NNAA[bond]['diangle']
                            )
                            if ele2[0] not in ('C', 'N', 'O', 'S'):
                                if ele2 == 'BRZ':
                                    ele2mod = 'BrZ'
                                else:
                                    ele2mod = ele2
                                atom_object = Atom(ele2mod, atom_coord, 0.0, 1.0, " ", " {}".format(ele2mod), 0,
                                                ele2[0:len(ele2) - 1])
                            else:
                                atom_object = Atom(ele2, atom_coord, 0.0, 1.0, " ", " {}".format(ele2), 0, ele2[0])
                            if ele2 not in atom_to_select:
                                atom_to_select[ele2] = atom_object
                                res.add(atom_object)
        # Saving the new structure
        io = PDBIO()
        io.set_structure(reference)
        io.save(self.mutated_file)
        self.add_connect(self.mutated_file)
        if Chem.MolFromPDBFile(self.mutated_file) is None:
            raise Exception
        clash_pairs = self.check_clashes(self.mutated_file, self.chainID, ignore_pairs=self.ignore_pairs)
        if clash_pairs:
            # print(clash_pairs)
            raise Exception

    # End of Mutation class
    ############################################################

########################################################################################
# Function derived from the PeptideBuilder module
########################################################################################

def calculateCoordinates(
    refA: Residue, refB: Residue, refC: Residue, L: float, ang: float, di: float
) -> np.ndarray:
    '''
    Function from PeptideBuilder to generate the coordinates of the NNAA based on the stored parameters
    '''
    AV = refA.get_vector()
    BV = refB.get_vector()
    CV = refC.get_vector()

    CA = AV - CV
    CB = BV - CV

    ##CA vector
    AX = CA[0]
    AY = CA[1]
    AZ = CA[2]

    ##CB vector
    BX = CB[0]
    BY = CB[1]
    BZ = CB[2]

    ##Plane Parameters
    A = (AY * BZ) - (AZ * BY)
    B = (AZ * BX) - (AX * BZ)
    G = (AX * BY) - (AY * BX)

    ##Dot Product Constant
    F = math.sqrt(BX * BX + BY * BY + BZ * BZ) * L * math.cos(ang * (math.pi / 180.0))

    ##Constants
    const = math.sqrt(
        math.pow((B * BZ - BY * G), 2)
        * (
            -(F * F) * (A * A + B * B + G * G)
            + (
                B * B * (BX * BX + BZ * BZ)
                + A * A * (BY * BY + BZ * BZ)
                - (2 * A * BX * BZ * G)
                + (BX * BX + BY * BY) * G * G
                - (2 * B * BY) * (A * BX + BZ * G)
            )
            * L
            * L
        )
    )
    denom = (
        (B * B) * (BX * BX + BZ * BZ)
        + (A * A) * (BY * BY + BZ * BZ)
        - (2 * A * BX * BZ * G)
        + (BX * BX + BY * BY) * (G * G)
        - (2 * B * BY) * (A * BX + BZ * G)
    )

    X = (
        (B * B * BX * F) - (A * B * BY * F) + (F * G) * (-A * BZ + BX * G) + const
    ) / denom

    if (B == 0 or BZ == 0) and (BY == 0 or G == 0):
        const1 = math.sqrt(
            G * G * (-A * A * X * X + (B * B + G * G) * (L - X) * (L + X))
        )
        Y = ((-A * B * X) + const1) / (B * B + G * G)
        Z = -(A * G * G * X + B * const1) / (G * (B * B + G * G))
    else:
        Y = (
            (A * A * BY * F) * (B * BZ - BY * G)
            + G * (-F * math.pow(B * BZ - BY * G, 2) + BX * const)
            - A * (B * B * BX * BZ * F - B * BX * BY * F * G + BZ * const)
        ) / ((B * BZ - BY * G) * denom)
        Z = (
            (A * A * BZ * F) * (B * BZ - BY * G)
            + (B * F) * math.pow(B * BZ - BY * G, 2)
            + (A * BX * F * G) * (-B * BZ + BY * G)
            - B * BX * const
            + A * BY * const
        ) / ((B * BZ - BY * G) * denom)

    # Get the new Vector from the origin
    D = Vector(X, Y, Z) + CV
    with warnings.catch_warnings():
        # ignore inconsequential warning
        warnings.simplefilter("ignore")
        temp = calc_dihedral(AV, BV, CV, D) * (180.0 / math.pi)

    di = di - temp
    rot = rotaxis(math.pi * (di / 180.0), CV - BV)
    D = (D - BV).left_multiply(rot) + BV

    return D.get_array()

############################################################
## End of modifications.py
############################################################