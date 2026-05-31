# v1
# - major update 20230502: can find the ace for methylated and ethylated ncaa

# v2
# 

# v3
# - major update 20230521: can do N-methylation, ethylation, and cyclization
# - major update 20230523: can do write virtual atom and bond connection to mol2 file

# to do list:
# - do N-cyclized AA! by alternative substracture matching. or refer to prolin ot Beta proline NV  VIRT VIRT
# - take any order of amino acids with any name, find the mainchain, find the ACE NME, rename them
# - add substructure matching patterns to do proline like and N-alkylations
#

import sys
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse



#mol2_file="2S3R-HAMP-Ggaba.mol2"

# Parse command line arguments
parser = argparse.ArgumentParser(description='Write M session to mol2 file')
parser.add_argument('mol2file', type=str, help='The mol2 file to process')
parser.add_argument('--overwrite', action='store_true', help='Overwrite existing M session')
args = parser.parse_args()

# Read the mol2 file
with open(args.mol2file, 'r') as f:
    lines = f.readlines()



print ("starting write_M_session_v2.py ...")


filename = args.mol2file
print ("filename:", filename)

def mol_from_mol2_file(mol2_file):
    with open(mol2_file, 'r') as f:
        mol2_contents = f.read()
    
    mol = AllChem.MolFromMol2Block(mol2_contents, sanitize=False)
    
    mol2_lines = mol2_contents.split('\n')
    
    # Find the atom section
    atom_section_start = None
    atom_section_end = None
    for i, line in enumerate(mol2_lines):
        if "@<TRIPOS>ATOM" in line:
            atom_section_start = i + 1
        elif "@<TRIPOS>BOND" in line:
            atom_section_end = i
            break
    
    # Process the atom section
    if atom_section_start is not None and atom_section_end is not None:
        for line in mol2_lines[atom_section_start:atom_section_end]:
            fields = line.split()
            
            if len(fields) > 0 and fields[0].isdigit():
                # to use GetAtomWithIdx, the original index has to -1 since RDKit index starts with 0
                RDkit_atom_index = int(fields[0]) - 1
                #print("atom_index:", atom_index)
                atom_name = fields[1]
                atom_type = fields[5]
                atom_partial_charge = fields[8]
                
                atom = mol.GetAtomWithIdx(RDkit_atom_index)
                atom.SetProp("mol2_atom_index", fields[0])
                atom.SetProp("mol2_atom_name", atom_name)
                atom.SetProp("mol2_atom_type", atom_type)
                atom.SetProp("mol2_atom_partial_charge", atom_partial_charge)

    # Find the bond section
    bond_section_start = None
    bond_section_end = None
    for i, line in enumerate(mol2_lines):
        if "@<TRIPOS>BOND" in line:
            bond_section_start = i + 1
        elif "@<TRIPOS>SUBSTRUCTURE" in line:
            bond_section_end = i
            break
    # openbabel mol2 has no SUBSTRUCTURE record; bond section runs to EOF
    if bond_section_start is not None and bond_section_end is None:
        bond_section_end = len(mol2_lines)

    # Process the bond section
    if bond_section_start is not None and bond_section_end is not None:
        for line in mol2_lines[bond_section_start:bond_section_end]:
            fields = line.split()
            #print("fields:",fields)
            if len(fields) > 0 and fields[0].isdigit():
                bond_index = int(fields[0])
                RDKit_bonding_atom1_index = int(fields[1]) - 1
                RDKit_bonding_atom2_index = int(fields[2]) - 1
                bond_type = fields[3]
                
                bond = mol.GetBondBetweenAtoms(RDKit_bonding_atom1_index, RDKit_bonding_atom2_index)
                bond.SetProp("mol2_bond_index", fields[0])
                bond.SetProp("mol2_bond_type", bond_type)
                #print("bond:",bond.mol2_bond_index)
    return mol, mol2_lines, atom_section_start, atom_section_end, bond_section_start, bond_section_end


# Read the mol2 file and create the molecule object
mol2_file = filename
try:
    mol, mol2_lines, atom_section_start, atom_section_end, bond_section_start, bond_section_end = mol_from_mol2_file(mol2_file)
except Exception:
    print("Warning: The mol2 file may already contain a virtual atom. use the original mol2 file if want to overwrite")

print ("atom_section_start:",atom_section_start)
print ("atom_section_end:",atom_section_end)
print ("bond_section_start:",bond_section_start)
print ("bond_section_end:",bond_section_end)

# write the SMILE string
M__POLY_SMILES = Chem.MolToSmiles(mol)
print (M__POLY_SMILES)





##############get net charge##########
M__POLY_CHG = 0.0
for atom in mol.GetAtoms():
    partial_charge = float(atom.GetProp("mol2_atom_partial_charge"))
    M__POLY_CHG += partial_charge

# debug
# print("partial charge sum:",M__POLY_CHG)
M__POLY_CHG = round(M__POLY_CHG)




# Custom Atom class
class CustomAtom:
    def __init__(self, atom_index, atom_name, atom_type, atom_element, atom_partial_charge):
        self.atom_index = atom_index
        self.atom_name = atom_name
        self.atom_type = atom_type
        self.atom_element = atom_element
        self.atom_partial_charge = atom_partial_charge

# Assuming `mol` is your RDKit molecule object
custom_atoms = []

for atom in mol.GetAtoms():
    atom_index = int(atom.GetProp("mol2_atom_index"))
    #print("atom_index:",atom_index)
    atom_name = atom.GetProp("mol2_atom_name")
    atom_type = atom.GetProp("mol2_atom_type")
    atom_element = atom.GetSymbol()
    atom_partial_charge = atom.GetProp("mol2_atom_partial_charge")
    custom_atom = CustomAtom(atom_index, atom_name, atom_type, atom_element, atom_partial_charge)  # Pass atom_element as a parameter
    custom_atoms.append(custom_atom)
    #print ("custom_atom:",custom_atom.atom_index,)

# Accessing the custom atom object's attributes
#print(custom_atoms[0].atom_name)
#print("all atom index",custom_atoms[-1].atom_index)
#print(custom_atoms[0].atom_element)


# Custom Bond class
class CustomBond:
    def __init__(self, bond_index, bonding_atom1_index, bonding_atom2_index, bond_type):
        self.bond_index = bond_index
        self.bonding_atom1_index = bonding_atom1_index
        self.bonding_atom2_index = bonding_atom2_index
        self.bond_type = bond_type

custom_bonds = []

for bond in mol.GetBonds():
    bond_index = int(bond.GetProp("mol2_bond_index"))
    RDKit_bonding_atom1_index = int(bond.GetBeginAtomIdx())
    RDKit_bonding_atom2_index = int(bond.GetEndAtomIdx())
    bonding_atom1_index = int(bond.GetBeginAtomIdx()) + 1
    bonding_atom2_index = int(bond.GetEndAtomIdx()) + 1
    bond_type = bond.GetProp("mol2_bond_type")
    custom_bond = CustomBond(bond_index, bonding_atom1_index, bonding_atom2_index, bond_type)
    custom_bonds.append(custom_bond)
    #print ("custom_bond:",custom_bond.bonding_atom1_index,custom_bond.bonding_atom2_index)



#################get mainchain atoms#######
# Get ACE and NME matches

ace_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('O=C(N[H])C([H])([H])[H]'))
nme_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('O=CN([H])C([H])([H])[H]'))


#raise errors if no ACE, N-methylation, N-ethylation or ACE is not unique,
n_cyclized = False
if len(ace_matches) == 0:
    print("Warning: ACE atoms not found! Now try N-methylated ACE")
    ace_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('O=C(N([*])(C([H])([H])[H]))C([H])([H])[H]'))
    if not ace_matches:
        print("Warning: ACE atoms not found in N-methylated match! Now try N-ethylated ACE")
        ace_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('O=C(N([*])C([H])([H])C([H])([H])[H])C([H])([H])[H]'))
        if not ace_matches:
            # now try to find if it is the N-cyclized ACE
            print("Warning: ACE atoms not found in N-ethylated match! Now try N-cyclized ACE")
            smarts_patterns = [
                'O=C([N]1~*~*1)C([H])([H])[H]',
                'O=C([N]1~*~*~*1)C([H])([H])[H]',
                'O=C([N]1~*~*~*~*1)C([H])([H])[H]',
                'O=C([N]1~*~*~*~*~*1)C([H])([H])[H]',
                'O=C([N]1~*~*~*~*~*~*1)C([H])([H])[H]'
            ]
            ace_matches = None
            for smarts_pattern in smarts_patterns:
                # Find matches in the molecule
                matches = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts_pattern))
                # If matches are found, assign them to ace_matches and break the loop
                if matches:
                    ace_matches = matches
                    break
            if ace_matches is None:
                print("No N-cyclized ACE matches were found.") 
            else:
                print("N-cyclized ACE Matches were found")
                n_cyclized = True
            
        elif len(ace_matches) == 1:
            print("Warning: ACE atoms is N-ethylated ACE")
    elif len(ace_matches) == 1:
        print("Warning: ACE atoms is N-methylated ACE")


# assign
#for index in ace_matches:
#    if index in not in M__POLY_IGNORE, M__POLY_LOWER, M__POLY_N_BB, M_POLY




elif len(ace_matches) > 1:
    raise ValueError("Multiple ACE groups found. Expected only one ACE group.")
if len(nme_matches) == 0:
    raise ValueError("No NME group found. Expected one NME group.")
if len(nme_matches) > 1:
    raise ValueError("Multiple NME groups found. Expected only one NME group.")



  
# Print ACE atom information  ######################
ace_matches_custom_indices = []
# debug
# print("ACE atom information:")
for index in ace_matches[0]:
    atom = custom_atoms[index]
    # debug
    print(f"Atom Index: {atom.atom_index}, Atom Name: {atom.atom_name}, Element Name: {atom.atom_element}")
    ace_matches_custom_indices.append(atom.atom_index)



M__POLY_PROPERTIES = []
M__POLY_N_Alkylation = []
# check ACE matches, N-methylation, N-ethylation
ace_atom_sequence = ["O", "C", "N", "H", "C", "H", "H", "H"]
ace_match_indices = ace_matches[0]
sequence_matches = all(custom_atoms[index].atom_element == ace_atom_sequence[i] for i, index in enumerate(ace_match_indices))
if not sequence_matches:
    print("Warning: ACE atoms not matches! Now try matching N-methylated ACE")
    ace_atom_sequence = ["O", "C", "N", "C", "C", "H", "H", "H","C", "H", "H", "H"]
    ace_match_indices = ace_matches[0]
    sequence_matches = all(custom_atoms[index].atom_element == ace_atom_sequence[i] for i, index in enumerate(ace_match_indices))
    if not sequence_matches:
        print("Warning: ACE atoms not matches! Now try matching N-ethylated ACE")
        ace_atom_sequence = ["O", "C", "N", "C", "C", "H", "H", "C", "H", "H", "H","C", "H", "H", "H"]
        ace_match_indices = ace_matches[0]
        sequence_matches = all(custom_atoms[index].atom_element == ace_atom_sequence[i] for i, index in enumerate(ace_match_indices))
        if not sequence_matches:
            if not ace_matches:
                raise ValueError("The ACE group does not match with any pattern including N-alkylated")
            else:
                print("Warning: ACE atoms likey matched with a N-cyclized ACE")
                M__POLY_PROPERTIES.append("N-CYCLIZATION")  
        else: 
            print("Warning: ACE atoms match with N-ethylated ACE")        
            M__POLY_N_Alkylation.extend([custom_atoms[ace_match_indices[i]].atom_index for i in [4]])
            M__POLY_N_Alkylation.extend([custom_atoms[ace_match_indices[i]].atom_index for i in [7]])
            M__POLY_PROPERTIES.append("N-ETHYLATION")  
    else:
        print("Warning: ACE atoms match with N-methylated ACE")
        M__POLY_N_Alkylation.extend([custom_atoms[ace_match_indices[i]].atom_index for i in [4]])
        M__POLY_PROPERTIES.append("N-METHYLATION")  








# Print NME atom information #######################
nme_matches_custom_indices = []
# print("NME atom information:")
for index in nme_matches[0]:
    atom = custom_atoms[index]
    # print(f"Atom Index: {atom.atom_index}, Atom Name: {atom.atom_name}, Element Name: {atom.atom_element}")
    nme_matches_custom_indices.append(atom.atom_index)

# check NME matches
nme_atom_sequence = ["O", "C", "N", "H", "C", "H", "H", "H"]
nme_match_indices = nme_matches[0]
nme_sequence_matches = all(custom_atoms[index].atom_element == nme_atom_sequence[i] for i, index in enumerate(nme_match_indices))
if not nme_sequence_matches:
    raise ValueError("The NME group does not match the expected sequence.")


# assign each NME and ACE atom to the corresponding variable

M__POLY_IGNORE = []


# for NME
M__POLY_UPPER = custom_atoms[nme_match_indices[2]].atom_index
M__POLY_IGNORE.extend([custom_atoms[nme_match_indices[i]].atom_index for i in [3, 4, 5, 6, 7]])

M__POLY_C_BB = custom_atoms[nme_match_indices[1]].atom_index
M__POLY_O_BB = custom_atoms[nme_match_indices[0]].atom_index

# for ACE, now the N-H has been writen in for future usage for patches
M__POLY_LOWER = custom_atoms[ace_match_indices[1]].atom_index
M__POLY_IGNORE.extend([custom_atoms[ace_match_indices[i]].atom_index for i in [0, -1, -2, -3, -4]])
M__ROOT = custom_atoms[ace_match_indices[2]].atom_index
M__POLY_N_BB = custom_atoms[ace_match_indices[2]].atom_index

if M__POLY_N_Alkylation is []:
    M__POLY_H_BB = custom_atoms[ace_match_indices[3]].atom_index
else:
    M__POLY_H_BB = None




################ now handle the mainchain###########

# note: the =O of NME and C=O of ACE are counted too
main_chain = None
main_chain_length = sys.maxsize  # Initialize with the maximum possible value
main_chain_atoms = []

for start in nme_matches:
    for end in ace_matches:
        # Find the shortest path between the two groups
        path = Chem.rdmolops.GetShortestPath(mol, start[0], end[0])
        path_length = len(path) - 1
        
        # Update the shortest path if the current path is shorter
        if path_length < main_chain_length:
            main_chain = (start, end)
            main_chain_length = path_length
            main_chain_atoms = path


mainchain_matches_custom_indices = []
for index in main_chain_atoms:
    atom = custom_atoms[index]
    # print(f"Mainchain Atom Index: {atom.atom_index}, Mainchain Atom Name: {atom.atom_name}")
    mainchain_matches_custom_indices.append(atom.atom_index)

print("mainchain index:",mainchain_matches_custom_indices)



def mainchain_sequence_match(sequence, main_chain_atoms, custom_atoms):
    return all(custom_atoms[index].atom_element == sequence[i] for i, index in enumerate(main_chain_atoms))

# Define the specific sequence to match
gamma_AA_mainchain_sequence = ["O", "C", "C", "C", "C", "N", "C", "O"]
beta_AA_mainchain_sequence = ["O", "C", "C", "C", "N", "C", "O"]
alpha_AA_mainchain_sequence = ["O", "C", "C", "N", "C", "O"]


# Check mainchain matches
if mainchain_sequence_match(gamma_AA_mainchain_sequence, main_chain_atoms, custom_atoms):
    print("Warning: The mainchain atoms sequence matches a gamma amino acid")
    M__POLY_CA_BB = custom_atoms[main_chain_atoms[2]].atom_index
    M__POLY_CB_BB = custom_atoms[main_chain_atoms[3]].atom_index
    M__POLY_CG_BB = custom_atoms[main_chain_atoms[4]].atom_index
    AA_type = "GAMMA_AA"
elif mainchain_sequence_match(beta_AA_mainchain_sequence, main_chain_atoms, custom_atoms):
    print("Warning: The mainchain atoms sequence matches a beta amino acid")
    M__POLY_CA_BB = custom_atoms[main_chain_atoms[2]].atom_index
    M__POLY_CB_BB = custom_atoms[main_chain_atoms[3]].atom_index
    AA_type = "BETA_AA"
elif mainchain_sequence_match(alpha_AA_mainchain_sequence, main_chain_atoms, custom_atoms):
    print("Warning: The mainchain atoms sequence matches an alpha amino acid")
    M__POLY_CA_BB = custom_atoms[main_chain_atoms[2]].atom_index
    AA_type = "ALPHA_AA"
else:
    raise ValueError("The mainchain atoms sequence does not match with any of gamma, beta, or alpha amino acids")


###################### write the N-cyclized ACE atoms to M__POLY_N_Alkylation  then change the mol2 file    #################
print ("now handling N-cyclization...")

#print ("ace_matches_custom_indices:",ace_matches_custom_indices)
#print ("mainchain_matches_custom_indices:",mainchain_matches_custom_indices)

M__POLY_N_Alkylation = [atom for atom in ace_matches_custom_indices 
                        if atom not in mainchain_matches_custom_indices 
                        and atom != M__POLY_LOWER 
                        and atom != M__POLY_H_BB
                        and atom not in M__POLY_IGNORE]

print ("M__POLY_H_BB:",M__POLY_H_BB)

print ("M__POLY_N_Alkylation:",M__POLY_N_Alkylation)

M__POLY_N_BONDED_ATOM = None

# Iterate over custom_bonds to find the N-bonded atom
for custom_bond in custom_bonds:
    # Check if atom indices match the specified criteria
    if (custom_bond.bonding_atom1_index in M__POLY_N_Alkylation and custom_bond.bonding_atom2_index == M__POLY_N_BB) or \
       (custom_bond.bonding_atom2_index in M__POLY_N_Alkylation and custom_bond.bonding_atom1_index == M__POLY_N_BB):
        print(f"Bonding information: Bond index - {custom_bond.bond_index}, Atom1 index - {custom_bond.bonding_atom1_index}, Atom2 index - {custom_bond.bonding_atom2_index}, Bond type - {custom_bond.bond_type}")
        # Store the atom index from M__POLY_N_Alkylation
        M__POLY_N_BONDED_ATOM = custom_bond.bonding_atom1_index if custom_bond.bonding_atom1_index in M__POLY_N_Alkylation else custom_bond.bonding_atom2_index
        print("M__POLY_N_BONDED_ATOM:",M__POLY_N_BONDED_ATOM)
        break  # Break the loop after finding the first match
# At this point, 'M__POLY_N_BONDED_ATOM' holds the atom index from M__POLY_N_Alkylation that is bonded to M__POLY_N_BB, or None if no such bond was found.


# do create and append the virtual atom, that is a duplciation Nbb
# Might give bugs for N-methylation, ethylation etc, will see
if M__POLY_N_BONDED_ATOM is not None and n_cyclized == True:
    # Find the line with the atom information for M__POLY_N_BONDED_ATOM
    for line in mol2_lines[atom_section_start:atom_section_end]:
        fields = line.split()
        #print("fields:",fields)
        if len(fields) > 0 and fields[0].isdigit() and int(fields[0]) == M__POLY_N_BONDED_ATOM:
            atom_info_line = line
            break
    
    # Modify the atom information line
    fields = atom_info_line.split()
    fields[0] = str(int(custom_atoms[-1].atom_index + 1))
    virtual_atomX_index = fields[0]
    print ("virtual_atomX_index:",virtual_atomX_index)
    fields[1] = '   X'  # Change the atom name to 'X'
    fields[5] = '   X'  # Change the atom type to 'X'
    fields[-1] = '0.000'  # Change the partial charge to '0.000'
    new_atom_info_line = ' '.join(fields)
    
    # Append the modified atom information line to the atom section
    mol2_lines.insert(atom_section_end, new_atom_info_line)

    # increase bond section + 1 since new line inserted
    bond_section_start+=1
    bond_section_end+=1

    # Combine the lines back into a single string
    new_mol2_contents = '\n'.join(mol2_lines)




# append the virtual bond information
N_bonded_line_index = None
virtual_bond_index = None

if M__POLY_N_BONDED_ATOM is not None and n_cyclized == True:
    # Find the line with the bond information for the bond between M__POLY_N_BONDED_ATOM and M__POLY_N_BB
    for i, line in enumerate(mol2_lines[bond_section_start:bond_section_end], bond_section_start):
        fields = line.split()
        print("fields:",fields)
        if len(fields) > 0 and fields[0].isdigit():
            bonding_atom1_index = int(fields[1])
            bonding_atom2_index = int(fields[2])
            
            if (bonding_atom1_index == M__POLY_N_BONDED_ATOM and bonding_atom2_index == M__POLY_N_BB) or \
               (bonding_atom2_index == M__POLY_N_BONDED_ATOM and bonding_atom1_index == M__POLY_N_BB):
                N_bonded_line_index = i
                virtual_bond_index = int(fields[0])
                break

    # Create a new bond line
    fields = mol2_lines[N_bonded_line_index].split()
    fields[0] = str(int(custom_bonds[-1].bond_index + 1))
    if bonding_atom1_index == M__POLY_N_BB:
        fields[1] = virtual_atomX_index
    else:
        fields[2] = virtual_atomX_index
    new_bond_line = ' '.join(fields)

    # Append the modified bond information line to the atom section
    mol2_lines.insert(bond_section_end, new_bond_line)

    # comment out the original N-bonded line
    mol2_lines[N_bonded_line_index] = '#' + mol2_lines[N_bonded_line_index]

    # Combine the lines back into a single string
    new_mol2_contents = '\n'.join(mol2_lines)
    print("new_mol2_contents:",new_mol2_contents)

    # Derive new filename
    # new_filename = mol2_file.rsplit('.', 1)[0] + "_virt.mol2"

    # Write to new file
    with open(filename, 'w') as f:
        f.write(new_mol2_contents)


#print ("M__POLY_N_Alkylation:",M__POLY_N_Alkylation)
#print ("M__POLY_N_BB:",M__POLY_N_BB)
##################################################   N-cyclization modification END ############################



############## find out properties #############

#M__POLY_PROPERTIES = []
M__POLY_PROPERTIES.append(AA_type)


if filename.startswith("L"):
    Chiral_type = "L_AA"
elif filename.startswith("D"):
    Chiral_type = "D_AA"
else:
    Chiral_type = "ACHIRAL_BACKBONE"


M__POLY_PROPERTIES.append(Chiral_type)


#print ("the ace custom indices:",ace_matches_custom_indices)
#print ("the nme custom indices:",nme_matches_custom_indices)
#print ("the mainchain custom indices:",mainchain_matches_custom_indices)
sidechain_heavy_atoms = []
sidechain_heavy_atom_custom_indices = []
#these includes the side chain heavy atoms. It will include mainchain Hs and H-N if not specify the element =! H
for atom in custom_atoms:
    if atom.atom_element != "H" and all(atom.atom_index not in match for match in (ace_matches_custom_indices, nme_matches_custom_indices, mainchain_matches_custom_indices)):
        # print(f" Sidechain Heavy Atom Index: {atom.atom_index}, Sidechain Heavy Atom Name: {atom.atom_name}")
        sidechain_heavy_atom_custom_indices.append(atom.atom_index)
        sidechain_heavy_atoms.append(atom)

sidechain_all_atoms = []
sidechain_all_atom_custom_indices = []
for atom in custom_atoms:
    if all(atom.atom_index not in match for match in (ace_matches_custom_indices, nme_matches_custom_indices, mainchain_matches_custom_indices)):
        # print(f" Sidechain All Atom Index: {atom.atom_index}, Sidechain All Atom Name: {atom.atom_name}")
        # print("Warning: backbone CHs are included too")
        sidechain_all_atom_custom_indices.append(atom.atom_index)
        sidechain_all_atoms.append(atom)

# print ("sidechain_all_atom_custom_indices:",sidechain_all_atom_custom_indices)


Polarity = "HYDROPHOBIC"  # Initial assumption
for atom in sidechain_all_atoms:
    if float(atom.atom_partial_charge) > 0.3:
        print("Warning: a polar hydrogen found")
        Polarity = "POLAR"
        break  # Once we find a polar atom, we can stop checking

M__POLY_PROPERTIES.append(Polarity)


Aliphatic = "ALIPHATIC"  # Start by assuming all atoms are aliphatic
for atom in sidechain_heavy_atoms:
    if atom.atom_element not in ["C", "H"] and atom.atom_type not in ["C.ar"]:
        Aliphatic = None  # If any atom is not "C" or "H", set Aliphatic to None
        break  # No need to continue checking other atoms

M__POLY_PROPERTIES.append(Aliphatic)


Metal = "METALBINDING"
M__POLY_PROPERTIES.append(Metal)




Aromaticity = ""
# for heteroaromatic cases
for atom in sidechain_heavy_atoms:
    if atom.atom_type in ["C.ar"] or (atom.atom_type in ["N.pl3"] and float(atom.atom_partial_charge) > -0.4) :
        Aromaticity = "AROMATIC"
        break  # Once we find aromatic carbon in the sidechain atoms, then it is aromatic

# for the Tyr and Phe cases
ar_count = sum(1 for bond in custom_bonds if bond.bond_type == "ar")
if ar_count > 4:
    Aromaticity = "AROMATIC"




M__POLY_PROPERTIES.append(Aromaticity)


Charged = ""

if M__POLY_CHG != 0:
    Charged = "CHARGED"
    print("Warning: the molecule is charged")
M__POLY_PROPERTIES.append(Charged)


Charged_status = ""
if M__POLY_CHG > 0:
    Charged_status = "POSITIVE_CHARGE"
    print("Warning: the molecule is positively charged")

if M__POLY_CHG < 0:
    Charged_status = "NEGATIVE_CHARGE"
    print("Warning: the molecule is negatively charged")

M__POLY_PROPERTIES.append(Charged_status)



#############write the M block################


# If there is an existing M session
if "M  END\n" in lines:
    # If the --overwrite option is present
    if args.overwrite:
        # Read the mol2 file and keep lines not starting with "M  "
        with open(args.mol2file, 'r') as f:
            lines = [line for line in f if not line.startswith("M  ")]
        # Write the lines back to the mol2 file
        with open(args.mol2file, 'w') as f:
            f.writelines(lines)
    else:
        print("Warning: The mol2 file already has a M session indicated by 'M  END'. Not writing anything")
        sys.exit()


print("Now writing M session..")
with open(filename, "a") as f:
    f.write("M  ROOT {}\n".format(M__ROOT))
    f.write("M  POLY_N_BB {}\n".format(M__POLY_N_BB))
    try:
        f.write("M  POLY_CG_BB {}\n".format(M__POLY_CG_BB))
    except NameError:
        print ("Warning: M  POLY_CG_BB has no name, it is not a GAMMA_AA")
    try:
        f.write("M  POLY_CB_BB {}\n".format(M__POLY_CB_BB))
    except NameError:
        print ("Warning: M  POLY_CB_BB has no name, it is not a BETA_AA")
    f.write("M  POLY_CA_BB {}\n".format(M__POLY_CA_BB))
    f.write("M  POLY_C_BB {}\n".format(M__POLY_C_BB))
    f.write("M  POLY_O_BB {}\n".format(M__POLY_O_BB))
    f.write("M  POLY_H_BB {}\n".format(M__POLY_H_BB))
    f.write("M  POLY_IGNORE {}\n".format(' '.join(map(str, M__POLY_IGNORE))))
    f.write("M  POLY_UPPER {}\n".format(M__POLY_UPPER))
    f.write("M  POLY_LOWER {}\n".format(M__POLY_LOWER))
    f.write("M  POLY_N-ALKYLATION {}\n".format(' '.join(map(str, M__POLY_N_Alkylation))))
    f.write("M  POLY_N_BONDED_ATOM {}\n".format(M__POLY_N_BONDED_ATOM)) 
    f.write("M  POLY_CHG {}\n".format(int(round(M__POLY_CHG))))
    f.write("M  POLY_SMILES {}\n".format(M__POLY_SMILES))
    f.write("M  POLY_PROPERTIES PROTEIN {}\n".format(' '.join(str(x) for x in M__POLY_PROPERTIES if x is not None)))
    f.write("M  END\n")


