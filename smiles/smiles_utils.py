from rdkit import Chem
from smiles.smiles_tokenizer import SmilesTokenizer


def cleanup_list_smiles(list_of_smiles):
    valid_smiles = []
    for smi in list_of_smiles:
        mol = Chem.MolFromSmiles(smi)
        if (mol is not None):
            valid_smiles.append(Chem.MolToSmiles(mol))
    return valid_smiles

def validate_mols(list_of_smiles):
    valid_mols = []
    for smi in list_of_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            valid_mols.append(mol)
    return valid_mols

def convert_mols_to_smiles(list_of_mols):
    valid_smiles = [Chem.MolToSmiles(mol) for mol in list_of_mols]
    return valid_smiles

def encode_list_smiles(list_of_smiles):
    st = SmilesTokenizer()
    encoded_smiles = []
    for smi in list_of_smiles:
        token = st.smiles_to_tokens(smi)
        encoded_smiles.append(st.tokens_to_embeddings(token))
    return encoded_smiles
