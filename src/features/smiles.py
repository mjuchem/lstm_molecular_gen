import numpy as np
from rdkit import Chem


class SmilesTokenizer(object):
    def __init__(self):
        atoms = ['Li', 'Na', 'Al', 'Si', 'Cl', 'Sc', 'Zn', 'As', 'Se', 'Br', 'Sn', 'Te', 'Cn', 'H', 'B', 'C', 'N', 'O', 'F', 'P', 'S', 'K', 'V', 'I', ]
        special = ['(', ')', '[', ']', '=', '#', '%', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '-', 'se', 'te', 'c', 'n', 'o', 's']
        padding = ['G', 'A', 'E']

        self.table = sorted(atoms, key=len, reverse=True) + special + padding
        self.table_len = len(self.table)

        self.table_2_chars = list(filter(lambda x: len(x) == 2, self.table))
        self.table_1_chars = list(filter(lambda x: len(x) == 1, self.table))

        self.one_hot_dict = {}
        for i, symbol in enumerate(self.table):
            vec = np.zeros(self.table_len, dtype=np.float32)
            vec[i] = 1
            self.one_hot_dict[symbol] = vec

    def tokenize(self, smiles):

        smiles = smiles + ' '
        
        N = len(smiles)
        
        token = []
        i = 0
        
        while (i < N):
            c1 = smiles[i]
            c2 = smiles[i : i+2]
            
            if (c2 in self.table_2_chars):
                token.append(c2)
                i = i + 1
                continue
                
            if (c1 in self.table_1_chars):
                token.append(c1)
                i = i + 1
                continue
                
            i = i + 1

        return token

    def one_hot_encode(self, tokenized_smiles):
        result = np.array(
            [self.one_hot_dict[symbol] for symbol in tokenized_smiles],
            dtype=np.float32)
        result = result.reshape(1, result.shape[0], result.shape[1])
        return result

    def embeddings(self, tokenized_smiles):
        result = [self.table.index(symbol) for symbol in tokenized_smiles]
        return result

    def embeddings_to_tokens(self, embeddings):
        embeddings = embeddings[embeddings != self.zero()]
        result = [self.table[s] for s in embeddings]
        return result

    def embeddings_to_smiles(self, embeddings):
        tokens = self.embeddings_to_tokens(embeddings)
        s = ""
        for t in tokens:
            s = s + t
        return s
    
    def zero(self):
        return self.table.index('A')


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
