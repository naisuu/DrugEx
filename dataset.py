import os
import pandas as pd
from rdkit import Chem
from rdkit import rdBase
from tqdm import tqdm
import gzip
import utils

rdBase.DisableLog('rdApp.info')


def corpus(input, output, is_sdf=False, requires_clean=True, is_isomerice=False):
    """ Constructing dataset with SMILES-based molecules, each molecules will be decomposed
        into a series of tokens. In the end, all the tokens will be put into one set as vocaulary.

        Arguments:
            input (string): The file path of input, either .sdf file or tab-delimited file

            output (string): The file path of output

            is_sdf (bool): Designate if the input file is sdf file or not

            requires_clean (bool): If the molecule is required to be clean, the charge metal will be
                    removed and only the largest fragment will be kept.

            is_isomerice (bool): If the molecules in the dataset keep conformational information. If not,
                    the conformational tokens (e.g. @@, @, \, /) will be removed.

    """
    if is_sdf:
        # deal with sdf file with RDkit
        inf = gzip.open(input)
        fsuppl = Chem.ForwardSDMolSupplier(inf)
        df = []
        for mol in fsuppl:
            try:
                df.append(Chem.MolToSmiles(mol, is_isomerice))
            except:
                print(mol)
    else:
        # deal with table file
        df = pd.read_table(input).Smiles.dropna()
    voc = utils.Voc()
    words = set()
    canons = []
    tokens = []
    if requires_clean:
        smiles = set()
        for smile in tqdm(df):
            try:
                smile = utils.clean_mol(smile, is_isomeric=is_isomerice)
                smiles.add(Chem.CanonSmiles(smile))
            except:
                print('Parsing Error:', smile)
    else:
        smiles = df.values
    for smile in tqdm(smiles):
        token = voc.tokenize(smile)
        # Only collect the organic molecules
        if {'C', 'c'}.isdisjoint(token):
            print('Warning:', smile)
            continue
        # Remove the metal tokens
        if not {'[Na]', '[Zn]'}.isdisjoint(token):
            print('Redudent', smile)
            continue
        # control the minimum and maximum of sequence length.
        if 10 < len(token) <= 100:
            words.update(token)
            canons.append(smile)
            tokens.append(' '.join(token))

    # output the vocabulary file
    log = open(output + '_voc.txt', 'w')
    log.write('\n'.join(sorted(words)))
    log.close()

    # output the dataset file as tab-delimited file
    log = pd.DataFrame()
    log['Smiles'] = canons
    log['Token'] = tokens
    log.drop_duplicates(subset='Smiles')
    log.to_csv(output + '_corpus.txt', sep='\t', index=False)


def _download_chembl_file():
    """
    Support function that downloads the necessary ChEMBL data file "chembl_26.sdf.gz" to the data folder.
    """
    import requests

    url = 'https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_26/chembl_26.sdf.gz'
    filename = url.split('/')[-1]
    with requests.get(url, stream=True, verify=False) as response:
        response.raise_for_status()
        file_size = int(response.headers.get('content-length', 0))
        pbar = tqdm(total=file_size, desc=f'ChEMBL data: {filename}', unit='B', unit_scale=True)
        with open(f'data/{filename}', 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
                pbar.update(len(chunk))
        pbar.close()


if __name__ == '__main__':
    if not os.path.exists('data/chembl_26.sdf.gz'):
        try:
            _download_chembl_file()
        except Exception as e:  # TODO: Expand exception clause after testing edge-cases
            print(f'Something went wrong: {e}')
            print('Try downloading the dataset manually from the following url:'
                  ' https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_26/chembl_26.sdf.gz')

    corpus('data/chembl_26.sdf.gz', 'data/chembl', is_sdf=True)
    corpus('data/LIGAND_RAW.tsv', 'data/ligand', is_sdf=False)
    # corpus('data/guacamol.smiles', 'data/guacamol', requires_clean=True, is_sdf=False)
