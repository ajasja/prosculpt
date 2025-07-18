"""Unit to calculate protein properties from protein sequence."""

import Bio
from Bio import SeqIO

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import Bio.SeqUtils as SeqUtils

# import truncator
import pandas as pd

kd_hydrophobicity = \
    {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
     'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
     'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
     'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

# taken form http://www.petercollingridge.co.uk/book/export/html/474
pKa = {'D': 3.9,
       'E': 4.3,
       'H': 6.1,
       'C': 8.3,
       'Y': 10.1,
       'K': 10.5,
       'R': 12,
       'N-terminus': 8,
       'C-terminus': 3.1}
charges = {'D': -1., 'E': -1., 'H': 1., 'C': -1., 'Y': -1., 'K': 1., 'R': 1., 'N-terminus': 1., 'C-terminus': -1.}


def amino_acid_charge(amino_acid, pH=7.5):
    aa_pKa = pKa.get(amino_acid, 0)
    if aa_pKa == 0:
        return 0
    ratio = 1. / (1. + 10 ** (pH - aa_pKa))

    if charges[amino_acid] == 1:
        return ratio
    else:
        return ratio - 1

standard_amino_acids = "ARNDBCEQZGHILKMFPSTWYV"
def get_charge_scale(pH=7.5):
    """Returns a dictionary for the charges of each standard amino acid at a given pH"""
    return {aa: amino_acid_charge(aa, pH=pH) for aa in standard_amino_acids}


def protein_charge(sequence, pH=7.5, blocked=False):
    protein_charge = 0
    if not blocked:
        protein_charge = amino_acid_charge('N-terminus', pH)
        protein_charge += amino_acid_charge('C-terminus', pH)

    for aa in list(pKa.keys()):
        protein_charge += sequence.count(aa) * amino_acid_charge(aa, pH)

    return protein_charge


def analyse_sequence(seq, name=None, pH=7.5, initial_dict=None):
    """Returns sequnce properties in a dictionary"""
    if initial_dict is None:
        res = {}
    else:
        res = initial_dict
    
    
    ana = ProteinAnalysis(seq)
    
    if not (name is None):
        res['name'] = name
    res['sequence'] = seq
    res['molecular_weight'] = SeqUtils.molecular_weight(seq, monoisotopic=False, seq_type='protein')
    res['molecular_weight_monoisotopic'] = SeqUtils.molecular_weight(seq, monoisotopic=True, seq_type='protein')
    res['extinction_280nm_reduced'] = ana.molar_extinction_coefficient()[0]
    res['extinction_280nm_cys_cys'] = ana.molar_extinction_coefficient()[1]
    res['Abs_1mg_ml_280nm_reduced'] = res['extinction_280nm_reduced'] / res['molecular_weight']
    res['Abs_1mg_ml_280nm_cys_cys'] = res['extinction_280nm_cys_cys'] / res['molecular_weight']
    res['isoelectric_point'] = ana.isoelectric_point()
    res[f'charge_pH{pH}'] = protein_charge(seq, pH)
    res['gravy'] = ana.gravy()
    
    
    return res

def extract_seq_from_vector(vector_file_name, label='GENE_PRODUCT'):
    """Loads a vector file, gets the DNA sequence and translates it protein. Returns a string"""
    vector =  Bio.SeqIO.read(vector_file_name, "genbank")
    vector.seq.alphabet=Bio.Alphabet.DNAAlphabet()
    gp = truncator.find_annotation(vector, label)
    full_sequence_AA = gp.extract(vector).seq.translate(cds=True)
    return str(full_sequence_AA)
    
#extract_seq_from_vector('06_DNA_vectors/ALAF02_A__pet28b+.gb')

#my_seq = "GSSQEEYVELLEQHERAVRELLRIAEEHKKGTENADELLRKLDTILDEAQKIIQTANKLLKESGSGTTEAIKRSEESVQQVETVIEIFQQSREKG"    
#analyse_sequence(my_seq, "Test")

def add_seq_fields(df, fields=None):
    """Adds the fields returned by analyse_sequence INPLACE"""
    if fields is None:
        fields = truncator.analyse_sequence('FAKEGFPY')
    for field in fields.keys():
        df[field] = fields[field]    
    return df

def fasta_to_dataframe(fasta, pH=7.5, ignore_seq='*EGDII'):
    """Coverts a fasta file to a data frame. Ignores the seq `ignore_seq` (for example from a ribosome binding site)"""
    from Bio import SeqIO
    fastas_dict = {}
    res = []
    for record in SeqIO.parse(fasta, "fasta"):
        seq = str(record.seq)
        seq = seq.replace(ignore_seq,'')
        if seq.endswith('**'):
            seq=seq[:-2]
        if seq.endswith('*'):
            seq=seq[:-1]
        fastas_dict[record.id]=seq
  
    
    for name in fastas_dict.keys():
        line = {}
        line['name'] = name
        line['sequence'] = fastas_dict[name]
        #print(fastas_dict[name])
        line = truncator.analyse_sequence(line['sequence'], pH=7.5, initial_dict=line)
        res.append(line)
    return pd.DataFrame(res, columns=res[0].keys())