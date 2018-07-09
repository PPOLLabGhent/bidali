# -*- coding: utf-8 -*-
"""Publication crawler and text mining module


Functions for creating dictionaries and starting up text mining explorations.
"""
import pandas as pd, numpy as np

# Utilities for making dictionaries
##CIDcounter class to create CID (concept id) according to Adil's requirement that it
##is the same as TID reference symbol
class CIDcounter:
    def __init__(self,startvalue=0):
        import itertools as it
        self.counter = it.count(startvalue)
        
    def next(self,increment):
        if increment:
            self.value = next(self.counter)
        else: next(self.counter) #this step ensures a CID == TID for a new reference concept
        return self.value

def makeCIDandTID(df,conceptCol='symbol',aliasCol='alias'):
    """Make text mining dictionary with TID and CID set

    According to Adil Salhi's text mining requirements.

    Args:
        df (pd.DataFrame): Dataframe to transform.
        conceptCol (str): Column name of main concept/symbol.
        aliasCol (str): Column name where all aliases for the concept
          are gathered in a list, including the concept itself as first list member.
    """
    from bidali.util import unfoldDFlistColumn
    df = unfoldDFlistColumn(df,aliasCol)
    df.reset_index(inplace=True,drop=True)
    df.index.name = 'TID' #term id
    c = CIDcounter()        
    df.insert(0, 'CID', (~df[conceptCol].duplicated()).apply(c.next)) #need to work with `not` duplicated column!
    return df

## get icd9 info
def get_icd9info(icd9code,type='disease',includeShortName=True):
    """Get info related to icd9code

    Info:
        https://clinicaltables.nlm.nih.gov/

    Reference:
        https://clinicaltables.nlm.nih.gov/apidoc/icd9cm_dx/v3/doc.html

    Args:
        icd9code (str): Code to query.
        type (option): `disease` or `procedure`.
        includeShortName: Also include the short name in output.
    """
    import requests, json
    if type == 'disease':
        url = 'https://clinicaltables.nlm.nih.gov/api/icd9cm_dx/v3/search?terms={term}'
    elif type == 'procedure':
        url = 'https://clinicaltables.nlm.nih.gov/api/icd9cm_sg/v3/search?terms={term}'
    else: raise Exception('Wrong type specified, should be disease or procedure, not',type)
    if includeShortName: url+='&ef=short_name'
    r = requests.get(url.format(term = icd9code))
    return json.loads(r.content)
    
# Dictionary functions
def get_biomart_gene_dictionary():
    from bidali.LSD.dealer.external.ensembl import get_biomart
    from bidali.genenames import fetchAliases
    bm = get_biomart(atts=[
        'ensembl_gene_id',
        'ensembl_peptide_id',
        'pdb',
        'entrezgene',
        'hgnc_symbol',
    ])
    HGNC_aliases = bm['HGNC symbol'].apply(lambda x: [(a,x) for a in fetchAliases(x,unknown_action='list')])
    HGNC_aliases = pd.concat([pd.Series(dict(a)) for a in HGNC_aliases])
    HGNC_aliases = pd.Series(HGNC_aliases.index.values, index=HGNC_aliases)
    HGNC_aliases.index.set_names('HGNC symbol', inplace = True)
    HGNC_aliases.name = 'HGNC alias'
    return bm.join(other = HGNC_aliases, on = 'HGNC symbol')

def get_gene_dictionary():
    """Gene symbol dictionary.
    Uses genenames.org data.
    """
    from bidali.LSD.dealer.external import ensembl
    from bidali.util import unfoldDFlistColumn
    gn = ensembl.get_genenames()
    gn['alias'] = gn.T.apply(
        lambda x: [x.symbol, x['name']] + #optionally add other names such as protein names here
        (
            x.alias_symbol.split('|') if x.alias_symbol is not np.nan else []
        )
    )
    print('Original columns',gn.columns)
    gn = gn[['hgnc_id','symbol','alias','uniprot_ids']].copy()
    print('Columns ketp:',gn.columns)
    gn = makeCIDandTID(gn)
    return gn
