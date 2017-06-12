import sys
from urllib import error,parse
import urllib.request
import json
import time

class EnsemblRestClient:
    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, headers=None, params=None):
        if headers is None:
            headers = {}

        if 'Content-Type' not in headers:
            headers['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + parse.urlencode(params)

        data = None
        if self.req_count >= self.reqs_per_sec: # check if we need to rate limit ourselves
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0
        try:
            request = urllib.request.Request(self.server + endpoint, headers=headers)
            response = urllib.request.urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content.decode())
            self.req_count += 1
        except error.HTTPError as e:
            if e.code == 429: # check if we are being rate limited by the server
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, headers, params)
            else:
                sys.stderr.write('Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))
        return data
    
    def get_variants(self, species, symbol):
        genes = self.perform_rest_action(
            '/xrefs/symbol/{0}/{1}'.format(species, symbol), 
            params={'object_type': 'gene'}
        )
        if genes:
            stable_id = genes[0]['id']
            variants = self.perform_rest_action(
                '/overlap/id/{0}'.format(stable_id),
                params={'feature': 'variation'}
            )
            return variants
        return None

def run(species, symbol):
    """
    Example:
    run('human', 'BRAF')
    """
    client = EnsemblRestClient()
    variants = client.get_variants(species, symbol)
    if variants:
        for v in variants:
            print('{seq_region_name}:{start}-{end}:{strand} ==> {id} ({consequence_type})'.format(**v))

def getSequences(geneset,file,expand_5prime=0,sequenceEnd=-1):
    """
    getSequences expects a list of gene symbols, to retrieve sequences
    and a file name, to which the sequences are stored for future
    requests
    expand_5prime: include x bp upstream
    sequenceEnd: only keep x bp; -1 -> whole sequence

    Raises exception if symbols in file do not match requested symbols
    TODO also raise if other params do not match
    """
    #Settings
    currentSettings = (sorted(geneset),expand_5prime,sequenceEnd)
    
    #Imports
    import pickle

    #Setup client
    client = EnsemblRestClient()

    #Get sequences
    try:
        settings,sequences = pickle.load(open(file,'rb'))
        if settings != currentSettings:
            raise Exception('File content and settings provided do not match')
    except FileNotFoundError:
        sequences = {}
        for gene in geneset:
            try:
                geneID = client.perform_rest_action('/xrefs/symbol/human/{}'.format(gene))[0]['id']
                sequence = client.perform_rest_action('/sequence/id/{}'.format(geneID),
                                                      params={'expand_5prime':expand_5prime})['seq']
                sequences[gene] = sequence[:sequenceEnd]
            except IndexError: pass
        pickle.dump((currentSettings,sequences),open(file,'wb'))
    return sequences

def getOrthologs(symbol,species='human',ospecies='mus_musculus'):
    #Setup client
    client = EnsemblRestClient()
    
    data = client.perform_rest_action('/homology/symbol/{}/{}'.format(
        species,parse.quote(symbol)))
    retrieved_homologies = []
    for homologies in data['data']:
        id = homologies['id']
        for homology in homologies['homologies']:
            if ospecies in homology['target']['species']:
                retrieved_homologies.append(homology)
    osymbols = []
    for homology in retrieved_homologies:
        data = client.perform_rest_action('/xrefs/id/{}'.format(
            parse.quote(homology['target']['id'])))
        for db in data:
            if db['dbname'] == 'EntrezGene':
                osymbols.append(db['display_id'])
                break
    return osymbols
