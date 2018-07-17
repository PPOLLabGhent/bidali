# -*- coding: utf-8 -*-
"""Ontologies from BioPortal

Reference: https://bioportal.bioontology.org/
"""
from bidali import LSD
from bidali.config import secrets
from bidali.LSD import retrieveSources,cacheableTable,processedDataStorage,datadir
import os, gzip, pandas as pd
from io import TextIOWrapper, StringIO

def list_ontologies():
    import urllib.request, urllib.error, urllib.parse
    import json
    import os
    from pprint import pprint
    
    REST_URL = "http://data.bioontology.org"
    API_KEY = secrets.get_secret() #need to register -> then https://bioportal.bioontology.org/account API key
    
    def get_json(url):
        opener = urllib.request.build_opener()
        opener.addheaders = [('Authorization', 'apikey token=' + API_KEY)]
        return json.loads(opener.open(url).read())
    
    # Get the available resources
    resources = get_json(REST_URL + "/")
    
    # Get the ontologies from the `ontologies` link
    ontologies = get_json(resources["links"]["ontologies"])
    
    # Get the name and ontology id from the returned list
    ontology_output = []
    for ontology in ontologies:
        ontology_output.append(f"{ontology['name']}\n{ontology['@id']}\n")
    
    # Print the first ontology in the list
    pprint(ontologies[0])
    
    # Print the names and ids
    print("\n\n")
    for ont in ontology_output:
        print(ont)

