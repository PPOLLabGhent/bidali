# -*- coding: utf-8 -*-
"""HGNC genenames

Reference: https://www.genenames.org/
"""
from urllib import error,parse
from urllib.request import Request,urlopen
import json

# Get the whole table
from bidali.LSD.dealer.external.ensembl import get_genenames

# genenames API
def fetchGenenamesJSON(symbol,alias=False):
    """
    Requires a gene symbol
    """
    if not alias:
        requesturl = 'http://rest.genenames.org/fetch/symbol/{}'
    else:
        requesturl = 'http://rest.genenames.org/fetch/alias_symbol/{}'
    requesturl = requesturl.format(parse.quote(symbol))
    request = Request(requesturl,
                      headers={'Accept':'application/json'})
    response = urlopen(request)
    content = response.read()
    if content:
        return json.loads(content.decode())
    
def fetchAliases(symbol,unknown_action='raise'):
    """
    Returns list of all aliases for a symbol.
    First alias in list is official symbol.

    Args:
        symbol (str): Symbol name.
        unknown_action (option):
          'raise': raise KeyError
          'list': [symbol]
          'none': None
    """
    symbolJSON = fetchGenenamesJSON(symbol)
    if not symbolJSON['response']['numFound']:
        symbolJSON = fetchGenenamesJSON(symbol,alias=True)
    try:
        symbols = [symbolJSON['response']['docs'][0]['symbol']
        ]+symbolJSON['response']['docs'][0]['alias_symbol']
        return symbols
    except KeyError:
        #Not a recognised symbol
        if unknown_action == 'raise': raise
        elif unknown_action == 'list': return [symbol]
        else: return None
