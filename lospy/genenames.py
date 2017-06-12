from urllib import error,parse
from urllib.request import Request,urlopen
import json

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
    
def fetchAliases(symbol):
    """
    Returns list of all aliases for a symbol.
    First alias in list is official symbol.
    """
    symbolJSON = fetchGenenamesJSON(symbol)
    if not symbolJSON['response']['numFound']:
        symbolJSON = fetchGenenamesJSON(symbol,alias=True)
    symbols = [symbolJSON['response']['docs'][0]['symbol']
    ]+symbolJSON['response']['docs'][0]['alias_symbol']
    return symbols
