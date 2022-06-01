import periodictable
import re
import json
import time
from rdkit import Chem

#################################### CLASS ####################################
class NotSupportError(Exception):
    pass
class CASNotValidError(Exception):
    pass

#################################### function definition ####################################
            
def CheckQueryType(query):
    if not query:
        raise ValueError('InputStringEmptyorNone')
        
    if re.match('^[0-9]+-[0-9]{2}-[0-9]{1}$',query):
        # CAS format checking
        l = query.split('-')
        validnum = int(l[2])
        tmp =(l[0]+l[1])[::-1]

        num = 0
        for idx, c in enumerate(tmp):
            num+=int(c)*(idx+1)
        
        if num%10 == validnum:
            # valid CAS
            return 'CAS'
        else:
            raise CASNotValidError('Invalid CAS: '+query)
    else:
        return 'Name'
    
def GetMWfromMolecularFormula(formula):
    try:
        periodictable.formula(formula).mass
        return periodictable.formula(formula).mass
    except:
        return None
    
def GetPropfromPUGVIEWJSON(cid,res,prop):
    # pass returned json into this function avoid keeping sending request
    retprop=None

    if prop.lower() == "smiles":
        for i in res:
            if i['TOCHeading'] == 'Names and Identifiers':
                for j in i['Section']:
                    if j['TOCHeading'] == "Computed Descriptors":
                        isosmi = cansmi = None
                        for k in j['Section']:
                            if k['TOCHeading'] == "Isomeric SMILES":
                                isosmi = k['Information'][0]['Value']['StringWithMarkup'][0]['String'].strip()
                            elif k['TOCHeading'] == "Canonical SMILES":
                                cansmi = k['Information'][0]['Value']['StringWithMarkup'][0]['String'].strip()
#                         return cansmi if not isosmi else isosmi
                        return [cansmi,isosmi]
                            
                        if retprop == None:
                            raise UnexpectedError('SMILES of cid '+str(cid)+' is not found!')
    elif prop.lower() == "cas":
        retprop=[]
        ll_0 = None
        for i in res:
            if i['TOCHeading'] == 'Names and Identifiers':
                for j in i['Section']:
                    # other identifiers may be not found
                    if j['TOCHeading'] == 'Other Identifiers':
                        ll_0 = j['Section']
                        break
                        
                # "other identifier" section not available
                if ll_0 == None:
                    return []
                
                
                ll_1 = None
                for k in ll_0:
                    if k['TOCHeading'] == "CAS":
                        ll_1 = k['Information']
                        break
                        
                # "CAS" section not available
                if ll_1 == None:
                    return []
                else:
                    # Record.Section[2].Section[2].Section[0].Information[1].Value.StringWithMarkup[0].String
                    for cas in ll_1:
                        retprop.append(cas['Value']['StringWithMarkup'][0]['String'].strip())
                    return retprop
    elif prop.lower() == 'mw':
        
        for i in res:
            if i['TOCHeading'] == 'Chemical and Physical Properties':
                for j in i['Section']:
                    if j['TOCHeading'] == "Computed Properties":
                        for k in j['Section']:
                            if k['TOCHeading'] == "Molecular Weight":
                                retprop = float(k['Information'][0]['Value']['StringWithMarkup'][0]['String'])
                                return retprop
    elif prop.lower() == 'synonyms':
        
        #Record.Section[1].Section[2].Section[0].Information[0].Value.StringWithMarkup[4]
        for i in res:
            if i['TOCHeading'] == 'Names and Identifiers':
                for j in i['Section']:
                    if j['TOCHeading'] == 'Synonyms':
                        for k in j['Section']:
                            if k['TOCHeading'] == 'Depositor-Supplied Synonyms':
                                for l in k['Information']:
                                    return [m['String'] for m in l['Value']['StringWithMarkup']]
                                

                                
def GetMolIdentifier(query,resolver,query_MW=None):
    import requests
    """
    valid Identifiers: Chemical Name, CAS
    Resolver: SciFindern, CommonChemistry, PubChem
    
    """
    
    
    # check input validity
    if (not query) | (not resolver):
        raise ValueError('Input is Empty or None!')
    else:
        if (not isinstance(query, str)) | (not isinstance(resolver, str)):
            raise ValueError('Input should be String!')
        else:
            query=query.strip()
            resolver=resolver.strip().lower()
    

    
    
    # check identifier type
    q_type = CheckQueryType(query)
    
    # start fetch data from websites
    if resolver == 'commonchemistry':
        
        #######################################################################################
        # CAS -> CAS,IsomericSMILES, MW
        # smiles may be missing but molecular formula wouldn't
        # so we apply periodtable to calculate MW from formula so as to compare to aspen database
        #######################################################################################
        if q_type == 'CAS':
            url = 'https://commonchemistry.cas.org/api/detail?cas_rn='+query
            res = requests.get(url)
            
            # search succeeded -> input query matches 1. CAS used currently / 2. deprecated
            if res.status_code == 200:
                ret_dict = res.json()
                # deprecated
                ret_cas = ret_dict['rn']
                if q_type in ret_dict['replacedRns']:
                    print('input CAS '+q_type+' is deprecated! Currently valid CAS '+ret_dict['rn']+' will be returned')
                    ret_cas = ret_dict['replacedRns']
                # current
                else:
                    # "smile" may be isomeric or canonical SMILES
                    if ret_dict['smile'].strip() == "":
                        ret_smi = None
                    else:
                        ret_smi = ret_dict['smile'].strip()
                        
                    try:
                        Chem.CanonSmiles(ret_smi)
                    except:
                        # invalid or None
                        ret_smi = None
                        
                    if ret_dict['molecularFormula'].strip() == "":
                        ret_mw = None
                    else:
                        #post-processing MolecularFormula
                        ret_dict['molecularFormula'] = ret_dict['molecularFormula'].replace('<sub>','').replace('</sub>','').replace('.','')
#                         ret_mw = GetMWfromMolecularFormula(ret_dict['molecularFormula'])
                        if ret_smi:
                            from rdkit.Chem import Descriptors
                            ret_mw = Descriptors.MolWt(Chem.MolFromSmiles(ret_smi))
                        else:
                            ret_mw = GetMWfromMolecularFormula(ret_dict['molecularFormula'])
                
                ret_cas = re.sub('[^0-9]','',ret_cas)
                ret_cas = "-".join([ret_cas[:-3],ret_cas[-3:-1],ret_cas[-1]])
                
                return [ret_cas, ret_smi, ret_mw]
            # CAS not found
            else: 
                return [None,None,None]
        else:
            raise NotSupportError('CommonChemistry does not accept Identifier other than CAS!')
            
    elif resolver == 'pubchem':
        
        import urllib
        from bs4 import BeautifulSoup, SoupStrainer
        
        ########################################################################################
        # CAS/name -> CAS,
        # leverage https://www.ncbi.nlm.nih.gov/pccompound/?term= to find cid
        # unlike https://pubchem.ncbi.nlm.nih.gov/#query=, matched results are not rendered by javascript
        # div with id messagearea can be used to check whether query is matched perfectly
        # if text of <div id="message"> is empty, select the first returned entries
        ########################################################################################
        
        
        url = 'https://www.ncbi.nlm.nih.gov/pccompound/?term='+ urllib.parse.quote(query)
        res = requests.get(url)
        soup = BeautifulSoup(res.text, 'lxml')
        
        if soup.select_one('div#messagearea').text.strip() == '':
            # match perfectly
            cid = soup.select_one('div.aux > div > dl > dd').text
        else:
            return [None,None,None,None]

        # get info from PUGVIEW
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/'+str(cid)+'/JSON/'
        res = requests.get(url).json()['Record']['Section']

        CASlist = GetPropfromPUGVIEWJSON(cid, res, 'CAS')
       
        # check CAS
        for cas in CASlist:
            if query==cas:
                return [query,cid,GetPropfromPUGVIEWJSON(cid, res, 'smiles'),GetPropfromPUGVIEWJSON(cid, res, 'mw')]

        # check synonyms
        synlist = GetPropfromPUGVIEWJSON(cid, res, 'synonyms')
        # synlist may be None
        
        if synlist != None:
            for syn in synlist:
                if q_type == 'CAS':
                    if query in syn:
                        return [query,cid,GetPropfromPUGVIEWJSON(cid, res, 'smiles'),GetPropfromPUGVIEWJSON(cid, res, 'mw')]
                elif q_type == 'Name':
                    if query.lower() == syn.lower():
                        return [CASlist,cid,GetPropfromPUGVIEWJSON(cid, res, 'smiles'),GetPropfromPUGVIEWJSON(cid, res, 'mw')]
        return [None,None,None,None]

            
stime = time.time()
query = '7732-18-5'
r = GetMolIdentifier(query,'pubchem')
print(r, time.time()-stime)



