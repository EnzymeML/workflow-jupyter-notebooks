'''
Created on 11.05.2021

@author: HD
'''

import os
import requests
import json
import numpy as np

class JSONHelper:

    __apiURL = None
    __enzmlJSON = None

    def __init__(self, apiURL):
        '''
        Helper Class to work with JSON-formatted EnzymeML document

        Args:
            apiURL (string): url to API
        '''
        self.__apiURL = apiURL

    def readRequest(self, omexPath):
        '''
        Send read request with EnzymeML document to API, saves JSON response.

        Args:
            omexPath (string): path to omex file
            omexName (string): name of omex file
        Returns:
            int: status code of API response (200 if everthing work, 500 if error occurred)
        '''
        omexName = os.path.basename(omexPath)
        endpoint_read = f"{self.__apiURL}/read"
        payload={}
        files=[('omex',(omexName,open(omexPath,'rb'),'application/octet-stream'))]
        headers={}
        response = requests.request("GET", endpoint_read, headers=headers, data=payload, files=files)
        if response.status_code==200:
            self.__enzmlJSON = EnzymeMLDocJSON(response.content)
        return response.status_code

    def getEnzmlJSON(self):
        '''
        Returns object with EnzymeML information
        '''
        return self.__enzmlJSON

class EnzymeMLDocJSON(object):

    protDic = {}
    reactantDic = {}
    reactionDic = {}
    typeDic = {}
    
    def __init__(self, JSON):
        '''
        Class containing all information of JSON-formatted EnzymeML document

        Args:
            JSON (json): JSON-formatted EnzymeML document
        '''
        self.__enzmlJSON = json.loads(JSON)
        self._setProteinDic()
        self._setReactantDic()
        self._setReactionDic()
        self.__str__()

    def __str__(self):
        '''
        Magic function return pretty string describing the object
        '''
        fin_string = ["Title: "+self.getTitle()]
        fin_string.append("Proteins: ")
        for id, name in self.protDic.items():
            fin_string.append('\tID: %s \t Name: %s' % ( id, name))
        fin_string.append("Reactants: ")
        for id, name in self.reactantDic.items():
            fin_string.append('\tID: %s \t Name: %s' % ( id, name))
        fin_string.append("Reactions: ")
        for id, reac in self.reactionDic.items():
            fin_string.append('Name: %s \t ID: %s' % (reac['name'], id))
            fin_string.append(">>> Educts: ")
            for educt in reac['educts']:
                fin_string.append('\tID: %s \t Name: %s' % ( educt['species'], self.reactantDic[educt['species']]))
                self.typeDic[educt['species']] = 'educts'
            fin_string.append(">>> Products: ")
            for product in reac['products']:
                fin_string.append('\tID: %s \t Name: %s' % ( product['species'], self.reactantDic[product['species']]))
                self.typeDic[product['species']] = 'products'
            fin_string.append(">>> Modifiers: ")
            for mod in reac['modifiers']:
                try:
                    fin_string.append('\tID: %s \t Name: %s' % ( mod['species'], self.reactantDic[mod['species']]))
                    self.typeDic[mod['species']] = 'modifiers'
                except KeyError:
                    pass
                try:
                    fin_string.append('\tID: %s \t Name: %s' % ( mod['species'], self.protDic[mod['species']]))
                    self.typeDic[mod['species']] = 'modifiers'
                except KeyError:
                    pass
        return "\n".join(fin_string)

    def getTitle(self):
        '''
        Returns Title of EnzymeML document
        '''
        return self.__enzmlJSON['name']

    def _setProteinDic(self):
        for p in self.__enzmlJSON['protein']:
            self.protDic[p['id']] = p['name']

    def _setReactantDic(self):
        for r in self.__enzmlJSON['reactant']:
            self.reactantDic[r['id']] = r['name']

    def _setReactionDic(self):
        for r in self.__enzmlJSON['reaction']:
            self.reactionDic[r['id']] = r
    
    def getData(self, reaction, species):
        time = []
        data = []
        for reactant in self.reactionDic[reaction][self.typeDic[species]]:
            if reactant['species'] == species:
                for rep in reactant['replicates']:
                    time = rep['time']
                    data.append(rep['data'])
        return (np.array(time), np.array(data))