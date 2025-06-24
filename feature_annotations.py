#!/usr/bin/env python3

import pandas as pd
import requests
import gget

#### With this file i can annotate lists 


def reg_to_list(regions, length):
    """
    Turns a list or tuples for regions into a full list for each AA
    All regions will be annotated with "1", the rest with 0
    Needed for function annot_IDR()
    """
    seq = [0]*length
    # print(regions)
    # print(length)

    if regions:
        for el in regions:
            for i in range(el[0]-1,el[1]):
                try:
                    seq[i] = 1
                except IndexError:
                    seq = [-1]*length
    return seq


def annot_PTMs(protein_names: list, source: list = ["uniprot"], return_failed_attempts: bool = False, show_progress: bool = False, verbose: bool = False) -> pd.DataFrame:
    """
    This function will find all post-translational modifications (PTMs) for the provided list of proteins from the chosen source and output them to you as a dataframe.

    Please provide the list of proteins as UniProt Protein IDs.

    Accepted sources are currently: "uniprot"
    """
    
    # Check that the provided source is valid
    if source != ["uniprot"]:
        print("The source (" + source + ") has not been implemented yet.")
        return pd.DataFrame()
    
    # Initialize empty lists for successful and failed requests, and a cache for results and a timeout
    result_list, failed_list, failed_list_e = [], [], []
    cache = {}
    timeout = 10
    
    # Loop through the provided list of proteins
    for index, prot in enumerate(protein_names):
        if show_progress:
            if index % 100 == 0:
                print(index, "/", len(protein_names),"done!")
        # print(prot)
        # If the protein's PTMs have already been retrieved, use cached results
        if prot in cache:
            payload = cache[prot]
        # Otherwise, query the UniProt API for the protein's PTMs
        else:
            url = f"https://www.ebi.ac.uk/proteins/api/features/{prot}?categories=PTM"
            try:
                response = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
                response.raise_for_status()
                payload = response.json()
                cache[prot] = payload
            # Handle any exceptions that may occur during the API request
            except requests.exceptions.RequestException as e:
                if verbose:
                    print(f"Error accessing UniProt API for protein {prot}: {e}")
                failed_list.append(prot)
                failed_list_e.append(e)
                continue
        
        # Parse the PTM information for the protein and store it in a list
        # print(payload)
        for ptm in payload['features']:
            # print(ptm)
            temp_dict={}
            temp_dict['protein_name'] = prot
            temp_dict['length'] = len(payload['sequence'])
            if 'evidences' not in ptm:
                continue
            try:
                temp_dict['num_evidences'] = len(ptm['evidences'])
                temp_dict2 = dict((k, ptm[k]) for k in ('type', 'description')) # 'begin', 'end', 'molecule' 
                temp_dict2["pos"] = ptm["begin"]
                temp_dict2["AA"] = payload['sequence'][int(ptm["begin"])-1]
                temp_dict2["ptm"] = ptm["description"].split(";")[0]
                temp_list= []
                if len(ptm["description"].split(";")) > 1:
                    if "by" in ptm["description"].split(";")[1]:
                        for by in ptm["description"].split(";")[1].replace(" and", ",").replace("by ", "").split(","):
                            temp_list.append(by.strip())
                temp_dict2["by"] = temp_list
                result_list.append({**temp_dict, **temp_dict2})
            except Exception as e:
                if verbose:
                    print(prot)
                    print(ptm)
                    print(ptm["begin"])
                    print(e)
                    print('__________________')
                continue
    
    # Convert the list of PTMs into a pandas dataframe
    results_df = pd.DataFrame(result_list)
    
    # If specified, return a tuple containing the list of failed requests and the PTM dataframe
    failed_df = pd.DataFrame({'failed_protein': failed_list, 'error': failed_list_e})

    if return_failed_attempts:
        return failed_df, results_df
    else:
        return results_df

def annot_IDR(protein_names: list, return_failed_attempts: bool = False, show_progress: bool = False, verbose: bool = False, cache_mode: bool = False) -> pd.DataFrame:
    """
    This will give find the IDRs for the provided list of proteins from the chosen source  and output them to you as a dataframe.\n
    Please provide the list of proteins as UniProt Protein IDs./n
    Accepted sources are currently: "iupred", "mobidb"
    """
    # Initialize empty lists for successful and failed requests, and a cache for results and a timeout
    result_list, failed_list, failed_list_e = [], [], []
    cache = [{},{}]
    timeout = 10
    for index,prot in enumerate(protein_names):
        if show_progress:
            if index % 100 == 0:
                print(index, "/", len(protein_names),"done!")
        if prot in cache:
            payload = cache[1][prot]
        else:
            url= f"https://mobidb.org/api/download_page?format=json&acc={prot}"
            try:
                response = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
                response.raise_for_status()
                payload = response.json()['data'][0]
                if cache_mode:
                    cache[1][prot] = payload
            # Handle any exceptions that may occur during the API request
            except (requests.exceptions.RequestException, ValueError, IndexError, KeyError) as e:
                # print(type(e))
                if verbose & (e == requests.exceptions.RequestException):
                    print(f"Error accessing UniProt API for protein {prot}: {e}")
                    print(index)
                elif verbose & (e == IndexError):
                    print("Index Error!!! Empty API call, missing data in mobidb lite")
                failed_list.append(prot)
                failed_list_e.append(e)
                continue
        # print(payload)
        # print(payload['data'])
        # print(payload[1])



        temp_dict={}
        temp_dict['protein_name'] = prot
        temp_dict['length'] = len(payload['sequence'])
        # if old_seq != payload['sequence']:
        #     print("there is no common sequence!!!")
        #     temp_dict['common_seq'] = False
        # else:
        #     temp_dict['common_seq'] = True
        # print(payload.keys())
        # print(payload)


        # print(temp_dict['length'])
        # print(len(payload['sequence']))
        temp_dict2 = {k: reg_to_list(payload.get(k, {}).get('regions'),len(payload['sequence']))  for k in ('curated-disorder-merge', 'prediction-disorder-mobidb_lite', 'prediction-disorder-alphafold')}
        # temp_dict2 = {k: reg_to_list(payload.get(k, {}).get('regions'),len(payload['sequence']))  for k in ('curated-disorder-merge', 'prediction-disorder-priority', 'prediction-disorder-alphafold')}
        # temp_dict2 = {k: reg_to_list(payload.get(k, {}).get('regions'),len(payload['sequence']))  for k in ('curated-disorder-merge', 'prediction-disorder-priority', 'prediction-disorder-alphafold')}
        # temp_dict2 = {k: reg_to_list(payload.get(k, {}).get('regions'),len(payload['sequence']))  for k in ('curated-disorder-merge', 'prediction-disorder-priority', 'prediction-disorder-alphafold')}
        temp_dict2 = {**temp_dict, **temp_dict2}
        temp_dict2['prediction-disorder-alphafold-score'] = payload.get('prediction-disorder-alphafold',{}).get('scores')
        result_list.append(temp_dict2)
    results_df = pd.DataFrame(result_list)
    failed_df = pd.DataFrame({'failed_protein': failed_list, 'error': failed_list_e})
    # If specified, return a tuple containing the list of failed requests and the PTM dataframe
    if return_failed_attempts:
        return failed_df, results_df
    else:
        return results_df

def annot_domains(protein_names: list, source: list = ["interpro"], return_failed_attempts: bool = False, show_progress: bool = False) -> pd.DataFrame:
    """
    This will find all domains for the provided list of proteins from the lsit of sources and output them to you as a dataframe.\n
    Please provide the list of proteins as UniProt Protein IDs./n
    Accepted sources are currently: "interpro"
    """
    # Initialize empty lists for successful and failed requests, and a cache for results and a timeout
    result_list, failed_list, failed_list_e = [], [], []   # lists to store successful and failed requests
    cache = {}   # dictionary to cache results
    timeout = 10   # timeout for requests in seconds
    
    for index, prot in enumerate(protein_names):
        # print(prot)
        if show_progress:
            if index % 100 == 0:
                print(index, "/", len(protein_names),"done!")
        if prot in cache:
            payload = cache[prot]
        else:
            url = f"https://www.ebi.ac.uk:443/interpro/api/entry/interpro/protein/uniprot/{prot}"
            try:
                # Make the API request to get the data for the protein
                response = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
                response.raise_for_status()  # raise an exception if response status is not 2xx
                payload = response.json()   # parse the JSON response into a Python dictionary
                cache[prot] = payload   # add the result to the cache
            # Handle any exceptions that may occur during the API request
            except requests.exceptions.RequestException as e:
                print(f"Error accessing UniProt API for protein {prot}: {e}")
                failed_list.append(prot)   # add the protein to the list of failed requests
                failed_list_e.append(e)
                continue

        # Extract domain information from the payload and add it to the result list
        for el in payload['results']:
            # print(el)
            temp_dict={}       
            if el['metadata']['type'] == 'domain':
                temp_dict['protein_name'] = prot   # add the protein name to the dictionary
                temp_dict['databases'] = list(el['metadata']['member_databases'].keys())   # add list of databases to the dictionary
                temp_dict2 = dict((k, el['metadata'][k]) for k in ('accession', 'name'))   # add the accession and name to the dictionary
                temp_dict2['GO_identifiers'] = []   # create an empty list for GO identifiers
                temp_dict2['GO_names'] = []   # create an empty list for GO names
                if el['metadata']['go_terms']:   # if there are GO terms, extract the identifier and name and add them to the corresponding lists
                    for id in el['metadata']['go_terms']:
                        temp_dict2['GO_identifiers'].append(id['identifier'])
                        temp_dict2['GO_names'].append(id['name'])
                if el['proteins'][0]['entry_protein_locations'] is None:
                    continue
                temp_dict3 = el['proteins'][0]['entry_protein_locations'][0]['fragments'][0]   # extract the domain location information
                result_list.append({**temp_dict, **temp_dict2, **temp_dict3})   # combine the 3 dictionaries and add to the result list

    results_df = pd.DataFrame(result_list)
    failed_df = pd.DataFrame({'failed_protein': failed_list, 'error': failed_list_e})
    # If specified, return a tuple containing the list of failed requests and the domain dataframe
    if return_failed_attempts:
        return failed_df, results_df
    else:
        return results_df

def annot_GO(protein_names: list, return_failed_attempts: bool = False, show_progress: bool = False, verbose: bool = False):
    # protein_names = pd.read_csv(r'/mnt/d/phd/scripts/0_motif_exploration/output/list_of_all_proteins.csv', header=None, names=["UniqueID"] )['UniqueID'].tolist()

    list_of_GO_terms_withRNAbinding = open(r'/mnt/d/phd/scripts/raw_data/GO_terms/GO_terms_RNAbinding.txt', "r").read().split("\n")[:-1]
    list_of_GO_terms_withNAbinding = open(r'/mnt/d/phd/scripts/raw_data/GO_terms/GO_terms_NAbinding.txt', "r").read().split("\n")[:-1]
    list_of_GO_terms_withDNAbinding = open(r'/mnt/d/phd/scripts/raw_data/GO_terms/GO_terms_DNAbinding.txt', "r").read().split("\n")[:-1]


    timeout = 2   # timeout for requests in seconds

    result_list, failed_list, failed_list_e = [], [], []    # lists to store successful and failed requests
    temp_dict = {}
    for index,prot in enumerate(protein_names):
        if show_progress:
            if index % 100 == 0:
                print(index, "/", len(protein_names),"done!")
        go_terms, go_aspects= [], []
        # print(prot)
        # print(i)
        # for prot in ['A0A009G0S2', 'P01138', 'P01594', 'P01833', 'P08134', 'P08754', 'P0CG29']:
        url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search?geneProductId={prot}"
        try:
            # Make the API request to get the data for the protein
            response = requests.get(url, headers={"Accept": "application/json"}, timeout=timeout)
            response.raise_for_status()  # raise an exception if response status is not 2xx
            payload = response.json()   # parse the JSON response into a Python dictionary
            # cache[prot] = payload   # add the result to the cache
        # Handle any exceptions that may occur during the API request
        except requests.exceptions.RequestException as e:
            if verbose:
                print(f"Error accessing UniProt API for protein {prot}: {e}")
            failed_list.append(prot)   # add the protein to the list of failed requests
            failed_list_e.append(e)
            continue
            # print(payload)
        # print(payload.keys())
        # print(len(payload["results"]))
        counts_RNA, counts_DNA, counts_NA = 0,0,0
        cache_small = []
        for el in payload["results"]:
            if el['goId'] in cache_small:
                continue
            cache_small.append(el['goId'])
            # print(el)
            # print(el['goId'])
            go_terms.append(el['goId'])
            go_aspects.append(el['goAspect'])
            if el['goId'] in list_of_GO_terms_withRNAbinding:
                counts_RNA += 1
            if el['goId'] in list_of_GO_terms_withDNAbinding:
                counts_DNA += 1
            if el['goId'] in list_of_GO_terms_withNAbinding:
                counts_NA += 1
                # print(prot)
            # print(el['goAspect'])
        # print(counts)
        # print("_")
        temp_dict = {}
        temp_dict["UniqueID"] = prot
        temp_dict["go_terms"] = go_terms
        temp_dict["go_aspects"] = go_aspects
        temp_dict["invs_RNAbind"] = counts_RNA
        temp_dict["invs_DNAbind"] = counts_DNA
        temp_dict["invs_NAbind"] = counts_NA
        # go_terms_all.append(go_terms)
        # go_aspects_all.append(go_aspects)
        # invs_RNAbind_all.append(counts)
        # print(temp_dict)
        # print(result_list)
        result_list.append(temp_dict)
        # print(result_list)
        # print("____")
    results_df = pd.DataFrame(result_list)
    failed_df = pd.DataFrame({'failed_protein': failed_list, 'error': failed_list_e})
        # If specified, return a tuple containing the list of failed requests and the domain dataframe
    if return_failed_attempts:
        return failed_df, results_df
    else:
        return results_df

def annot_ELM(protein_names: list, source: list = ["ELM"], return_failed_attempts: bool = False, show_progress: bool = False, verbose: bool = False) -> pd.DataFrame:
    """
    This function creates a large dataframe that contains motifs in a sequence, found by the ELM
    THis data is only filtered by the true_positive filter in the InstanceLogicn columns
    here we use the gget elm package, builds on ELM, but much quicker and efficient.
    """
    
    # Check that the provided source is valid
    if source != ["ELM"]:
        print("The source (" + source + ") has not been implemented yet.")
        return pd.DataFrame()
    
    failed_list = []
    regex_dfs = []
    gget.setup("elm")
    for index, prot in enumerate(protein_names):
        ortholog_df, regex_df = gget.elm(prot, uniprot=True, verbose=False)
        if regex_df.empty:
            failed_list.append(prot)
        else:

            regex_df['protein_name'] = [prot]*len(regex_df)
            regex_dfs.append(regex_df)

    results_df = pd.concat(regex_dfs, ignore_index=True)
    results_df = results_df[results_df["InstanceLogic"] == "true positive"]
    
    # If specified, return a tuple containing the list of failed requests and the PTM dataframe
    failed_df = pd.DataFrame({'failed_protein': failed_list, 'error': 'na'})

    if return_failed_attempts:
        return failed_df, results_df
    else:
        return results_df


if __name__ == "__main__":
    # Code to be executed when the file is run directly goes here
    print(annot_ELM(["Q99102", "Q8NAF0"]))
    # print(annot_PTMs(['O94911', 'P01138', 'P01594', 'P01833', 'P08134', 'P08754', 'P0CG29']))
    # print(annot_domains(["A6NFQ2", "P04637", "Q9Y6I9"]))
