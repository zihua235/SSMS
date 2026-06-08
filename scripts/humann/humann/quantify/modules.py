import os
import shutil
import math
import re
import sys
sys.path.append('/Users/zihuahuang/Documents/humann')
import subprocess
import logging

import chi2cdf
from humann import utilities
from humann import config
from humann import store

# name global logging instance
logger=logging.getLogger(__name__)

def minpath_command(reactions_file,metacyc_datafile):
    """
    Return the minpath command and the name of the output file
    """
    
    # Create temp files for the results
    tmpfile=utilities.unnamed_temp_file()
    tmpfile2=utilities.unnamed_temp_file()
    tmpfile3=utilities.unnamed_temp_file()     
    tmpfile4=utilities.unnamed_temp_file()  

    minpath_script=os.path.join(os.path.dirname(os.path.abspath(__file__)),
        config.minpath_script)
    args=[minpath_script,"-any",reactions_file]
    args+=["-map",metacyc_datafile,"-report",tmpfile3]
    args+=["-details",tmpfile,"-mps",tmpfile2]
        
    command=[sys.executable,args,[reactions_file,metacyc_datafile],[],tmpfile4,None,True,None]
    
    return tmpfile, command


def xipe_command(infile):
    """
    Return the xipe command and the name of the output files
    """
    
    xipe_exe=os.path.join(os.path.dirname(os.path.abspath(__file__)),
    config.xipe_script)
            
    args=[xipe_exe,"--file1",infile,"--file2",config.xipe_percent]
            
    stdout_file=utilities.unnamed_temp_file()
    stderr_file=utilities.unnamed_temp_file()
    
    command=[sys.executable,args,[infile],[],stdout_file,None,True,stderr_file]
    
    return stdout_file, stderr_file, command


def harmonic_mean(values):
    """
    Return the harmonic mean for the values
    """
    
    # If there are no values or if one of the values is zero, then the harmonic mean is zero
    mean=0
    if values and min(values) > 0:
        reciprocal_sum=sum((1.0/v) for v in values)
        mean=len(values)/reciprocal_sum
    
    return mean


def gap_fill(key_reactions, reaction_scores):
    """
    If all but one of the key reactions have abundance scores, then fill gap
    Boost the lowest abundance score
    """
    
    reaction_scores_gap_filled=reaction_scores.copy()
    
    # do not apply gap fill, if set to off
    if config.gap_fill_toggle == "off":
        return reaction_scores_gap_filled
    
    # get the scores for all of the key reactions
    key_reactions_nonzero_scores=[]
    for reaction in key_reactions:
        score=reaction_scores.get(reaction,0)
        if score > 0:
            key_reactions_nonzero_scores.append(score)

    if len(key_reactions)-len(key_reactions_nonzero_scores) == 1:
        # fill single zero gap with lowest key reaction score
        min_score=min(key_reactions_nonzero_scores)
        for reaction in key_reactions:
            score=reaction_scores.get(reaction,0)
            if score == 0:
                reaction_scores_gap_filled[reaction]=min_score
    elif len(key_reactions)-len(key_reactions_nonzero_scores) == 0:
        # boost lowest abundance score
        sorted_key_reactions_nonzero_scores=sorted(key_reactions_nonzero_scores)
        for reaction in key_reactions:
            score=reaction_scores.get(reaction,0)
            if score == sorted_key_reactions_nonzero_scores[0]:
                try:
                    reaction_scores_gap_filled[reaction]=sorted_key_reactions_nonzero_scores[1]
                except IndexError:
                    pass

    return reaction_scores_gap_filled

def identify_reactions_and_pathways_for_EC(EC_abundance, pathways_database):
    """
    Identify the reactions and then pathways from the hits found
    """

    #map ec to reactions
    bug = "all"
            
    if config.minpath_toggle == "on":
        # Write a flat reactions to pathways file
        logger.debug("Write flat reactions to pathways file for Minpath")
        pathways_database_file=utilities.unnamed_temp_file()
        file_handle=open(pathways_database_file,"w")
        file_handle.write(pathways_database.get_database())
        file_handle.close()
    
    # Create a store for the pathways and reactions by bug
    pathways_and_reactions_store=store.PathwaysAndReactions()
    reactions={}
    reactions[bug]={}
    reactions_store=store.Reactions()

    # set unmapped as unaligned reads count
    reactions_store.unmapped=EC_abundance['UNMAPPED']
    minpath_results={}
    minpath_commands=[]         
    reactions_file_lines=[]

    if reactions_database:
        for reaction in sorted(reactions_database.reaction_list()):
            ec_list=reactions_database.find_ec(reaction)
            
            #print(EC_abundance.keys())
            if ec_list in EC_abundance.keys():
                abundance=EC_abundance.get(ec_list)
                reactions[bug][reaction]=abundance
                reactions_store.add(bug,reaction,abundance)
                reactions_file_lines.append(reaction+config.output_file_column_delimiter
                        +str(abundance)+"\n")
            else:
                pass
    else:
        pass

    # Run minpath if toggle on and also if there is more than one reaction   
    if config.minpath_toggle == "on" and len(reactions_file_lines)>3:   
    
        # Create a temp file for the reactions results
        reactions_file=utilities.unnamed_temp_file()
        file_handle=open(reactions_file,"w")
        file_handle.write("".join(reactions_file_lines))
        file_handle.close()
        
        # Run minpath to identify the pathways
        logger.info("Run MinPath on ")
                
        tmpfile, command=minpath_command(reactions_file, pathways_database_file)
        minpath_results=tmpfile
        minpath_commands.append(command)
            
    # add the unintegrated reaction abundance for this bug to the total
    try:
        reactions_store.unintegrated_total=reactions_store.unintegrated["all"]
        del reactions_store.unintegrated["all"]
    except KeyError:
        reactions_store.unintegrated_total=0

    # Run through the minpath commands if minpath is to be run
    if minpath_commands:
        utilities.command_threading(config.threads,minpath_commands)
    
    # Link the pathways to reactions
    pathways={}
    tmpfile=minpath_results
    # Process the minpath results
    if os.path.isfile(tmpfile):
        
        file_handle_read=open(tmpfile, "rt")
        line=file_handle_read.readline()
               
        while line:
            data=line.strip().split(config.minpath_pathway_delimiter)
            if re.search(config.minpath_pathway_identifier,line):
                current_pathway=data[config.minpath_pathway_index]
            else:
                current_reaction=data[config.minpath_reaction_index]
                # store the pathway and reaction
                pathways[current_reaction]=pathways.get(
                    current_reaction,[]) + [current_pathway]      
            line=file_handle_read.readline()          
        file_handle_read.close()
      
        
    for current_reaction in reactions.get(bug,{}):
            # Find the pathways associated with reaction
            for current_pathway in pathways.get(current_reaction,[""]):
                # Only store data for items with pathway names
                if current_pathway:
                    #print("reactions[bug][current_reaction]",reactions[bug][current_reaction])
                    
                    pathways_and_reactions_store.add(bug,current_reaction, current_pathway, 
                        reactions[bug][current_reaction])     

    return pathways_and_reactions_store



#new function
def identify_reactions_and_pathways_for_EC(EC_abundance, pathways_database):
    """
    Identify the reactions and then pathways from the hits found
    """

    #map ec to reactions
    bug = "all"
            
    if config.minpath_toggle == "on":
        # Write a flat reactions to pathways file
        logger.debug("Write flat reactions to pathways file for Minpath")
        pathways_database_file=utilities.unnamed_temp_file()
        file_handle=open(pathways_database_file,"w")
        file_handle.write(pathways_database.get_database())
        file_handle.close()
    
    # Create a store for the pathways and reactions by bug
    pathways_and_reactions_store=store.PathwaysAndReactions()
    reactions={}
    reactions[bug]={}
    reactions_store=store.Reactions()

    # set unmapped as unaligned reads count
    reactions_store.unmapped=EC_abundance['UNMAPPED']
    minpath_results={}
    minpath_commands=[]         
    reactions_file_lines=[]

    if reactions_database:
        for reaction in sorted(reactions_database.reaction_list()):
            ec_list=reactions_database.find_ec(reaction)
            
            #print(EC_abundance.keys())
            if ec_list in EC_abundance.keys():
                abundance=EC_abundance.get(ec_list)
                reactions[bug][reaction]=abundance
                reactions_store.add(bug,reaction,abundance)
                reactions_file_lines.append(reaction+config.output_file_column_delimiter
                        +str(abundance)+"\n")
            else:
                pass
    else:
        pass

    # Run minpath if toggle on and also if there is more than one reaction   
    if config.minpath_toggle == "on" and len(reactions_file_lines)>3:   
    
        # Create a temp file for the reactions results
        reactions_file=utilities.unnamed_temp_file()
        file_handle=open(reactions_file,"w")
        file_handle.write("".join(reactions_file_lines))
        file_handle.close()
        
        # Run minpath to identify the pathways
        logger.info("Run MinPath on ")
                
        tmpfile, command=minpath_command(reactions_file, pathways_database_file)
        minpath_results=tmpfile
        minpath_commands.append(command)
            
    # add the unintegrated reaction abundance for this bug to the total
    try:
        reactions_store.unintegrated_total=reactions_store.unintegrated["all"]
        del reactions_store.unintegrated["all"]
    except KeyError:
        reactions_store.unintegrated_total=0

    # Run through the minpath commands if minpath is to be run
    if minpath_commands:
        utilities.command_threading(config.threads,minpath_commands)
    
    # Link the pathways to reactions
    pathways={}
    tmpfile=minpath_results
    # Process the minpath results
    if os.path.isfile(tmpfile):
        
        file_handle_read=open(tmpfile, "rt")
        line=file_handle_read.readline()
               
        while line:
            data=line.strip().split(config.minpath_pathway_delimiter)
            if re.search(config.minpath_pathway_identifier,line):
                current_pathway=data[config.minpath_pathway_index]
            else:
                current_reaction=data[config.minpath_reaction_index]
                # store the pathway and reaction
                pathways[current_reaction]=pathways.get(
                    current_reaction,[]) + [current_pathway]      
            line=file_handle_read.readline()          
        file_handle_read.close()
      
        
    for current_reaction in reactions.get(bug,{}):
            # Find the pathways associated with reaction
            for current_pathway in pathways.get(current_reaction,[""]):
                # Only store data for items with pathway names
                if current_pathway:
                    #print("reactions[bug][current_reaction]",reactions[bug][current_reaction])
                    
                    pathways_and_reactions_store.add(bug,current_reaction, current_pathway, 
                        reactions[bug][current_reaction])     

    return pathways_and_reactions_store


def compute_structured_pathway_abundance_or_coverage_EC(structure, key_reactions, reaction_scores, 
    coverage_computation, median_value):
    """
    Compute the abundance or coverage for a structured pathway
    """
    
    # Process through the structure to compute the abundance
    required_reaction_abundances=[]
    optional_reaction_abundances=[]
    # Select the join instead of removing from the list to not alter the list for
    # the calling function
    join=structure[0]
    for item in structure[1:]:
        if isinstance(item, list):
            required_reaction_abundances.append(compute_structured_pathway_abundance_or_coverage_EC(item, 
                key_reactions, reaction_scores, coverage_computation, median_value))
        else:
            score=reaction_scores.get(item,0)
                
            # Update the score for the reaction if this is a coverage computation
            if coverage_computation:
                score=chi2cdf.chi2cdf(score,median_value)

            # Check if this is an optional reaction
            if item in key_reactions:
                required_reaction_abundances.append(score)
            else:
                optional_reaction_abundances.append(score)
    
    # If this is an OR join then use the max of all of the reaction abundances
    if join == config.pathway_OR:
        all_reaction_abundances=required_reaction_abundances + optional_reaction_abundances
        abundance=0
        if all_reaction_abundances:
            abundance = max(all_reaction_abundances)
        
    else:
        # If this is not an OR, then take the harmonic mean of the reactions
        abundance=harmonic_mean(required_reaction_abundances)
        # Add the optional reactions if they are present
        if optional_reaction_abundances:
            # Filter the optional abundances to only include those that are greater than the abundance
            # from the required reactions
            optional_reaction_abundances_filtered=[value for value in optional_reaction_abundances if value > abundance]
            abundance=harmonic_mean(required_reaction_abundances + optional_reaction_abundances_filtered)
            
        
    return abundance

def compute_pathways_abundance_EC(pathways_and_reactions_store, pathways_database):
    """
    Compute the abundance of pathways for each bug
    Also find the set of the reactions with abundance in all pathways present
    """
    bug="all"
    # Store the reactions which have abundance in the pathways with abundance
    reactions_in_pathways_present={}
    reactions_in_pathways_present[bug]=set() 
    # Process through each pathway for each bug to compute abundance
    pathways_abundance_store=store.Pathways()
  
    for pathway in pathways_and_reactions_store.pathway_list(bug):
        reaction_scores=pathways_and_reactions_store.reaction_scores(bug, pathway)
        #print("reaction_score", reaction_scores)
            # Check if the pathways database is structured
        if pathways_database.is_structured():
            structure=pathways_database.get_structure_for_pathway(pathway)
            key_reactions=pathways_database.get_key_reactions_for_pathway(pathway)
                # Apply gap fill
            reaction_scores_gap_filled=gap_fill(key_reactions, reaction_scores)
    
                # Compute the structured pathway abundance
            abundance=compute_structured_pathway_abundance_or_coverage_EC(structure,
                key_reactions,reaction_scores_gap_filled,False,0)
            #print("reaction_scores_gap_filled", reaction_scores_gap_filled)
            
        else:
                # Initialize any reactions in the pathway not found to 0
            for reaction in pathways_database.find_reactions(pathway):
                reaction_scores.setdefault(reaction, 0)
                    
                # Sort the scores for all of the reactions in the pathway from low to high
            sorted_reaction_scores=sorted(reaction_scores.values())
                    
                # Select the second half of the list of reaction scores
            abundance_set=sorted_reaction_scores[int(len(sorted_reaction_scores)/ 2):]
                
                # Compute abundance
            abundance=sum(abundance_set)/len(abundance_set)
                
            # If this pathway is present, store those reactions with abundance
        if abundance > 0:
            for reaction,score in reaction_scores.items():
                if score > 0:
                    reactions_in_pathways_present[bug].add(reaction)
            
            # Store the abundance

        pathways_abundance_store.add('all',pathway, abundance)
    
    return pathways_abundance_store, reactions_in_pathways_present



utility_mapping_database = os.path.abspath("/Users/zihuahuang/Documents/humann") 
metacyc_gene_to_reactions=os.path.abspath(os.path.join(utility_mapping_database,"metacyc_reactions_level4ec_only.uniref.bz2"))
metacyc_reactions_to_pathways=os.path.abspath(os.path.join(utility_mapping_database,"metacyc_pathways_structured_filtered_v24_subreactions"))
pathways_database_part1=metacyc_gene_to_reactions
pathways_database_part2=metacyc_reactions_to_pathways
reactions_database=store.ReactionsDatabase(config.pathways_database_part1)
pathways_database=store.PathwaysDatabase(config.pathways_database_part2, reactions_database)

##################
import pandas as pd
#df_ec_abundance = pd.read_csv("/Users/zihuahuang/Documents/SSMS/data/imputed_ec/spieceasi/predicted_level4ec_spiec_Italy1asref_Italy2.csv", sep="\t", index_col=0)  # 行是 EC，列是样本名
df_ec_abundance = pd.read_csv("/Users/zihuahuang/Documents/SSMS/data/level4ec_and_pathway/ssms_level4ec_relab_unstratified/Italy2_500K_level4ec_relab_unstratified.tsv", sep="\t", index_col=0)  # 行是 EC，列是样本名
pathway_matrix = {}  # 最后将生成 pathway -> {sample: abundance}

for sample in df_ec_abundance.columns:
    print("sample", sample)
    ec_series = df_ec_abundance[sample]
    ec_dict = ec_series.to_dict()

    # 运行你已有的 pipeline
    pathways_and_reactions_store = identify_reactions_and_pathways_for_EC(ec_dict, pathways_database)
    pathways_abundance_store, _ = compute_pathways_abundance_EC(pathways_and_reactions_store, pathways_database)

    #print(pathways_abundance_store.get_pathways_list())
    # 提取 pathway abundance
    for pathway in pathways_abundance_store.get_pathways_list():
        abundance = pathways_abundance_store.get_score(pathway)
        #print("abundance", abundance)
        if pathway not in pathway_matrix:
            pathway_matrix[pathway] = {}
        pathway_matrix[pathway][sample] = abundance

# 转成 DataFrame：行为 pathway，列为样本
#print(pathway_matrix)
df_pathway_abundance = pd.DataFrame.from_dict(pathway_matrix, orient="index")
df_pathway_abundance = df_pathway_abundance.fillna(0)  # 缺失值补 0

# 输出查看
print(df_pathway_abundance.head())
df_pathway_abundance.to_csv("/Users/zihuahuang/Documents/SSMS/data/imputed_pw/predicted_Italy2_500k_level4ec_relab_unstratified.csv")
