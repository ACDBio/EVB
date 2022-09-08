import pandas as pd
import numpy as np
import pickle
from tqdm.autonotebook import tqdm
import sys
'''The script remaps the dataframe to HUGO and saves is as mapped_df.csv (HUGO_mappings.pkl should be in the same folder).'''
def remap_dataframe(dataframe_filepath, dataframe_genecol, leave_only_genenames_and_scorecol=False, scorecol=None, leave_only_genenames=False, index_column='na', drop_old_genecol=True, drop_unmapped=True, aggregate_data_bygene=True, aggregation_function='max', save=False, savepath='./mapped_df.csv', mappingfilepath='./HUGO_mappings.pkl', toreturn=False, sep=','):
            print(f'Unifying gene nomenclature for {dataframe_filepath}')
            if index_column==0:
                dataframe=pd.read_csv(dataframe_filepath, index_col=0, sep=sep)
            else:
                dataframe=pd.read_csv(dataframe_filepath, sep=sep, error_bad_lines=False)
            #dataframe=dataframe.reset_index()
            mappingfile=open(mappingfilepath, "rb")
            mappings=pickle.load(mappingfile)
            genes= dataframe[dataframe_genecol]
            mapped_col=[]
            map_values=list(mappings.values())
            map_values=[item for sublist in map_values for item in sublist]
            map_values=np.array(map_values)
            for g in tqdm(genes):
                if g in mappings.keys():
                    mapped_col.append(g)
                elif g in map_values:
                    target_gene=[key for key in mappings if g in mappings[key]]
                    mapped_col.append(target_gene[0])
                else:
                    mapped_col.append('No_data')
            dataframe['pipe_genesymbol']=mapped_col
            print(f"Number of unmapped values: {dataframe[dataframe['pipe_genesymbol']=='No_data'].shape[0]}")
            if drop_old_genecol==True:
                dataframe=dataframe.drop(dataframe_genecol, axis=1)
            if drop_unmapped==True:
                dataframe=dataframe[dataframe['pipe_genesymbol']!='No_data']

            if leave_only_genenames==True:
                dataframe=dataframe[['pipe_genesymbol']]
                if save==True:
                    dataframe.to_csv(savepath, index=False)
                print(f'Gene symbol unification for {dataframe_filepath} - done.')
                if toreturn==False:
                    return
                else:
                    return dataframe
            elif leave_only_genenames_and_scorecol==True:
                dataframe=dataframe[['pipe_genesymbol', scorecol]]
                dataframe=dataframe.rename(columns={scorecol:'score'})
                if save==True:
                    dataframe.to_csv(savepath, header=False, index=False)
                print(f'Gene symbol unification for {dataframe_filepath} - done.')
                if toreturn==False:
                    return
                else:
                    return dataframe
            else:
                if aggregate_data_bygene==True:
                    dataframe=dataframe.groupby('pipe_genesymbol').agg(aggregation_function)
                if save==True:
                    dataframe.to_csv(savepath, sep=sep)
                print(f'Gene symbol unification for {dataframe_filepath} - done.')
                if toreturn==True:
                    return dataframe
                else:
                    return
remap_dataframe(dataframe_filepath='genedf.csv',dataframe_genecol='genes', save=True, aggregate_data_bygene=False, mappingfilepath='./Data/HUGO_mappings.pkl')
