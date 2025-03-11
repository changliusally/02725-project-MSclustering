import numpy as np
import pandas as pd
from pyteomics import mzxml
from pyteomics.mzxml import MzXML
from pyteomics import mzml
from pyteomics.mzml import MzML
from pyteomics import mgf
from pyteomics import pepxml

# class MZMLAnalysis:

#     def __init__(self, mzml_filepath: str):
#         self.process_time: list = []
#         self.mz: list = []
#         self.intensity: list = []
#         self.data: MzML = self.read_mzml_file(mzml_filepath)
#         self.save_name = mzml_filepath.split('.')[0] + '.xlsx'
#         self.get_top_mz_and_intensity()

#     @staticmethod
#     def read_mzml_file(file_path: str) -> MzML:
#         """
#         get file data
#         :param file_path: mzml file path
#         :return: MzML object
#         """
#         if not file_path.endswith('mzML'):
#             raise TypeError(
#                 f'file must be .mzML, your name is {file_path}'
#             )
#         # use_index必须指定为True 如果不指定m/z array 和 intensity array得到的都是None
#         data = mzml.read(file_path, use_index=True)
#         return data
#     def get_top_mz_and_intensity(self) -> None:
#             """
#             质谱仪数据提取 质荷比 丰度
#             :return:
#             """

#             for spectrum in mzml.read(self.data, use_index=True):
#                 if spectrum.get('msLevel') == '2':
#                     intensity = spectrum.get('intensity array')
#                     if len(intensity) > 50:
#                         threshold = intensity.sort(key=lambda x: x[1], reverse=True)[50]  # Sort by intensity
#                         index = [i for i, value in enumerate(intensity) if value >= threshold]
#                         self.mz.append(spectrum.get('m/z array')[index])
#                         self.intensity.append(spectrum.get('intensity array')[index])
#                     else:
#                          continue

# class MZXMLAnalysis:

#     def __init__(self, mzxml_filepath: str):
#         self.process_time: list = []
#         self.mz: list = []
#         self.intensity: list = []
#         self.data: MzML = self.read_mzxml_file(mzxml_filepath)
#         self.save_name = mzxml_filepath.split('.')[0] + '.xlsx'
#         self.get_all_mz_and_intensity()

#     @staticmethod
#     def read_mzxml_file(file_path: str) -> MzXML:
#         """
#         get file data
#         :param file_path: mzxml file path
#         :return: MzXML object
#         """
#         if not file_path.endswith('mzXML'):
#             raise TypeError(
#                 f'file must be .mzXML, your name is {file_path}'
#             )
#         # use_index必须指定为True 如果不指定m/z array 和 intensity array得到的都是None
#         data = mzxml.read(file_path, use_index=True)
#         return data
#     def get_all_mz_and_intensity(self) -> None:
#             """
#             质谱仪数据提取 采集时间 质荷比 丰度

#             :return:
#             """

#             for spectrum in mzxml.read(self.data, use_index=True):
#             # print(spectrum.get('centroid spectrum'))
#                 #self.process_time.append(spectrum.get('scanList').get('scan')[0].get('scan start time'))
#                 self.mz.append(spectrum.get('m/z array'))
#                 self.intensity.append(spectrum.get('intensity array'))

# class MGFAnalysis:

#     def __init__(self, mgf_filepath: str):
#         self.process_time: list = []
#         self.mz: list = []
#         self.intensity: list = []
#         self.data: MzML = self.read_mzml_file(mgf_filepath)
#         self.save_name = mgf_filepath.split('.')[0] + '.xlsx'
#         self.get_all_mz_and_intensity()

#     @staticmethod
#     def read_mgf_file(file_path: str) -> mgf:
#         """
#         get file data
#         :param file_path: mgf file path
#         :return: mgf object
#         """
#         if not file_path.endswith('mgf'):
#             raise TypeError(
#                 f'file must be .mgf, your name is {file_path}'
#             )
#         # use_index必须指定为True 如果不指定m/z array 和 intensity array得到的都是None
#         data = mgf.read(file_path, use_index=True)
#         return data
#     def get_all_mz_and_intensity(self) -> None:
#             """
#             质谱仪数据提取 采集时间 质荷比 丰度

#             :return:
#             """

#             for spectrum in mgf.read(self.data, use_index=True):
#             # print(spectrum.get('centroid spectrum'))
#                 #self.process_time.append(spectrum.get('scanList').get('scan')[0].get('scan start time'))
#                 self.mz.append(spectrum.get('m/z array'))
#                 self.intensity.append(spectrum.get('intensity array'))


# file_path = './b1945_293T_proteinID_09B_QE3_122212.mzXML'
# analysis = MZXMLAnalysis(file_path)
# mz_values = analysis.mz
# intensities = analysis.intensity


# print("This is the peaks' intensity", intensities.values[:10])

file_path = './b1906_293T_proteinID_01A_QE3_122212.pepXML'
psms = []

import pprint
# Just for inspection purposes
pp = pprint.PrettyPrinter(indent=4)

with pepxml.read(file_path) as reader:
    for i, spectrum in enumerate(reader):
        if i > 5:  # Limiting to first 5 spectra to avoid too much output
            break
        pp.pprint(spectrum)



with pepxml.read(file_path) as reader:
    for spectrum in reader:
        if 'search_hit' in spectrum:
            hits = spectrum['search_hit']
            for hit in hits:
                psm = {
                    'spectrum': spectrum['spectrum'],
                    'charge': spectrum['assumed_charge'],
                    'peptide': hit['peptide'],
                    'protein': hit['proteins'][0].get('protein'),
                    # 'peptide_prob': hit.get('peptideprophet_result', {}).get('probability', None),
                    'hyperscore': hit['search_score']['hyperscore']
                }
                psms.append(psm)

# Convert list to DataFrame
df_psm = pd.DataFrame(psms)
print(df_psm.head())
print(len(df_psm['protein'].unique()))

# Assuming 'peptide_prob' as the probability column
fdr_threshold = 0.01  # You need to adjust this based on your specific FDR calculation
filtered_df = df_psm[df_psm['peptide_prob'] >= fdr_threshold]
print(filtered_df.head())




# -----------------------------below are functions needed for the project
# import numpy as np
# import pandas as pd
# from pyteomics import mzxml
# from pyteomics.mzxml import MzXML
# class MZXMLAnalysis:

#     def __init__(self, mzxml_filepath: str):
#         self.process_time: list = []
#         self.mz: list = []
#         self.intensity: list = []
#         self.precursorMs : list = []
#         self.data: MzXML = self.read_mzxml_file(mzxml_filepath)
#         self.save_name = mzxml_filepath.split('.')[0] + '.xlsx'
#         self.get_top_mz_and_intensity()

#     @staticmethod
#     def read_mzxml_file(file_path: str) -> MzXML:
#         """
#         get file data
#         :param file_path: mzml file path
#         :return: MzML object
#         """
#         if not file_path.endswith('mzXML'):
#             raise TypeError(
#                 f'file must be .mzXML, your name is {file_path}'
#             )
#         # use_index必须指定为True 如果不指定m/z array 和 intensity array得到的都是None
#         data = mzxml.read(file_path, use_index=True)
#         return data
#     def get_top_mz_and_intensity(self) -> None:
#             """
#             质谱仪数据提取 质荷比 丰度
#             :return:
#             """

#             for spectrum in mzxml.read(self.data, use_index=True):
#                 if spectrum.get('msLevel') == 2:
#                     if spectrum.get('precursorMz')[0]['precursorCharge'] == 2:
#                         mz = np.array(spectrum.get('m/z array'))
#                         intensity = np.array(spectrum.get('intensity array'))
#                         # Filter the top 50 peaks
#                         sorted_idx = np.argsort(intensity)[::-1][:50]
#                         self.mz.append(mz[sorted_idx])
#                         self.intensity.append(intensity[sorted_idx])
#                         self.precursorMs.append(spectrum.get('precursorMz')[0]['precursorMz'])

# # this file contain 17752 MS2 spectra, all has polarity.
# file_path = './b1945_293T_proteinID_09B_QE3_122212.mzXML'
# analysis = MZXMLAnalysis(file_path)
# # they are list of top 50 peaks of all spectrum (list of float lists)
# mz_values = analysis.mz[:100]
# intensities = analysis.intensity[:100]
# print(len(mz_values[0]))
# print(len(analysis.precursorMs))

