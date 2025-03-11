# calculate the cosine distance in pair
import numpy as np
import pandas as pd
from pyteomics import mzxml
from pyteomics.mzxml import MzXML
import time

# we record the time the program used to process the ms database, the time used to calculate distance matrix, the time it finished the hierarchical clustering.
times = []
start_time = time.time()


class MZXMLAnalysis:

    def __init__(self, mzxml_filepath: str):
        self.process_time: list = []
        self.mz: list = []
        self.intensity: list = []
        self.precursorMs : list = []
        self.data: MzXML = self.read_mzxml_file(mzxml_filepath)
        self.save_name = mzxml_filepath.split('.')[0] + '.xlsx'
        self.get_top_mz_and_intensity()

    @staticmethod
    def read_mzxml_file(file_path: str) -> MzXML:
        """
        get file data
        :param file_path: mzml file path
        :return: MzML object
        """
        if not file_path.endswith('mzXML'):
            raise TypeError(
                f'file must be .mzXML, your name is {file_path}'
            )
        # use_index必须指定为True 如果不指定m/z array 和 intensity array得到的都是None
        data = mzxml.read(file_path, use_index=True)
        return data
    def get_top_mz_and_intensity(self) -> None:
            """
            质谱仪数据提取 质荷比 丰度
            :return:
            """

            for spectrum in mzxml.read(self.data, use_index=True):
                if spectrum.get('msLevel') == 2:
                    if spectrum.get('precursorMz')[0]['precursorCharge'] == 2:
                        mz = np.array(spectrum.get('m/z array'))
                        if len(mz) > 5: 
                            intensity = np.array(spectrum.get('intensity array'))
                            # Filter the top 50 peaks
                            sorted_idx = np.argsort(intensity)[::-1][:50]
                            self.mz.append(mz[sorted_idx])
                            self.intensity.append(intensity[sorted_idx])
                            self.precursorMs.append(spectrum.get('precursorMz')[0]['precursorMz'])

# this file contain 17752 MS2 spectra, and we use 10894 spectra whose charge is 2.
file_path = './b1906_293T_proteinID_01A_QE3_122212.mzXML'
analysis = MZXMLAnalysis(file_path)
# they are list of top 50 peaks of all spectrum (list of float lists)
mz_values = analysis.mz
intensities = analysis.intensity
intensity = {key: [] for key, _ in enumerate(mz_values)}


# normalize the intensity
for i, ints in enumerate(intensities):
    # Normalize the spectrum
    norm_factor = max(ints)
    normalized_intensities = np.array(ints) / norm_factor * 1000
    intensity[i] = normalized_intensities

# spectra index with their top 50 peaks (m/z values)
spectra = {key: spectrum for key, spectrum in enumerate(mz_values)}

end_time1 = time.time()
times.append(end_time1-start_time)

# 2. calculating the distance
# given n points, there are n * (n - 1) / 2 pairwise distances

# index pair
from scipy.spatial.distance import pdist, squareform
# pair = (spec1, spec2)
# def custom_distance(pair):
#     all_mz = mz_values[pair[0]].copy()
#     vector1 = intensity[pair[0]].copy
#     vector2 = [0 for _ in vector1]
#     for mz in mz_values[pair[1]]:
#         for known_mz in mz_values[pair[0]]:
#             if abs(mz-known_mz) < 0.5:
#                 index1 = mz_values[pair[0]].index(known_mz)
#                 index2 = mz_values[pair[1]].index(mz)
#                 vector2[index1] = intensity[pair[1]][index2]
#                 continue
#             else:
#                 all_mz.append(mz)
#                 vector1.append(0)
#                 index2 = mz_values[pair[1]].index(mz)
#                 vector2.append(intensity[pair[1]][index2])
#     vectors = [vector1, vector2]
#     # Calculate pairwise distances using an appropriate metric, e.g., cosine
#     return pdist(vectors, 'cosine')[0]

def align_mz_values(mz1, int1, mz2, int2, tolerance=0.5):
    vector1 = []
    vector2 = []
    
    i, j = 0, 0
    while i < len(mz1) and j < len(mz2):
        if abs(mz1[i] - mz2[j]) <= tolerance:
            # If m/z values are within the tolerance, consider them as the same
            vector1.append(int1[i])
            vector2.append(int2[j])
            i += 1
            j += 1
        elif mz1[i] < mz2[j]:
            vector1.append(int1[i])
            vector2.append(0)  # No matching m/z in mz2
            i += 1
        else:
            vector1.append(0)  # No matching m/z in mz1
            vector2.append(int2[j])
            j += 1

    # Append remaining elements if any
    while i < len(mz1):
        vector1.append(int1[i])
        vector2.append(0)
        i += 1
    while j < len(mz2):
        vector1.append(0)
        vector2.append(int2[j])
        j += 1

    return vector1, vector2

def custom_distance(spec1, spec2, mz_values, intensities):
    mz1 = mz_values[spec1]
    int1 = intensities[spec1]
    mz2 = mz_values[spec2]
    int2 = intensities[spec2]
    
    # Align the m/z values with tolerance
    vector1, vector2 = align_mz_values(mz1, int1, mz2, int2)
    
    # Calculate pairwise distances using an appropriate metric, e.g., cosine
    return pdist([vector1, vector2], 'cosine')[0]

import itertools

# Assuming `num_spectra` is the number of spectra you have
num_spectra = len(spectra)  # For example
spectra_pairs = list(itertools.combinations(range(num_spectra), 2))

distance_mat = []
# Example to calculate distances for all pairs
for spec1, spec2 in spectra_pairs:
    distance = custom_distance(spec1, spec2, mz_values, intensity)
    distance_mat.append(distance)
    

end_time2 = time.time()
times.append(end_time2-start_time)


# 3. hierarchical clustering
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
Z = linkage(distance_mat, 'average')

k = 10000
clusters = fcluster(Z, k, criterion='maxclust')

# Alternatively, cut the dendrogram at a certain distance
distance_threshold = 0.5
clusters_by_dist = fcluster(Z, distance_threshold, criterion='distance')

# Assume you have some identifiers for each sample if not just use a range index
sample_ids = range(len(clusters))

# Create a DataFrame to hold your data
df = pd.DataFrame({
    'Sample ID': sample_ids,
    'Cluster Maxclust': clusters,
    'Cluster Distance': clusters_by_dist,
    'Precursor Mass': analysis.precursorMs
})

# Write the DataFrame to a CSV file
df.to_csv('cluster_output_normal_1.csv', index=False)



# plt.figure(figsize=(10, 8))
# dendrogram(Z)
# plt.title('Hierarchical Clustering without sparse index')
# plt.xlabel('Sample index')
# plt.ylabel('Distance')
# plt.savefig('dendrogram_normal.png')  # Save the plot before showing it
# plt.show()
# plt.close()

end_time3 = time.time()
times.append(end_time3-start_time)

print(times)

from sklearn.metrics import silhouette_score

# Calculate the silhouette score
silhouette_avg = silhouette_score(distance_mat, clusters)
print('The average silhouette score for hierarchical clustering is:', silhouette_avg)