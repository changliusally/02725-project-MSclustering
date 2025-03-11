
# 1. read spectrum information and do preprocessing in parallel, k core, using numpy, pandas
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
# spectra index with their top 50 peaks (m/z values)
spectra = {key: spectrum for key, spectrum in enumerate(mz_values)}
# this is a unique peak list of the whole MS Database (which generate bins considering the mz range)
#  {bin range: {spectra_index: indensities}}

# # spectrual binning
def binning(mz_values):
    """
    Efficiently bins m/z values into fixed-width bins using numpy.
    """
    mz_flat = np.concatenate(mz_values)  # Flatten the list of lists
    min_val, max_val = mz_flat.min(), mz_flat.max()  # Compute global min and max
    bin_edges = np.arange(min_val - 0.5, max_val + 1.5, 1)
    return bin_edges

bins = binning(mz_values)

def create_sparse_index(mz_values, intensities, bin_edges):
    """
    Creates a sparse index dictionary for MS data.
    """
    unique_peaks = {tuple([bin_edges[i], bin_edges[i+1]]): {} for i in range(len(bin_edges)-1)}

    # Process each spectrum
    for i, (mzs, ints) in enumerate(zip(mz_values, intensities)):
        # Normalize the spectrum
        norm_factor = max(ints)
        normalized_intensities = np.array(ints) / norm_factor * 1000

        # Digitize the m/z values to find corresponding bins
        bin_indices = np.digitize(mzs, bin_edges) - 1

        for mz, intensity, bin_idx in zip(mzs, normalized_intensities, bin_indices):
            if bin_idx < len(bin_edges) - 1:  # Ensure the index is within the range
                bin_key = (bin_edges[bin_idx], bin_edges[bin_idx+1])
                if i not in unique_peaks[bin_key]:
                    unique_peaks[bin_key][i] = intensity
                else:
                    # Store the maximum intensity for this bin
                    unique_peaks[bin_key][i] = max(unique_peaks[bin_key][i], intensity)

    keys_to_remove = [key for key, value in unique_peaks.items() if not value]
    for key in keys_to_remove:
        del unique_peaks[key]
    return unique_peaks

unique_peaks = create_sparse_index(mz_values,intensities,bins)

# print(unique_peaks)
end_time1 = time.time()
times.append(end_time1-start_time)

    
# 2. transforming the m/z array and intensity array into sparse vectors.
# extracting the top k peaks from both spectrum and store in a set form, assign their intensities to corresponding peak position.
bins_index = unique_peaks.keys()
sparse_vector = {}

for spectrum_index in range(len(spectra.keys())):
    # Create a sparse vector with default intensity of 0 for bins where the spectrum does not appear
    spectrum_vector = []
    for bin_key in bins_index:
        # Check if the current spectrum has an entry in the current bin
        if spectrum_index in unique_peaks[bin_key]:
            spectrum_vector.append(unique_peaks[bin_key][spectrum_index])
        else:
            spectrum_vector.append(0)
    sparse_vector[spectrum_index] = np.array(spectrum_vector)

# print(sparse_vector[0])


# 3. Implement hierarchical clustering on the sparse vectors.
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
import matplotlib.pyplot as plt
data = list(sparse_vector.values())

  
# Compute the pairwise distances only once
condensed_distance = pdist(data, 'cosine')  # this is the condensed distance matrix

end_time2 = time.time()
times.append(end_time2-start_time)

k = 12000
Z = linkage(condensed_distance, 'average')
clusters = fcluster(Z, k, criterion='maxclust')

# Alternatively, cut the dendrogram at a certain distance
distance_threshold = 0.7
clusters_by_dist = fcluster(Z, distance_threshold, criterion='distance')

# Assume you have some identifiers for each sample if not just use a range index
sample_ids = range(len(clusters))

# Create a DataFrame to hold your data
df = pd.DataFrame({
    'Sample ID': sample_ids,
    'Cluster Maxclust': clusters,
    'Cluster Distance': clusters_by_dist
})

# Write the DataFrame to a CSV file
df.to_csv('cluster_output_sparse.csv', index=False)


# Plotting dendrogram
# plt.figure(figsize=(10, 8))
# dendrogram(Z)
# plt.title('Hierarchical Clustering Dendrogram with sparse index(Cosine Distance)')
# plt.xlabel('Sample index')
# plt.ylabel('Distance')
# plt.savefig('dendrogram_sparse_index_100.png')  # Save the plot before showing it
# plt.show()
# plt.close()

end_time3 = time.time()
times.append(end_time3-start_time)

print(times)


from sklearn.metrics import silhouette_score

# Calculate the silhouette score
silhouette_avg = silhouette_score(condensed_distance, clusters)
print('The average silhouette score for hierarchical clustering is:', silhouette_avg)