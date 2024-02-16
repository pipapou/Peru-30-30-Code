import subprocess
import xml.etree.ElementTree as ET
from matplotlib import pyplot as plt
from fonctions_spatial import read_geotiff, write_geotiff
import numpy as np
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import rasterio
import geopandas as gp
from rasterio.features import shapes
from rasterio.crs import CRS

def cp():

    # Parsing XML
    tree = ET.parse("results/parametres.xml")
    root = tree.getroot()

    # Extract DIM, X, Y and nb_protect from XML
    rat_protect = root.find("rat_protect").text
    nb_protect = root.find("nb_protect").text
    ratio = root.find("ratio").text
    X = root.find("X").text
    Y = root.find("Y").text

    # right type
    rat_protect = float(rat_protect)
    nb_protect = int(nb_protect)
    ratio = float(ratio)
    X = int(X)
    Y = int(Y)
    nb_protect_cp = str(round(nb_protect*(1-ratio)))

    jar_file = 'restopt-2.1.0.jar'

    # args
    habitat_file = 'results/Protect_data_10km_ilp_'+str(int(rat_protect*100))+'.tif'
    accessible_file = 'results/Masked_data_10km_ilp.tif'
    accessible_val = "-1"
    no_data_val = "-1"
    output_raster = 'results/Protect_data_10km_ilp_cp_'+str(int(rat_protect*100))+'.tif'
    output_csv = 'results/Protect_data_10km_ilp_cp_'+str(int(rat_protect*100))+'.csv'
    region = 'data/10km/10km_SERNANP_ecoregions.tif'
    minpu = nb_protect_cp
    maxpu = nb_protect_cp
    dom = 'NEIGHBORHOOD'
    opt = 'MAX_MESH_CWA_FRAC' # MIN_NP ; MAX_MESH
    timeout_seconds = "180s"
    buf = '1'
    nbker = '4'
#
    # commande
    command = ['java', '-Xss1g', '-jar', jar_file,  '-t', timeout_seconds, '-v',
               '-minpu', minpu, '-maxpu', maxpu, '-opt', opt, '-p', '-c', nbker, '-b', buf, '-r', region,
               habitat_file, accessible_file, accessible_val, no_data_val, output_raster, output_csv]

    command_str = ' '.join(command)
    print(command_str)

    # Execution subprocess.Popen for real time
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, universal_newlines=True)

    # real time execution
    while True:
        output_line = process.stdout.readline()
        if output_line == '' and process.poll() is not None:
            break
        if output_line:
            print(output_line.strip())

    process.wait()

    # output
    stdout, stderr = process.communicate()

    ### Graphic view

    result, model_tif_ds = read_geotiff("./results/Protect_data_10km_ilp_cp_"+str(int(rat_protect*100))+".tif")
    result = result.astype(float)

    result_bin, model_tif_ds = read_geotiff("results/Protect_bin_10km_ilp_"+str(int(rat_protect*100))+".tif")

    # Creation of numpy array

    for i in range(X):
        for j in range(Y):
            if result[i,j] < 0:
                result[i,j] = np.nan
            if result[i,j] == 1:
                result[i,j] = 0
            if (result[i,j] == 2) & (result_bin[i,j] == 0):
                result[i,j] = result[i,j]-1

                # graphics
    colors = ["#30678d", "#35b778", "#FFD700", "#B87333"]
    cmap = ListedColormap(colors)
    im = plt.imshow(result, cmap=cmap)
    values = [0, 1, 2, 3]
    colors = [ im.cmap(im.norm(value)) for value in values]
    labels = ["Unprotected", "Current PCA", "To Conserve (Step 1)", "To Conserve (Step 2)"]
    # create a patch for every color
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(values))]
    # put those patched as legend-handles into the legend
    plt.axis('off')
    plt.legend(handles=patches, bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.title("Output_ilp_cp", multialignment='left')
    plt.savefig('./results/Output_ilp_cp.png')
    plt.show()

    result_bin[result==3]=1

    write_geotiff("results/Protect_bin_10km_ilp_cp_"+str(int(rat_protect*100))+".tif", result_bin, model_tif_ds)

    mask = None
    with rasterio.Env():
        with rasterio.open('results/Protect_bin_10km_ilp_cp_'+str(int(rat_protect*100))+'.tif') as src:
            image = src.read(1)  # first band
            mask = image == 1
            results = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for i, (s, v)
                in enumerate(
                shapes(image, mask=mask, transform=src.transform)))
    geoms = list(results)
    gpd_polygonized_raster = gp.GeoDataFrame.from_features(geoms)
    target_crs = CRS.from_epsg(32718)

    # Define CRS for GeoDataFrame
    gpd_polygonized_raster.crs = target_crs
    gpd_polygonized_raster.to_file("data/shp/New_PA")
