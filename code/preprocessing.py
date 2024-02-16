import glob
import numpy as np
from pulp import *
from osgeo import gdal
import xml.etree.ElementTree as ET
from fonctions_spatial import read_geotiff, write_geotiff
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

def proprocess():
    ############################ Parameters initialization #################################

    print("*Parameters initialization*")

    ### Import data

    DIM = "10"  # pixels size

    Rasters = {}

    for rast in glob.glob("./data/" + DIM + "km/" + DIM + "km*"):
        temp = gdal.Open(rast)
        name = re.search(DIM + 'km_(.+?).tif', rast).group(1)
        Rasters["{0}".format(name)] = np.array(temp.ReadAsArray())

    # Maps dimensions

    X = Rasters['DECVAR'].shape[0]
    Y = Rasters['DECVAR'].shape[1]

    # Creation of the mask of not relevant pixels

    to_mask = [
        (Rasters['DECVAR'] == 0), # out of peru
        (Rasters['MODIS_IGBP_landuse'] == 12), # intensive agri
        (Rasters['MODIS_IGBP_landuse'] == 13), # built areas
        (Rasters['Gouv_PA'] != 0)
    ]

    Mask = np.any(to_mask, axis=0).astype(int)

    # percentage of peru to protect (actual + proposed)

    rat_protect = 0.3
    per_protect = rat_protect - Rasters['Gouv_PA'].sum() / Rasters['DECVAR'].sum()

    # Number of pixels to protect (to reach 30% of Peru)
    nb_protect = round(per_protect * Rasters['DECVAR'].sum())

    #ratio phase 1 / phase 2
    ratio = 0.8

    #parameters ecoregions
    ecosrc = "SERNANP"  # EPA SERNANP ; type of ecoregions map
    pourc_eco = 0.17 - (1/5) * per_protect # we remove the % of ecoregs equally allocated during phase 2

    # "control" information about the number of data
    print("pixels img: ", X * Y, "\npixels peru:", Rasters['DECVAR'].sum(), "\npixels masked: ", Mask.sum(),
          "\npixels selected: ",
          X * Y - Mask.sum(), "\npixels protected:", Rasters['Gouv_PA'].sum(), "\npixels to protect: ", nb_protect)

    # root element for XML storage
    root = ET.Element("Parametres")

    # storage parameters DIM, X, Y, nb_protect, etc.
    dim_element = ET.SubElement(root, "DIM")
    dim_element.text = str(DIM)

    x_element = ET.SubElement(root, "X")
    x_element.text = str(X)

    y_element = ET.SubElement(root, "Y")
    y_element.text = str(Y)

    rat_protect_element = ET.SubElement(root, "rat_protect")
    rat_protect_element.text = str(rat_protect)

    nb_protect_element = ET.SubElement(root, "nb_protect")
    nb_protect_element.text = str(nb_protect)

    ratio_element = ET.SubElement(root, "ratio")
    ratio_element.text = str(ratio)

    ecosrc_element = ET.SubElement(root, "ecosrc")
    ecosrc_element.text = str(ecosrc)

    pourc_eco_element = ET.SubElement(root, "pourc_eco")
    pourc_eco_element.text = str(pourc_eco)


    # ElementTree object and saving
    tree = ET.ElementTree(root)
    tree.write("results/parametres.xml")

    model_tif_arr, model_tif_ds = read_geotiff("./data/" + DIM + "km/" + DIM + "km_DECVAR.tif")
    write_geotiff("results/Masked_data_" + DIM + "km.tif", Mask, model_tif_ds)

    #graphic view existing PA

    result = Rasters['DECVAR'].copy() - 1

    for i in range(X):
        for j in range(Y):
            if result[i,j] < 0:
                result[i, j] = np.nan
            if Rasters['Gouv_PA'][i,j] == 1:
                result[i,j] = 1

    # graphics

    colors = ["#30678d", "#35b778"]
    cmap = ListedColormap(colors)
    im = plt.imshow(result, cmap=cmap)
    values = [0, 1]
    labels = ["Unprotected", "PCA"]
    # create a patch for every color
    patches = [mpatches.Patch(color=colors[i], label=labels[i] ) for i in range(len(values))]
    # put those patched as legend-handles into the legend
    plt.axis('off')
    plt.legend(handles=patches, bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.title("input_pa", multialignment='left')
    plt.savefig('./results/input_pa.png')
    plt.show()