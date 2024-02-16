import glob
import numpy as np
from pulp import *
from osgeo import gdal
from matplotlib import pyplot as plt
from fonctions_spatial import read_geotiff, write_geotiff
import xml.etree.ElementTree as ET
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import math

def ref_val():
    ############################ Parameters initialization #################################
    print("*Parameters initialization*")

    ### Import data
    # Parsing XML
    tree = ET.parse("results/parametres.xml")
    root = tree.getroot()

    # Extracting DIM, X, Y and nb_protect from XML
    DIM = root.find("DIM").text
    X = root.find("X").text
    Y = root.find("Y").text
    nb_protect = root.find("nb_protect").text
    ratio = root.find("ratio").text
    ecosrc = root.find("ecosrc").text
    pourc_eco = root.find("pourc_eco").text

    # python variables in the right type
    X = int(X)
    Y = int(Y)
    nb_protect = int(nb_protect)
    ratio = float(ratio)
    nb_protect_ilp = int(nb_protect*ratio)
    ecosrc = str(ecosrc)
    pourc_eco = float(pourc_eco)

    # import rasters
    Rasters = {}

    for rast in glob.glob("./data/" + DIM + "km/" + DIM + "km*"):
        temp = gdal.Open(rast)
        name = re.search(DIM + 'km_(.+?).tif', rast).group(1)
        Rasters["{0}".format(name)] = np.array(temp.ReadAsArray())

    # Creation of the mask of not relevant pixels
    Mask, Mask_ds = read_geotiff(
        "results/Masked_data_10km.tif")

    #setting weights
    vars_dict = {}
    vars_dict['CN_biodiversity'] = 1  # biodiversity
    vars_dict['SERNANP_biological_importance'] = 1 # biodiversity
    vars_dict['WCS_fonct_connect_A'] = 1  # biodiversity
    vars_dict['CN_carbon'] = 1  # carbon
    vars_dict['CN_water'] = 1 # water
    vars_dict['MINAM_water_retention_index'] = 1  # water
    vars_dict['PNUD_human_footprint'] = -1 # human pressure
    vars_dict['bosques_degradated_forest'] = -1  # human pressure

    for var in ['CN_biodiversity', 'CN_carbon', 'CN_water', 'PNUD_human_footprint', 'bosques_degradated_forest',
            'SERNANP_biological_importance', 'WCS_fonct_connect_A', 'MINAM_water_retention_index']:

            ############################### Model optimization #####################################################

        for value in ['opti', 'worst']:
            print(var, value)
            print("*Model optimization building*")

            ### Creating the problem for bests and worst values
            if value == 'opti':
                model = pulp.LpProblem('Conservation Planning Peru', LpMaximize)

            if value == 'worst':
                model = pulp.LpProblem('Conservation Planning Peru', LpMinimize)

            ### Definition of decision variables
            protected_pix = pulp.LpVariable.dicts("protected_pix", [(i,j) for i in range(X) for j in range(Y)], cat=LpBinary)

            ### Definition of the objective function
            model += pulp.lpSum([protected_pix[(i,j)] * Rasters[var][i,j] * vars_dict[var] for i in range(X) for j in range(Y)])

            ### Constraints
            model += pulp.lpSum([protected_pix[(i,j)] for i in range(X) for j in range(Y)]) == nb_protect_ilp  # number of pixels to protect + already protected = 30%

            model += pulp.lpSum([protected_pix[(i,j)] * Mask[i,j] for i in range(X) for j in range(Y)]) == 0 # masked pixels never choosen to be protected

            # at least 17% of pixels to be protected foreach ecoregion

            for val in np.unique(Rasters[ecosrc + '_ecoregions'][Rasters[ecosrc + '_ecoregions']>0]):
                ecoreg_nbmin_to_protect = math.ceil((pourc_eco - np.sum((Rasters[ecosrc + '_ecoregions'] == val) & (Rasters['Gouv_PA'] == 1)) /
                                            np.count_nonzero((Rasters[ecosrc + '_ecoregions'] == val)))*np.count_nonzero((Rasters[ecosrc + '_ecoregions'] == val)))
                mask_temp = np.where(Rasters[ecosrc + '_ecoregions'] == val, 1, 0)
                model += pulp.lpSum([protected_pix[(i, j)] * mask_temp[i, j] for i in range(X) for j in range(Y)]) >= ecoreg_nbmin_to_protect

            # structure PA constraints
            for i in range(1, X-1):
                for j in range(1, Y-1):

                    # list of neighbors of protected_pix[i,j]
                    neighbors = [protected_pix[(i-1,j)], protected_pix[(i+1,j)], protected_pix[(i,j-1)], protected_pix[(i,j+1)]]#

                    # list of neighbors of protected_pix[i,j] already PAs
                    neighbors_pa = [Rasters['Gouv_PA'][i-1,j], Rasters['Gouv_PA'][i+1,j], Rasters['Gouv_PA'][i,j-1], Rasters['Gouv_PA'][i,j+1]]

                    # no new protected pixel alone
                    model += protected_pix[(i,j)] <= pulp.lpSum(neighbors) + pulp.lpSum(neighbors_pa)

            ### Solving the problem
            print("*Solving*")

            #Choice of solver

            #solver = PULP_CBC_CMD()
            path_to_glpk = r'.\winglpk-4.65\glpk-4.65\w64\glpsol.exe'
            solver = GLPK_CMD(path=path_to_glpk, timeLimit = 1800)

            model.solve(solver)

            print("*Graphic view and export*")

            ### Graphic view

            # Creation of numpy array

            result = Rasters['DECVAR'].copy() - 1

            for i in range(X):
                for j in range(Y):
                    if result[i,j] < 0:
                        result[i, j] = np.nan
                    if Rasters['Gouv_PA'][i,j] == 1:
                        result[i,j] = 1
                    if protected_pix[(i,j)].varValue == 1:
                        result[i,j] = 2

            # graphics
            colors = ["#30678d", "#35b778", "#FFD700"]
            cmap = ListedColormap(colors)
            im = plt.imshow(result, cmap=cmap)
            values = [0, 1, 2]
            colors = [im.cmap(im.norm(value)) for value in values]
            labels = ["Peru", "Protected", "To Conserve (stage 1)"]
            # create a patch for every color
            patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(values))]
            # put those patched as legend-handles into the legend
            plt.axis('off')
            plt.legend(handles=patches, bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0., fontsize=10)
            plt.title("Output_ilp", multialignment='left')
            #plt.savefig('./results/Output_ilp.png')
            #plt.show()

            ### Raw output data

            # Creation of outputs

            protected_pix_np = np.array([protected_pix[(i,j)].varValue for i in range(X) for j in range(Y)]).reshape((X, Y))
            opt_val = protected_pix_np*Rasters[var]
            opt_val = np.sum(opt_val)

            # Créer un nouvel élément XML pour max_val et l'ajouter au fichier
            opt_val_element = ET.Element(value + '_' +var)
            opt_val_element.text = str(opt_val)
            root.append(opt_val_element)

            # Enregistrer le fichier XML mis à jour
            tree.write('./results/parametres.xml')