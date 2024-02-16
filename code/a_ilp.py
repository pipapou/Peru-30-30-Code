import glob
import numpy as np
from pulp import *
from osgeo import gdal
from matplotlib import pyplot as plt
from scipy.ndimage import label
from fonctions_spatial import read_geotiff, write_geotiff
import xml.etree.ElementTree as ET
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import math

def ilp():
    ############################ Parameters initialization #################################

    print("*Parameters initialization*")

    ### Import data
    # Parsez le fichier XML
    tree = ET.parse("results/parametres.xml")
    root = tree.getroot()

    # Extract DIM, X, Y et nb_protect from XML
    DIM = root.find("DIM").text
    X = root.find("X").text
    Y = root.find("Y").text
    rat_protect = root.find("rat_protect").text
    nb_protect = root.find("nb_protect").text
    ratio = root.find("ratio").text
    ecosrc = root.find("ecosrc").text
    pourc_eco = root.find("pourc_eco").text

    # values reference point import
    ref_pt = {}
    for var in ['CN_biodiversity', 'CN_carbon', 'CN_water', 'PNUD_human_footprint', 'bosques_degradated_forest',
        'SERNANP_biological_importance', 'WCS_fonct_connect_A', 'MINAM_water_retention_index']:
        for value in ['opti', 'worst']:
            ref_pt[value + '_' + var] = root.find(value + '_' + var).text
            ref_pt[value + '_' + var] = round(float(ref_pt[value + '_' + var]))
        ref_pt['diff_' + var] = abs(ref_pt['opti' + '_' + var] - ref_pt['worst' + '_' + var])

    # right type
    X = int(X)
    Y = int(Y)
    nb_protect = int(nb_protect)
    ratio = float(ratio)
    rat_protect = float(rat_protect)
    nb_protect_ilp = round(nb_protect*ratio)
    ecosrc = str(ecosrc)
    pourc_eco = float(pourc_eco)

    #import rasters
    Rasters = {}

    for rast in glob.glob("./data/" + DIM + "km/" + DIM + "km*"):
        temp = gdal.Open(rast)
        name = re.search(DIM + 'km_(.+?).tif', rast).group(1)
        Rasters["{0}".format(name)] = np.array(temp.ReadAsArray())

    # import of the mask of not relevant pixels
    Mask, Mask_ds = read_geotiff(
        "results/Masked_data_10km.tif")

    ############################### Model optimization #####################################################

    print("*Model optimization building*")
    ### Creating the problem
    model = pulp.LpProblem('Conservation Planning Peru', LpMaximize)

    ### Definition of decision variables
    protected_pix = pulp.LpVariable.dicts("protected_pix", [(i,j) for i in range(X) for j in range(Y)], cat=LpBinary)

    ### Definition of the objective function

    # variables included: *decomment variables you want and choose weight*
    vars_dict = {}
    vars_dict['CN_biodiversity'] = 1/3 # biodiversity
    vars_dict['SERNANP_biological_importance'] = 1/3 # biodiversity
    vars_dict['WCS_fonct_connect_A'] = 1/3 # biodiversity
    vars_dict['CN_carbon'] = 1 # carbon
    vars_dict['CN_water'] = 1/2 # water
    vars_dict['MINAM_water_retention_index'] = 1/2 # water
    vars_dict['PNUD_human_footprint'] = -0.80 # human pressure
    vars_dict['bosques_degradated_forest'] = -0.20 # human pressure

    # ref point method
    model += pulp.lpSum([vars_dict[var] * (pulp.lpSum(([protected_pix[(i,j)] * Rasters[var][i,j] for i in range(X) for j in range(Y)])) - ref_pt['opti_' + var])/ref_pt['diff_' + var] for var in list(vars_dict.keys())])

    # classical weighted sum
    #model += pulp.lpSum([protected_pix[(i,j)] * pulp.lpSum([Rasters[var][i,j] * vars_dict[var] for var in list(vars_dict.keys())]) for i in range(X) for j in range(Y)])

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
            neighbors = [protected_pix[(i-1,j)], protected_pix[(i+1,j)], protected_pix[(i,j-1)], protected_pix[(i,j+1)]]

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
    labels = ["Unprotected", "Current PCA", "To Conserve (Step 1)"]
    # create a patch for every color
    patches = [mpatches.Patch(color=colors[i], label=labels[i]) for i in range(len(values))]
    # put those patched as legend-handles into the legend
    plt.axis('off')
    plt.legend(handles=patches, bbox_to_anchor=(0.9, 1), loc=2, borderaxespad=0., fontsize=10)
    plt.title("Output_ilp", multialignment='left')
    plt.savefig('./results/Output_ilp.png')
    plt.show()

    ### Raw output data

    # Creation of outputs

    protected_pix_np = np.array([protected_pix[(i,j)].varValue for i in range(X) for j in range(Y)]).reshape((X, Y))

    result_t = (Rasters['DECVAR'] - 1) + protected_pix_np + Rasters['Gouv_PA']
    result_bin = protected_pix_np
    Mask_output = Mask + protected_pix_np
    Mask_output = Mask_output - 1

    # connected components
    structure = np.ones((3, 3))
    labeled, ncomponents = label(result_t, structure)
    print("connected components: ",ncomponents)

    #Export of outputs
    model_tif_arr, model_tif_ds = read_geotiff("./data/" + DIM + "km/" + DIM + "km_DECVAR.tif")
    write_geotiff("results/Protect_data_"+DIM+"km_ilp_"+str(int(rat_protect*100))+".tif", result_t, model_tif_ds)
    write_geotiff("results/Protect_bin_" + DIM + "km_ilp_"+str(int(rat_protect*100))+".tif", result_bin, model_tif_ds)
    write_geotiff("results/Masked_data_"+DIM+"km_ilp.tif", Mask_output, model_tif_ds)