import pandas
import pathlib2
import numpy
import matplotlib.pyplot as plt
import cartopy




AADR_VERSION = "v54.1"
EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}"
ANNO_FILE = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.table")
GENOTYPING = {
    'capture_only' : ['1240K', '1240k'],
    'capture_SG' : ['1240K', 'Shotgun', '1240k'],
}
ALL_GENOTYPES_SET = set([val for dic in GENOTYPING.values() for val in dic])
MAP_FILES = {}
for geno in GENOTYPING.keys():
    MAP_FILES[geno] = pathlib2.Path (f"{EXTRACTED_PREFIX}_{geno}_map.pdf")




def plotMap (annoFrame, mapFile):

    plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})
    # fig = plt.figure(figsize=(18,18),layout="constrained")
    fig = plt.figure(figsize=(5,5),layout="constrained")

    # some might not have lat or long, so don't plot them
    latLongMask = (annoFrame['Long.'] != '..') & (annoFrame['Lat.'] != '..')
    
    print (annoFrame.shape)
    print (numpy.sum (~latLongMask))
    print ("Individuals without coordinates:")
    print (annoFrame[~latLongMask][['Genetic_ID', 'Date_mean', 'Locality', 'Publication']])

    plotIndividuals = annoFrame.loc[latLongMask]

    age_years = numpy.array(plotIndividuals['Date_mean']).astype(float)
    longs = numpy.array (plotIndividuals['Long.']).astype(float)
    lats = numpy.array (plotIndividuals['Lat.']).astype(float)
    # use master id because of potential duplicates
    labels = numpy.array (plotIndividuals['Master_ID'])

    ax = plt.axes (projection=cartopy.crs.PlateCarree())

    # map boundaries, in W-E S-N
    ax.set_extent([-9, 4, 49, 62], crs=cartopy.crs.PlateCarree())
    # ax.set_global()

    ax.coastlines (lw=1)
    # needs cartopy 0.21
    im = ax.scatter (longs, 
            lats, 
            c=age_years,
            transform=cartopy.crs.PlateCarree(), s=30, alpha=0.75)

    # # labelling helpful for island removal
    # for (idx, thisLabel) in enumerate(labels):
    #     ax.annotate (thisLabel, (longs[idx], lats[idx]), transform=cartopy.crs.PlateCarree())

    ax.set_aspect('auto')
    # ax.outline_patch.set_visible(False)
    ax.add_feature(cartopy.feature.BORDERS, linestyle=':', linewidth=1.5)

    fig.canvas.draw()
    plt.colorbar(im, label='Years before present', shrink=.75)

    #plt.title ("Location of ancient individuals")
    #plt.tight_layout()

    print (f"save to: {mapFile}")
    plt.savefig (mapFile)
    plt.clf()


def main():

    # load anno file
    # one combined file for both genotyping methods
    anno = pandas.read_csv (ANNO_FILE, sep='\t', low_memory=False)
    assert (set(anno['Data_source']).issubset(ALL_GENOTYPES_SET))

    # go through two different modes of genotyping
    for (genoLabel, genoList) in GENOTYPING.items():

        # get a mask for only these individuals
        genoMask = [(x in genoList) for x in anno['Data_source']]

        # plot map for only the individuals with right genotyping
        plotMap (anno[genoMask], MAP_FILES[genoLabel])


if __name__ == "__main__":
    main()