import numpy
import matplotlib.pyplot as plt
import pandas
import pathlib2
from seaborn import color_palette as snscolor
import colorcet as cc


numpy.random.seed (4713)

plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})


AADR_VERSION = "v54.1"
EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}"
# all of them
ANNO_FILE = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.table")




DEFAULT_COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
def getDefaultColor (idx):
    return DEFAULT_COLORS[idx % len(DEFAULT_COLORS)]


DEFAULT_MARKERS = ["o", "v", "^", "<", ">", "s", "p", "P", "*", "X", "D"]
def getDefaultMarker (idx):
    return DEFAULT_MARKERS[idx % len(DEFAULT_MARKERS)]


def plotPCA (theTitle, pcs, evals, categories, pca_plot_file, selectMask=None, plotCentroid=False,axtext=""):


    color_palette = snscolor(cc.glasbey, n_colors=len(categories)).as_hex()
    # for now, otherwise we might get some problems
    assert (selectMask == None)

    # default: select all
    if (selectMask is None):
        selectMask = numpy.ones(pcs.shape[0], dtype=bool)

    # plot PC 1 vs PC 2
    firstAxis = 0
    secondAxis = 1
    # firstAxis = 2
    # secondAxis = 3
    pcOne = pcs[:,firstAxis]
    pcTwo = pcs[:,secondAxis]
    (minOne, maxOne) = (pcOne[selectMask].min(), pcOne[selectMask].max())
    (minTwo, maxTwo) = (pcTwo[selectMask].min(), pcTwo[selectMask].max())
    varOne = evals[firstAxis]/numpy.sum(evals)
    varTwo = evals[secondAxis]/numpy.sum(evals)

    theMarker = numpy.zeros (pcs.shape[0], dtype=str)
    theColor = numpy.array (numpy.repeat ("#000000", pcs.shape[0]))
    # print (theMarker)
    # print (theColor)
    theLegend = []
    legendColor = []
    legendMarker = []
    # categories = numpy.flip (categories)
    for c_i, (catMask, catName, catMarker, catColor) in enumerate(categories):
        catColor = color_palette[c_i]
        print(catColor)
        thisMask = catMask & selectMask
        # do we have anything to plot in this category?
        if (numpy.sum(thisMask) > 0):
            theMarker[thisMask] = catMarker
            theColor[thisMask] = catColor
            # print (catName, catColor)
            theLegend.append (catName.replace("_CAP_", "."))
            legendColor.append (catColor)
            legendMarker.append (catMarker)

    # print (theMarker)
    # print (theColor)

    # plot points in random order
    theOrder = numpy.random.permutation (numpy.arange(pcs.shape[0]))
    if (plotCentroid):
        theAlpha = 0.25
    else:
        theAlpha = 1
    for idx in theOrder:
        plt.plot (pcOne[idx], pcTwo[idx], theMarker[idx], color=theColor[idx], markersize=1, alpha=theAlpha)

    # plot centroids after this
    for c_i, (catMask, catName, catMarker, catColor) in enumerate(categories):
        catColor = color_palette[c_i]
        thisMask = catMask & selectMask
        # do we have anything to plot in this category?
        if (numpy.sum(thisMask) > 0):
            # plot centroid already here
            if (plotCentroid):
                centOne = numpy.mean(pcOne[thisMask])
                centTwo = numpy.mean(pcTwo[thisMask])
                plt.plot (centOne, centTwo, catMarker, color=catColor, markersize=3.75, alpha=1)

    plt.legend (theLegend, fontsize=7, labelspacing=.2, handlelength=1.5, handleheight=.5, handletextpad=.4, borderpad=.2, borderaxespad=.2)
    ax = plt.gca()
    ax.text(-.2, .97, rf"$\bf{{{axtext}}}$", fontsize=13, transform=ax.transAxes)
    leg = ax.get_legend()
    for i in numpy.arange(len(legendColor)):
        leg.legendHandles[i].set_color(legendColor[i])
        leg.legendHandles[i].set_marker(legendMarker[i])
        leg.legendHandles[i].set_alpha(1)
        leg.legendHandles[i].set_markersize(4)

    plt.xlabel (f"PC1 ({varOne*100:.2f}%)")
    plt.ylabel (f"PC2 ({varTwo*100:.2f}%)")

    # coord grid
    if ((minOne < 0) and (maxOne > 0)):
        plt.vlines (0, minTwo, maxTwo, linestyles='--',colors='black', lw=0.5)
    if ((minTwo < 0) and (maxTwo > 0)):
        plt.hlines (0, pcOne.min(), pcOne.max(), linestyles='--',colors='black', lw=0.5)

    #plt.title (theTitle)
    plt.savefig (pca_plot_file, bbox_inches="tight")
    plt.clf()


def allPcaPlots():

    # make plots a bit bigger for now

    # we need the publications for refined labels
    anno = pandas.read_csv (ANNO_FILE, sep='\t', low_memory=False)
    # compress names of publications
    short_pubs = numpy.array ([x.split()[0][0] + x.split()[0][-2:] for x in anno['Publication']])

    # go through all combinations
    # for reference in ['broad', 'europe']:
    for reference in ['broad', 'europe', 'gbr_ceu']:
        for to_shrink in ['shrinkage', 'noshrinkage']:
            for genotyping in ['capture_only', 'capture_SG']:
                axtext = ""
                if to_shrink == "shrinkage" and genotyping == "capture_only":
                    if reference == "broad":
                        axtext = "A"
                    elif reference == "europe":
                        axtext = "B"
                plt.figure(figsize=(3.1, 3.1), layout="constrained")
                # load specific pca files
                eval_file = pathlib2.Path (f"{EXTRACTED_PREFIX}_{genotyping}_{reference}_pca_{to_shrink}.eval")
                evec_file = pathlib2.Path (f"{EXTRACTED_PREFIX}_{genotyping}_{reference}_pca_{to_shrink}.evec")

                evals = numpy.array (pandas.read_csv (eval_file, header=None)).flatten()
                evecFrame = pandas.read_csv (evec_file, delim_whitespace=True)

                # get the data in convenient format
                pcs = numpy.array(evecFrame.iloc[:,:-1])
                popLabels = numpy.array (evecFrame.iloc[:,-1])
                ids = numpy.array (evecFrame.index)

                # leverage pub labels from just ancients to whole pca sample
                tmp_anno_ids = list(anno['Genetic_ID'])
                leveraged_pubs = numpy.array(['' if (x not in tmp_anno_ids) else (short_pubs[tmp_anno_ids.index(x)]) for x in ids])

                # have some categories for plotting
                categories = []

                # publication categories
                colorIdx = 0
                for pop in sorted(set(popLabels)):
                    thisPopMask = (popLabels == pop)
                    if (pop in ["FOCAL_ANCIENT_SG", "FOCAL_ANCIENT_CAPTURE"]):
                        short_pop = "ANC_CAP" if "CAPTURE" in pop else "ANC_SG"
                        # special
                        for pub in set(short_pubs):
                            thisPubMask = (leveraged_pubs == pub)
                            if (thisPopMask & thisPubMask).sum() > 0:
                                categories.append((thisPopMask & thisPubMask, short_pop + "_" + pub,
                                                   getDefaultMarker(colorIdx), getDefaultColor(colorIdx)))
                                colorIdx += 1
                    else:
                        # regular
                        categories.append ((thisPopMask, pop, getDefaultMarker (colorIdx), getDefaultColor (colorIdx)))
                        colorIdx += 1

                # plot a pca with all samples
                pca_plot_file = pathlib2.Path (f"{EXTRACTED_PREFIX}_{genotyping}_{reference}_pca_{to_shrink}.pdf")
                if (reference == 'broad'):
                    plotCentroid = False
                else:
                    plotCentroid = True
                plotPCA (f"PCA: {genotyping}, {reference}, {to_shrink}", pcs, evals, categories, pca_plot_file, plotCentroid=plotCentroid, axtext=axtext)

                # # and for the broad analysis, we also plot a zoom on ancient + europe
                # if (reference == 'broad'):
                #     zoomed_pops = ['FOCAL_ANCIENT', 'CEU.SG', 'GBR.SG', 'FIN.SG', 'TSI.SG', 'IBS.SG']
                #     selectMask = numpy.array ([x in zoomed_pops for x in popLabels])
                #     # try to omit outliers (certainly for the zoom)
                #     selectMask[(pcs[:,0] > 0) & (popLabels == 'FOCAL_ANCIENT')] = False
                #     print (ids[(pcs[:,0] > 0) & (popLabels == 'FOCAL_ANCIENT')])
                #     pca_plot_file = pathlib2.Path (f"{EXTRACTED_PREFIX}_{reference}_pca_{to_shrink}_zoom.pdf")
                #     plotPCA (f"PCA: {reference}, {to_shrink}, zoom", pcs, evals, categories, pca_plot_file, selectMask=selectMask)


def main():

    # make some plots
    allPcaPlots ()


if __name__ == "__main__":
    main()