import pandas
import pathlib2
import numpy
import matplotlib.pyplot as plt
import collections




AADR_VERSION = "v54.1"
EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}"
ANNO_FILE = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.table")
GENOTYPING = {
    'capture_only' : ['1240K', '1240k'],
    'capture_SG' : ['1240K', 'Shotgun', '1240k'],
}
ALL_GENOTYPES_SET = set([val for dic in GENOTYPING.values() for val in dic])
SAMPLE_SIZES_TABLES = {}
SAMPLE_SIZES_PDFS = {}
for geno in GENOTYPING.keys():
    SAMPLE_SIZES_TABLES[geno] = pathlib2.Path (f"{EXTRACTED_PREFIX}_{geno}_sample_sizes.table")
    SAMPLE_SIZES_PDFS[geno] = pathlib2.Path (f"{EXTRACTED_PREFIX}_{geno}_sample_sizes.pdf")




def plotSampleSizes (annoFrame, sampleSizeTable, sampleSizePdf):

    plt.rcParams.update({'font.size': 9,
                     'text.usetex': False,
                     'font.family': 'serif',
                     'font.serif': 'cmr10',
                     'mathtext.fontset': 'cm',
                     'axes.unicode_minus': False,
                     'axes.formatter.use_mathtext': True,})
    plt.locator_params(axis='x', nbins=5)

    # when are the samples, in generations (Moorjani et al., 2016)
    YEARS_PER_GEN = 28
    sampleCounter = collections.Counter ((numpy.array(annoFrame['Date_mean'])/YEARS_PER_GEN).astype(int))
    # plot it
    sampleSizes = pandas.DataFrame (columns=['time', 'sample_size'])
    for (k,v) in sampleCounter.items():
        # plot it
        # plt.plot ([k, k], [0, v], '-')
        ypb = -k*YEARS_PER_GEN
        plt.plot ([ypb, ypb], [0, v], '-', lw=2.5)
        # and record it
        newRow = {'time': k, 'sample_size': v}
        # sampleSizes = sampleSizes.append (newRow, ignore_index = True)
        # needs pandas 2.0
        sampleSizes = pandas.concat([sampleSizes, pandas.DataFrame([newRow])], ignore_index=True)
    # plt.yscale ("log")
    plt.xlim(-max(list(sampleCounter.keys()))*YEARS_PER_GEN+YEARS_PER_GEN, 3*YEARS_PER_GEN)
    plt.xlabel ("Years before present")
    plt.ylabel ("Sample size")
    plt.title ("Sample sizes for ancient dataset")

    plt.tight_layout()

    plt.savefig (sampleSizePdf)
    plt.clf()

    # also save all the sizes (including present)
    sampleSizes = sampleSizes.sort_values ('time')
    sampleSizes.to_csv (sampleSizeTable, sep='\t', index=False)

def main():

    # load anno file
    anno = pandas.read_csv (ANNO_FILE, sep='\t', low_memory=False)
    assert (set(anno['Data_source']).issubset(ALL_GENOTYPES_SET))

    # go through two different modes of genotyping
    for (genoLabel, genoList) in GENOTYPING.items():

        # get a mask for only these individuals
        genoMask = [(x in genoList) for x in anno['Data_source']]

        # plot smaples sizes for only the individuals with right genotyping
        plotSampleSizes (anno[genoMask], SAMPLE_SIZES_TABLES[genoLabel], SAMPLE_SIZES_PDFS[genoLabel])


if __name__ == "__main__":
    main()