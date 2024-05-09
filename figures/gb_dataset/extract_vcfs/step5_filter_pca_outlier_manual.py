import numpy
import pandas
import pathlib2



AADR_VERSION = "v54.1"
EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}_capture_SG"
INPUT_TABLE = pathlib2.Path (f'{EXTRACTED_PREFIX}_pre_pca_inds.table')
OUTPUT_TABLE = pathlib2.Path (f'{EXTRACTED_PREFIX}_inds.table')
# # these ones look like outliers in the PCA plot
# GENETIC_IDS_TO_REMOVE = ['I11570']
# for now no removal
GENETIC_IDS_TO_REMOVE = []




def removeIds():
    
    # load old anno
    anno = pandas.read_csv (INPUT_TABLE, sep='\t', low_memory=False)
    print (anno.shape)

    # remove some genetic ids
    removeMask = numpy.array ([x in GENETIC_IDS_TO_REMOVE for x in anno['Genetic_ID']])
    newAnno = anno.loc[~removeMask]
    print (newAnno.shape)

    # store new anno
    newAnno.to_csv (OUTPUT_TABLE, sep='\t', index=False)
    
    
def main():

    removeIds()


if __name__ == "__main__":
    main()