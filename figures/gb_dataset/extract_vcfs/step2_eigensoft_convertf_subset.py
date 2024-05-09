import pathlib2
import numpy
import pandas
import subprocess




AADR_VERSION = "v54.1"
# use HO for PCA
AADR_GENO = pathlib2.Path (f"AADR/{AADR_VERSION}_HO/{AADR_VERSION}_HO_public.geno")
AADR_SNP = pathlib2.Path (f"AADR/{AADR_VERSION}_HO/{AADR_VERSION}_HO_public.snp")

EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}"
# some mock inds
MOCK_FILE = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.ind")
MOCK_VERSIONS = {
    'capture_only' : ['FOCAL_ANCIENT_CAPTURE'],
    'capture_SG' : ['FOCAL_ANCIENT_CAPTURE', 'FOCAL_ANCIENT_SG'],
}
# reference pops for different pcas
REFERENCE_POPS = {
    'gbr_ceu' : ['GBR.SG', 'CEU.SG'],
    'europe' : ['GBR.SG', 'FIN.SG', 'IBS.SG', 'TSI.SG', 'CEU.SG'],
    'broad' : ['GBR.SG', 'FIN.SG', 'IBS.SG', 'TSI.SG', 'CEU.SG', 'YRI.SG', 'LWK.SG', 'CHB.SG', 'JPT.SG'],
}




def callEigensoftConvert ():

    # make sure the mock-file is as expected
    allowedLabels = numpy.concatenate ((REFERENCE_POPS['broad'], ['NO'], MOCK_VERSIONS['capture_SG']))
    mockIndFrame = pandas.read_csv (MOCK_FILE, sep='\t', header=None)
    assert (set(mockIndFrame[2]) == set(allowedLabels)), (set(mockIndFrame[2]), set(allowedLabels))

    # iterate through references
    for (thisLabel, thisReference) in REFERENCE_POPS.items():

        # iterate through versions of mock individuals
        for (thisMockLabel, thisMockPopulations) in MOCK_VERSIONS.items():

            print (f"===== ({thisMockLabel}, {thisLabel})")

            # some specific files
            mockPrefix = f"{EXTRACTED_PREFIX}_{thisMockLabel}_{thisLabel}"
            to_keep_file = pathlib2.Path (f"{mockPrefix}_pops_to_keep.tsv")
            convertf_param_file = pathlib2.Path (f"{mockPrefix}_convertf.param")
            pca_geno_file = pathlib2.Path (f"{mockPrefix}_pca.geno")
            pca_snp_file = pathlib2.Path (f"{mockPrefix}_pca.snp")
            pca_ind_file = pathlib2.Path (f"{mockPrefix}_pca.ind")

            # prepare the files of populations that we want to keep
            ofs = open (to_keep_file, 'w')
            for pop in numpy.concatenate ((thisReference, thisMockPopulations)):
                ofs.write (pop + '\n')
            ofs.close()

            # then prepare the parameter file for convertf
            ofs = open(convertf_param_file, 'w')
            ofs.write (
                f"genotypename: {AADR_GENO}\n" + \
                f"snpname: {AADR_SNP}\n" + \
                f"indivname: {MOCK_FILE}\n" + \
                "outputformat: EIGENSTRAT\n" + \
                f"indivoutname: {pca_ind_file}\n" + \
                f"snpoutname: {pca_snp_file}\n" + \
                f"genotypeoutname: {pca_geno_file}\n" + \
                f"poplistname: {to_keep_file}\n" + \
                "numoutlieriter: 0\n"
            )
            ofs.close()

            # and run convertf
            convertfCmd = f'convertf -p {convertf_param_file}'
            print (convertfCmd)
            # run the command
            subprocess.run([convertfCmd], shell=True, text=True)


def main():
    callEigensoftConvert ()


if __name__ == "__main__":
    main()