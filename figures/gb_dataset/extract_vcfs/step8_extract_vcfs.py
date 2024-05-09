import pandas
import pathlib2
import collections
import os
import subprocess




AADR_VERSION = "v54.1"
# we take the 1240K version here
AADR_ROOT = pathlib2.Path (f'AADR/{AADR_VERSION}_1240K/{AADR_VERSION}_1240K_public')
AADR_SNP = pathlib2.Path (f'{AADR_ROOT}.snp')
GENETIC_ID = 'Genetic_ID'

# file specifiying individuals to extract
EXTRACTED_PREFIX = f"extracted/GB_{AADR_VERSION}_capture_SG"
INDIVIDUALS_FILE = pathlib2.Path (f'{EXTRACTED_PREFIX}_inds.table')

# external scripts
EIGENSTRAT_CONVERSION_SCRIPT = pathlib2.Path ('gdc/eigenstrat2vcf.py')
CORRECT_HAPLO_ENCODING_SCRIPT = pathlib2.Path ('diplo_to_haplo_vcf.py')




def extractTimeSeries():

    # load the individuals and make an individuals file for conversion
    conversionIndFile = pathlib2.Path (f'{EXTRACTED_PREFIX}_inds.conv')
    individualsFrame = pandas.read_csv (INDIVIDUALS_FILE, sep='\t')
    ofs = open (conversionIndFile, 'w')
    for thisInd in individualsFrame[GENETIC_ID]:
        ofs.write (f"{thisInd}\n")
    ofs.close()

    # load all the snps
    snpFrame = pandas.read_csv (AADR_SNP, delim_whitespace=True)
    print (snpFrame.shape)
    chromHist = collections.Counter (snpFrame.iloc[:,1])
    # make sure all chromosomes accounted for
    print (chromHist.keys())
    assert (len(chromHist.keys()) <= 24), len(chromHist.keys())
    assert (min(chromHist.keys()) == 1), min(chromHist.keys())
    assert (max(chromHist.keys()) == 24), max(chromHist.keys())

    # one vcf-file for each chromosome
    for c in chromHist.keys():

        chromName = str(c)
        # 23 is X
        if (c == 23):
            chromName = 'X'
        # 24 is Y
        elif (c == 24):
            chromName = 'Y'
        else:
            pass

        # make a snpfile for eigenstrat
        conversionSnpFile = pathlib2.Path (f'{EXTRACTED_PREFIX}_c{chromName}.snps')
        ofs = open (conversionSnpFile, 'w')
        thisSnpFrame = snpFrame.loc[snpFrame.iloc[:,1] == c]
        for thisSnp in thisSnpFrame.iloc[:,0]:
            ofs.write (f"{thisSnp}\n")
        ofs.close()

        # prepare the output file
        outputDiploVCF = pathlib2.Path (f'{EXTRACTED_PREFIX}_c{chromName}.diplo_vcf')

        # put eigentstrat command together
        stratCmd = f"python {EIGENSTRAT_CONVERSION_SCRIPT} -r {AADR_ROOT} -i {conversionIndFile} -s {conversionSnpFile} >{outputDiploVCF}"
        print (stratCmd)

        # extract it
        subprocess.run([stratCmd], shell=True, text=True)

        # clean up
        os.remove (conversionSnpFile)

        # and make it a proper vcf with haploid calls encoded correctly
        outputVCF = pathlib2.Path (f'{EXTRACTED_PREFIX}_c{chromName}.vcf')

        haploCmd = f"python {CORRECT_HAPLO_ENCODING_SCRIPT} {outputDiploVCF} > {outputVCF}"
        print (haploCmd)

        # convert it to pseudo-haploid
        subprocess.run([haploCmd], shell=True, text=True)

        # clean up
        os.remove (outputDiploVCF)


def main():
    extractTimeSeries ()


if __name__ == "__main__":
    main()