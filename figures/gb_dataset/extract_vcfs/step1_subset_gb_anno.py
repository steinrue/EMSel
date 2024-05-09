import pandas
import numpy
import collections
import pathlib2



# where is the AADR?
AADR_VERSION = "v54.1"
ANNO_FILE_1240K = pathlib2.Path (f"AADR/{AADR_VERSION}_1240K/{AADR_VERSION}_1240K_public.anno")
ANNO_FILE_HO = pathlib2.Path (f"AADR/{AADR_VERSION}_HO/{AADR_VERSION}_HO_public.anno")
IND_FILE_1240K = pathlib2.Path (f"AADR/{AADR_VERSION}_1240K/{AADR_VERSION}_1240K_public.ind")
IND_FILE_HO = pathlib2.Path (f"AADR/{AADR_VERSION}_HO/{AADR_VERSION}_HO_public.ind")
# where to put the output
OUTPUT_DIR = pathlib2.Path ('extracted')
EXTRACTED_PREFIX = f"{OUTPUT_DIR}/GB_{AADR_VERSION}"
# we want to extract capture + SG as a vcf, we can then filter further at a later step
OUTPUT_TABLE_SG = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.table")
# the mock files for PCA also one for capture + SG
OUTPUT_IND_SG = pathlib2.Path (f"{EXTRACTED_PREFIX}_capture_SG_pre_pca_inds.ind")
# reference populations, take the SG versions
REFERENCE_POPS = ['GBR.SG', 'FIN.SG', 'IBS.SG', 'TSI.SG', 'CEU.SG', 'YRI.SG', 'LWK.SG', 'CHB.SG', 'JPT.SG']




# some constants for the field names
# v54
POLITICAL_ENTITY = 'Political_Entity'
GENETIC_ID = 'Genetic_ID'
TARGET_COVERAGE = '1240k_coverage'
SNPS_HIT = 'SNPs_hit_1240k'
NO_RELATIVE_ENTRIES = ['n/a (no relatives detected)', 'n/a(norelativesdetected)', 'n/a (No relatives detected)']




def loadAnnotationFile (annoFile):

    # read it
    anno = pandas.read_csv (annoFile, sep='\t', low_memory=False)

    # reduce the column names
    # fuse first two names, or take the only one
    newColumns = []
    for c in anno.columns:
        fields = c.split()
        if (len(fields) <= 1):
            # just this one
            newColumns.append (fields[0])
        elif ((fields[0] == "SNPs") and (fields[1] == "hit")):
            # this one a bit special
            if any([x == '1240k' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}_1240k")
            elif any([x == 'HO' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}_HO")
            elif any([x == '3.2M' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}_3.2M")
            elif any([x == 'non-padded' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}_non-padded")
            else:
                # this should be ok for version 50.0?
                newColumns.append (f"{fields[0]}_{fields[1]}")
        elif ((fields[0] == "Y") and (fields[1] == "haplogroup")):
            # this one's special too
            if any([x == 'terminal' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}_terminal")
            elif any([x == 'ISOGG' for x in fields]):
                if (any([x == 'curation' for x in fields])):
                    newColumns.append (f"{fields[0]}_{fields[1]}_ISOGG_curation")
                else:
                    newColumns.append (f"{fields[0]}_{fields[1]}_ISOGG")
            else:
                assert (False)
        elif ((fields[0] == "Xcontam") and (fields[1] == "ANGSD")):
            # many special snowflakes
            if any([x == 'SNPs' for x in fields]):
                newColumns.append (f"{fields[0]}_{fields[1]}SNPs")
            elif any([x == 'MOM' for x in fields]):
                if (any([x == 'point' for x in fields])):
                    newColumns.append (f"{fields[0]}_{fields[1]}_point")
                elif (any([x == 'Z-score' for x in fields])):
                    newColumns.append (f"{fields[0]}_{fields[1]}_Z-score")
                elif (any([x == 'CI' for x in fields])):
                    newColumns.append (f"{fields[0]}_{fields[1]}_CI")
                else:
                    assert (False)
            else:
                assert (False)
        else:
            # fuse first two
            newColumns.append (f"{fields[0]}_{fields[1]}")

    # set the new names
    anno.columns = newColumns

    # make sure we don't have duplicates in the new names
    # print(anno.columns)
    assert (numpy.all ([x <= 1 for x in collections.Counter(anno.columns).values()]))

    return anno


def removeDuplicates (annoFrame):

    # remove duplicates
    # make sure that we don't have duplicated Master_ID's
    counterMasterID = collections.Counter (annoFrame['Master_ID'])
    print ("[REMOVE_DUPLICATES]")
    print (annoFrame.shape)
    print ("=====")
    genIdsToKeep = []
    for (k,v) in counterMasterID.items():
        # see if this one is unique
        if (v == 1):
            theseIds = list(annoFrame.loc[annoFrame['Master_ID'] == k][GENETIC_ID])
            assert (len(theseIds) == 1)
            genIdsToKeep.append (theseIds[0])
        if (v != 1):
            # not unique, so show the culprit
            print (f"{k} [Master_ID] is present {v} times.")
            theseInds = annoFrame.loc[annoFrame['Master_ID'] == k]
            theseIds = list(theseInds[GENETIC_ID])
            theseCoverage = [-1 if x == ".." else float(x) for x in list(theseInds[TARGET_COVERAGE])]
            greatestHits = list(theseInds[SNPS_HIT])
            print (f"Genetic_ID's: {theseIds}")
            print (f"coverage: {theseCoverage}")
            print (f"1240K snps covered: {greatestHits}")
            print (f"dates: {list(theseInds['Date_mean'])}")

            # decide based on most SNPs covered
            idxToKeep = numpy.argmax(greatestHits)
            print (f"{theseIds[idxToKeep]} taken cause more SNPs covered")

            genIdsToKeep.append (theseIds[idxToKeep])
            print ("=====")

    # take only unique ones
    finalIndividuals = annoFrame.loc[[ x in genIdsToKeep for x in annoFrame[GENETIC_ID]]]

    # see if all unique now
    counterMasterID = collections.Counter (finalIndividuals['Master_ID'])
    assert ((numpy.array(list(counterMasterID.values())) == 1).all())
    print ("[DUPLICATES_REMOVED]")
    print ("[FINAL_COUNT]")
    print (finalIndividuals.shape)

    return finalIndividuals


def filterTime (annoFrame, minTimeBPyears, maxTimeBPyears):

    print (f"{annoFrame.shape} individuals")

    assert (minTimeBPyears  >= 0)
    assert (minTimeBPyears <= maxTimeBPyears)

    timeMask = (minTimeBPyears <= annoFrame['Date_mean']) & (annoFrame['Date_mean'] <= maxTimeBPyears)

    finalFrame = annoFrame.loc[timeMask]

    print (f"{finalFrame.shape} individuals in time frame [{minTimeBPyears}, {maxTimeBPyears}]")
    presentMask = (annoFrame['Date_mean'] == 0)
    ancientMask = (annoFrame['Date_mean'] > 0)
    print (f"\t{numpy.sum(timeMask & presentMask)} at present")
    print (f"\t{numpy.sum(timeMask & ancientMask)} in past")

    return finalFrame


def keepOnlyCapture (annoFrame):

    print (collections.Counter(annoFrame['Data_source']))
    print (annoFrame.shape)

    # both spellings exist
    captureMask = numpy.isin (annoFrame['Data_source'], ["1240K", "1240k"])
    # check whether names are somewhat consitent with labeling
    assert (not numpy.any ([((".SG" in x) or (".DG" in x)) for x in annoFrame[captureMask][GENETIC_ID]]))
    # this should be all the other ones
    shotgunMask = numpy.isin (annoFrame['Data_source'], ['Shotgun', 'Shotgun.diploid'])
    # check whether names are somewhat consitent with labeling
    assert (numpy.all ([((".SG" in x) or (".DG" in x)) for x in annoFrame[shotgunMask][GENETIC_ID]]))
    # make sure we accounted for all individuals
    assert  (annoFrame.shape[0] == numpy.sum (captureMask) + numpy.sum (shotgunMask))

    # only keep the capture individuals
    captureFrame = annoFrame[captureMask]

    print (captureFrame.shape)

    return captureFrame


def removeIslands (annoFrame):

    print ("[REMOVE_ISLANDS]")

    # remove island individuals
    # these ids were obtained by manual inspection of the time transect in v54.1
    islandIds = ['VK171', 'I5367', 'KD005', 'I2824', 'I2655']
    # and a lot of individuals north of 58.8 degrees (mostly orkney)
    islandMask = numpy.array ([x in islandIds for x in annoFrame['Master_ID']]) | (numpy.array([float(x) if x != '..' else 0 for x in annoFrame['Lat.']]) >= 58.5)
    print ("Individuals on islands:")
    print (annoFrame[islandMask][['Master_ID', 'Genetic_ID', 'Date_mean', 'Locality', 'Publication']])
    newAnno = annoFrame[~islandMask]
    # print (newAnno[numpy.array([float(x) if x != '..' else 0 for x in newAnno['Lat.']]) > 59][['Genetic_ID', 'Lat.']])

    print ("[ISLANDS_REMOVED]")

    return newAnno


def extractGB (annoFrame):

    print ("[EXTRACT]")
    # take all that have political entity UK (subsequently referenced as GB)
    entityMask = (annoFrame[POLITICAL_ENTITY] == 'United Kingdom')

    # that have PASS for the filter column
    # we take strict PASS only here
    # this actually also removes the diploid version of the 1000g, since they don't have a clear PASS
    passMask = (annoFrame['ASSESSMENT'] == 'PASS')

    # remove relatives
    noRelativeMask = numpy.isin (annoFrame['Family_ID'], NO_RELATIVE_ENTRIES)

    # take these individuals for now
    GBMask = entityMask & passMask & noRelativeMask
    annoFrame = annoFrame.loc[GBMask]

    print (f"{numpy.sum(GBMask)} individuals in GB that pass filter (and no relatives)")

    # remove individuals outside of time frame
    annoFrame = filterTime (annoFrame, 0, 4450)

    print (annoFrame.shape)

    # APPLY THIS FILTERING LATER
    # # keep only capture data
    # annoFrame = keepOnlyCapture (annoFrame)

    # print (annoFrame.shape)

    # remove duplicates
    annoFrame = removeDuplicates (annoFrame)

    print (annoFrame.shape)
    
    # remove islands
    annoFrame = removeIslands (annoFrame)

    print (annoFrame.shape)

    print ("[EXTRACT_DONE]")

    return annoFrame


def writeAnnotation (annoGB, output_file):
    # just write it
    annoGB.to_csv (output_file, sep='\t', index=False)


def createMockIndFile (inputHoIndFile, outputMockFile, focalFrame):

    print ("[CREATE_MOCK_IND_FILE]")

    # read the original ind file
    aadrHoIndFrame = pandas.read_csv (inputHoIndFile, delim_whitespace=True, header=None, names=['ID', 'sex', 'location'])

    # clear all previous group data
    # except for the reference populations that are suposed to stay
    referenceMask = numpy.array([x in REFERENCE_POPS for x in aadrHoIndFrame['location']])
    # don't include anything but the specific reference populations we want
    aadrHoIndFrame['location'][~referenceMask] = 'NO'
    # we fill some ancients later

    # make sure data sources are good
    assert (set(focalFrame['Data_source']).issubset(set(['Shotgun', 'Shotgun.diploid', '1240K', '1240k'])))

    # and get the relevant masks
    shotgunMask = numpy.isin (focalFrame['Data_source'], ['Shotgun', 'Shotgun.diploid'])
    captureMask = numpy.isin (focalFrame['Data_source'], ["1240K", "1240k"])
    # make sure nothing missed
    assert (not numpy.any(shotgunMask & captureMask))
    assert (numpy.all (shotgunMask | captureMask))

    presentGBRMask = (focalFrame['Group_ID'] == 'GBR.SG')
    # see that these are the only guys at present
    assert (numpy.all (numpy.isclose (focalFrame['Date_mean'][presentGBRMask], 0)))
    assert (numpy.all (focalFrame['Date_mean'][~presentGBRMask] > 0))
    # and also all shotgun
    assert (numpy.all (shotgunMask | ~presentGBRMask))

    # and all part of the reference
    refFrame = aadrHoIndFrame[referenceMask]
    assert (set(refFrame[refFrame['location'] == 'GBR.SG']['ID']) == set(focalFrame[presentGBRMask][GENETIC_ID]))

    # now fill in entries for our focal individuals
    # but in two versions
    captureIds = numpy.array(focalFrame[captureMask][GENETIC_ID])
    sgIds = numpy.array(focalFrame[shotgunMask][GENETIC_ID])
    assert (len(sgIds) + len(captureIds) == focalFrame.shape[0]), (len(sgIds), len(captureIds), focalFrame.shape[0])
    assert (set(captureIds).intersection(set(sgIds)) == set([]))

    # get the right masks
    captureMask = [(x in captureIds) for x in aadrHoIndFrame['ID']]
    sgMask = [(x in sgIds) for x in aadrHoIndFrame['ID']]
    # make sure we have them all
    assert (numpy.sum(captureMask) == len(captureIds)), (numpy.sum(captureMask), len(captureIds))
    assert (numpy.sum(sgMask) == len(sgIds)), (numpy.sum(sgMask), len(sgIds))

    # and rewrite the frames, but leave the reference untouched (which should include GBR.SG)
    aadrHoIndFrame['location'][captureMask & ~referenceMask] = 'FOCAL_ANCIENT_CAPTURE'
    aadrHoIndFrame['location'][sgMask & ~referenceMask] = 'FOCAL_ANCIENT_SG'

    print (collections.Counter (aadrHoIndFrame['location']))

    # write the modified ind-file
    aadrHoIndFrame.to_csv (outputMockFile, sep='\t', header=False, index=False)

    print ("[CREATE_MOCK_IND_FILE_DONE]")


def main ():

    # load the 1240K anno file
    anno = loadAnnotationFile (ANNO_FILE_1240K)

    # and also the HO file
    annoHO = loadAnnotationFile (ANNO_FILE_HO)
    # and make sure that 1240K anno is in HO anno
    assert (set(anno[GENETIC_ID]).issubset(set(annoHO[GENETIC_ID])))

    # also check whether the individuals files match
    indFrame = pandas.read_csv (IND_FILE_1240K, delim_whitespace=True, header=None, names=['ID', 'sex', 'location'])
    indFrameHO = pandas.read_csv (IND_FILE_HO, delim_whitespace=True, header=None, names=['ID', 'sex', 'location'])
    # could also check whether data matches, but for now only check if invidivuals are there
    assert (set(indFrame["ID"]).issubset(set(indFrameHO["ID"])))

    # extract the GB individuals that we want to focus on
    annoGB = extractGB (anno)

    # ACTUALLY NOT TRUE ANYMORE
    # # for this particular dataset, we should have only ancient individuals
    # assert (annoGB['Date_mean'].min() > 0)

    # just check some things
    print (f"{annoGB.shape[0]} GB individuals in time frame")
    presentMask = (annoGB['Date_mean'] <= 0)
    ancientMask = (annoGB['Date_mean'] > 0)
    shotgunMask = numpy.isin (annoGB['Data_source'], ['Shotgun', 'Shotgun.diploid'])
    captureMask = numpy.isin (annoGB['Data_source'], ["1240K", "1240k"])
    print (collections.Counter(annoGB['Publication']))
    print (f"{numpy.sum (presentMask)} individuals at present")
    # all at present should be shotgun
    assert (numpy.all (presentMask == (presentMask & shotgunMask)))
    print ('\tall shotgun')
    print (f"{numpy.sum (ancientMask)} ancient individuals")
    print (f"\t{annoGB[ancientMask & shotgunMask].shape[0]} shotgun; {annoGB[ancientMask & captureMask].shape[0]} capture")

    # # create the directory for the output
    # os.makedirs (OUTPUT_DIR)

    # and write the new table file
    writeAnnotation (annoGB, OUTPUT_TABLE_SG)

    # now, we also need to modify the *.ind file
    createMockIndFile (IND_FILE_HO, OUTPUT_IND_SG, annoGB)


if __name__ == "__main__":
    main()