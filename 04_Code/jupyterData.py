import glob
import os
import pathlib
import pandas as pd
import numpy as np

def scanFolders(paths):

    # Initialize folders list of same size as paths list
    folders = list(range(len(paths)))

    # Loop over all folder paths and extract folder names
    for i, path in enumerate(paths):
        folders[i] = pathlib.PurePath(path).name
        
    # Sort list of folders alphanumerically
    folders.sort()

    # Create dictionary (keys and corresponding values inserted later)
    csvFiles = {}

    # Create list of files (files inserted later)
    fileList = []

    # Loop over all paths and corresponding foldernames
    for path, folder in zip(paths, folders):
        # Empty file list at each new iteration
        fileList = []
    
        # Loop over all CSV-files in the corresponding CSV subfolder
        for file in glob.glob(os.path.join(path, '01_CSV', '*.csv')):
            fileList.append(file)
        
        # Sort current file list
        fileList = sorted(fileList)
        
        # Store the current file list under a key (= current foldername) in the dictionary "csvFiles"
        csvFiles[folder] = fileList
        
    print('The following folders have been detected:', '\n', '\n', folders, '\n')
    print('The CSV-files of each folder are stored in the dictionary "csvFiles"!', '\n')
    print('Note: The corresponding foldername serves as a key, i.e. csvFiles[folder] = "CSV files inside folder".', '\n')

    return [folders, csvFiles]


def extractDataDLS(folders, csvFiles):
    # Create dictionary for data of all CSV-files (correlation & intensity data)
    corrData = {}
    intData = {}

    allFiles = []

    # Loop over all folders and over all contained CSV-files
    for folder in folders:
        for csvFile in csvFiles[folder]:
            # Store name of CSV-file (without the .csv ending)
            csvName = pathlib.PurePath(csvFile).name[:-4]
            # Add name of CSV-file to list of all filenames
            allFiles.append(csvName)
            # Read CSV-file
            csvDataFrame = pd.read_csv(csvFile, sep=',')
            # Separate contents into correlation and intensity data
            corrDataFrame = csvDataFrame.iloc[:, 0:4].dropna()
            intDataFrame = csvDataFrame.iloc[:, 4:8].dropna()
            # Save data into dictionaries
            corrData[csvName] = corrDataFrame.values
            intData[csvName] = intDataFrame.values
            
    print('The correlation and intensity data from the CSV-files are stored ' \
        'in the dictionaries "corrData" and "intData"', '\n')
    print('with respective filenames as keys, i.e. corrData[filename] = "correlation data inside given file"!', '\n')

    return [allFiles, corrData, intData]


def extractDataMDI(folders, csvFiles):
    # Create dictionary for data of all CSV-files (correlation & intensity data)
    countData = {}

    allFiles = []

    # Loop over all folders and over all contained CSV-files
    for folder in folders:
        for csvFile in csvFiles[folder]:
            # Store name of CSV-file (without the .csv ending)
            csvName = pathlib.PurePath(csvFile).name[:-4]
            # Add name of CSV-file to list of all filenames
            allFiles.append(csvName)
            # Read CSV-file
            csvDataFrame = pd.read_csv(csvFile, sep=',')
            # Extract MDI Data
            DataFrame = csvDataFrame.iloc[0:, 1:8]
            # Save data into a dictionary (multiply by factor 20 to get # per mL, as 50 uL used)
            countData[csvName] = 20 * DataFrame.values.transpose()
                    
    print('The correlation and intensity data from the CSV-files are stored ' \
        'in the dictionaries "corrData" and "intData"', '\n')
    print('with respective filenames as keys, i.e. corrData[filename] = "correlation data inside given file"!', '\n')

    return [allFiles, countData]


def mergeDataDLS(allFiles, corrData, intData):
    # Sort list of all files
    allFiles = sorted(allFiles)

    # Find total number of files
    numFiles = len(allFiles)

    # Initialize symmetric matrix of same experiments (different name = 0, same name = 1)
    sameExp = np.zeros((numFiles, numFiles))

    collFiles = []

    # Part below still inefficient: would be better to only consider upper triangular matrix in the first place!
    
    for i, firstFile in enumerate(allFiles):
        for j, secondFile in enumerate(allFiles):
        
            # Remove exp. number and date (last 11 characters) from filename
            firstPrefix = firstFile[:-11]
            secondPrefix = secondFile[:-11]
            
            # Initialize dictionaries as zero
            corrData[firstPrefix] = np.zeros((1, 1))
            intData[firstPrefix] = np.zeros((1, 1))
            
            # If prefix matches, set value of sameExp matrix to 1 (same experimental conditions)
            if firstPrefix == secondPrefix:
                sameExp[i,j] = 1

    # Only keep upper triangular matrix (of main diagonal), as matrix is symmetric
    sameExpUpper = np.triu(sameExp, k=0)

    # Set diagonal to zero (trivial)
    np.fill_diagonal(sameExpUpper, 0)
                
    # Find indices of nonzero elements in sameExp matrix, returns tuple of arrays of coordinates
    samesame = np.nonzero(sameExpUpper) 

    # Initialize array of previously used indices
    prevCorrIdx = []
    prevIntIdx = []

    # Create set of same name indices (easier to check whether certain index is present)
    sameSet0 = set(samesame[0])
    sameSet1 = set(samesame[1])
    sameSet = sameSet0 | sameSet1

    for i, firstFile in enumerate(allFiles):
        for j, secondFile in enumerate(allFiles):
            
            # Remove exp. number and date (last 11 characters) from filename
            firstPrefix = firstFile[:-11]
            secondPrefix = secondFile[:-11]
            
            # Create collection files (collecting all data from experiments with the same exp. conditions)
            if (firstPrefix not in collFiles):
                collFiles.append(firstPrefix)
                
            if (i not in sameSet) and (i not in prevCorrIdx) and (j == 0):
                # Add new entries to dictionary (w/o exp. number and date)
                corrData[firstPrefix] = corrData[firstFile]
                intData[firstPrefix] = intData[firstFile]

                # Add index to list of previous indices
                prevCorrIdx.append(i)
                prevIntIdx.append(i)
                
            elif (i in sameSet0) and (j in sameSet1):
                # Find coordinates of index pairs in samesame arrays
                coord = np.argwhere((samesame[0] == i) & (samesame[1] == j))
                
                # If indices match, append 
                if (coord.size > 0):
                    
                    if i not in prevCorrIdx:
                        corrData[firstPrefix] = corrData[firstFile]
                        prevCorrIdx.append(i)
                    
                    if j not in prevCorrIdx:
                        corrData[firstPrefix] = np.append(corrData[firstPrefix], corrData[secondFile][:,1:], axis=1)
                        prevCorrIdx.append(j)
                    
                    if (i not in prevIntIdx) and intData[firstFile].any():
                        #print(i, j, 'firstFileInt')
                        intData[firstPrefix] = intData[firstFile]
                        prevIntIdx.append(i)    
                    
                    if (j not in prevIntIdx) and not intData[firstPrefix].any() and intData[secondFile].any():
                        #print(i, j, 'secondFileInt')
                        intData[firstPrefix] = intData[secondFile]
                        prevIntIdx.append(j)
                        
                    elif (j not in prevIntIdx)and intData[firstPrefix].any() and intData[secondFile].any():
                        #print(i, j, 'firstandsecondFileInt')
                        intData[firstPrefix] = np.append(intData[firstPrefix], intData[secondFile][:,1:], axis=1)
                        prevIntIdx.append(j)
                            

    # Export all files with same experimental conditions as average and as collection CSV-files
    print('The following collections of same experimental conditions are generated (with different suffixes): \n')

    for collFile in collFiles:
        print(collFile)
        
        fileNameCorrAvg = collFile + '_AvgCorr' + '.csv'
        fileNameIntAvg = collFile + '_AvgInt' + '.csv'
        fileNameCorrAll = collFile + '_AllCorr' + '.csv'
        fileNameIntAll = collFile + '_AllInt' + '.csv'
        
        filePathCorrAvg = os.path.join('01_Data', '01_DLS', '02_Average', '01_Correlation', fileNameCorrAvg)
        filePathIntAvg = os.path.join('01_Data', '01_DLS', '02_Average', '02_Intensity', fileNameIntAvg)
        filePathCorrAll = os.path.join('01_Data', '01_DLS', '03_Collection', '01_Correlation', fileNameCorrAll)
        filePathIntAll = os.path.join('01_Data', '01_DLS', '03_Collection', '02_Intensity', fileNameIntAll)
            
        firstColCorr = corrData[collFile][:,0]
        meanCorr = np.mean(corrData[collFile][:, 1:], axis=1)
        stDevCorr = np.std(corrData[collFile][:, 1:], axis=1, ddof=1)
        
        fileContentCorrAvg = np.transpose(np.array([firstColCorr, meanCorr,stDevCorr]))
        fileContentCorrAll = corrData[collFile]
        
        np.savetxt(filePathCorrAvg, fileContentCorrAvg, delimiter=",", header='Lag Time [us], Avg. Corr. [-], Std. Dev. Corr. [-]', comments='')
        np.savetxt(filePathCorrAll, fileContentCorrAll, delimiter=",", header='Lag Time [us], Corr. 1.1 [-], Corr 1.2 [-], Corr. 1.3 [-], Corr. 2.1 [-], Corr. 2.2 [-], Corr. 2.3 [-], ...', comments='')
        
        if intData[collFile].any():
            firstColInt = intData[collFile][:,0]
            meanInt = np.mean(intData[collFile][:, 1:], axis=1)
            stDevInt = np.std(intData[collFile][:, 1:], axis=1)
            
            fileContentIntAvg = np.transpose(np.array([firstColInt, meanInt, stDevInt]))
            fileContentIntAll = intData[collFile]

            np.savetxt(filePathIntAvg, fileContentIntAvg, delimiter=",", header='Hydrodyn. Diameter [nm], Avg. Int. [-], Std. Dev. Int. [-]', comments='')
            np.savetxt(filePathIntAll, fileContentIntAll, delimiter=",", header='Hydrodyn. Diameter [nm], Int. 1.1 [-], Int. 1.2 [-], Int. 1.3 [-], Int. 2.1 [-], Int. 2.2 [-], Int. 2.3 [-], ...')

    return [corrData, intData, collFiles]


def mergeDataMDI(allFiles, countData):
    # Sort list of all files
    allFiles = sorted(allFiles)

    # Find total number of files
    numFiles = len(allFiles)

    # Initialize symmetric matrix of same experiments (different name = 0, same name = 1)
    sameExp = np.zeros((numFiles, numFiles))

    #print(allFiles)

    collFiles = []

    # Part below still inefficient: would be better to only consider upper triangular matrix in the first place!
    
    for i, firstFile in enumerate(allFiles):
        for j, secondFile in enumerate(allFiles):
        
            # Remove exp. number and date (last 11 characters) from filename
            firstPrefix = firstFile[:-11]
            secondPrefix = secondFile[:-11]
            
            # Initialize dictionaries as zero
            countData[firstPrefix] = np.zeros((1, 1))
            
            # If prefix matches, set value of sameExp matrix to 1 (same experimental conditions)
            if firstPrefix == secondPrefix:
                sameExp[i,j] = 1

    # Only keep upper triangular matrix (of main diagonal), as matrix is symmetric
    sameExpUpper = np.triu(sameExp, k=0)

    # Set diagonal to zero (trivial)
    np.fill_diagonal(sameExpUpper, 0)
                
    # Find indices of nonzero elements in sameExp matrix, returns tuple of arrays of coordinates
    samesame = np.nonzero(sameExpUpper) 

    # Initialize array of previously used indices
    prevCountIdx = []

    # Create set of same name indices (easier to check whether certain index is present)
    sameSet0 = set(samesame[0])
    sameSet1 = set(samesame[1])
    sameSet = sameSet0 | sameSet1

    for i, firstFile in enumerate(allFiles):
        for j, secondFile in enumerate(allFiles):
            
            # Remove exp. number and date (last 11 characters) from filename
            firstPrefix = firstFile[:-11]
            secondPrefix = secondFile[:-11]
            
            # Create collection files (collecting all data from experiments with the same exp. conditions)
            if (firstPrefix not in collFiles):
                collFiles.append(firstPrefix)
                
            if (i not in sameSet) and (i not in prevCountIdx) and (j == 0):
                # Add new entries to dictionary (w/o exp. number and date)
                countData[firstPrefix] = countData[firstFile]

                # Add index to list of previous indices
                prevCountIdx.append(i)
                
            elif (i in sameSet0) and (j in sameSet1):
                # Find coordinates of index pairs in samesame arrays
                coord = np.argwhere((samesame[0] == i) & (samesame[1] == j))
                
                # If indices match, append 
                if (coord.size > 0):
                    
                    if i not in prevCountIdx:
                        countData[firstPrefix] = countData[firstFile]
                        prevCountIdx.append(i)
                                    
                    if j not in prevCountIdx:
                        countData[firstPrefix] = np.append(countData[firstPrefix], countData[secondFile], axis=1)
                        prevCountIdx.append(j)
                    

    # Export all files with same experimental conditions as average and as collection CSV-files
    print('The following collections of same experimental conditions are generated (with different suffixes): \n')

    for collFile in collFiles:
        
        print(collFile)

        fileNameCountAvg = collFile + '_AvgCount' + '.csv'
        fileNameCountAll = collFile + '_AllCount' + '.csv'
        
        filePathCountAvg = os.path.join('01_Data', '02_MDI', '02_Average', fileNameCountAvg)
        filePathCountAll = os.path.join('01_Data', '02_MDI', '03_Collection', fileNameCountAll)
            
        if countData[collFile].shape[1] > 1:
            meanCount = np.mean(countData[collFile], axis=1)
            stdDevCount = np.std(countData[collFile], axis=1, ddof=1)
        else: 
            meanCount = countData[collFile][:, 0]
            stdDevCount = np.zeros(countData[collFile].shape[0])
        
        fileContentCountAvg = np.transpose(np.array([meanCount,stdDevCount]))
        fileContentCountAll = countData[collFile]
        
        np.savetxt(filePathCountAvg, fileContentCountAvg, delimiter=",", header='Avg. Num. Dens. [1/mL], Std. Dev. Num. Dens. [1/mL]', footer='Rows: ECD 1-2; ECD 2-5; ECD 5-10, ECD 10-25; ECD 25-70; ECD ge70; MFD ge70', comments='')
        np.savetxt(filePathCountAll, fileContentCountAll, delimiter=",", header='Num. Dens. 1.1 [1/mL], Num. Dens. 1.2 [1/mL], Num. Dens. 2.1 [1/mL], Num. Dens. 2.2 [1/mL], ...', footer='Rows: ECD 1-2; ECD 2-5; ECD 5-10, ECD 10-25; ECD 25-70; MFD ge70', comments='')

    return [countData, collFiles]
