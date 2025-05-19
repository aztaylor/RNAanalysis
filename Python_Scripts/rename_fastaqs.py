import os

cwd = os.getcwd()

data_fp = os.path.join(cwd, 'FASTQ_Generation_2023-09-06_06_48_25Z-691456767')

run = False

if run:
    for folder in os.listdir(data_fp):
        if folder == '.DS_Store':
            continue
        for file in os.listdir(os.path.join(data_fp, folder)):
            split1 = file.split('_S')
            split2 = split1[1].split('_L')
            os.rename(os.path.join(data_fp,folder,file), os.path.join(data_fp, folder,split1[0].replace('-','_')+'_L'+split2[1]))