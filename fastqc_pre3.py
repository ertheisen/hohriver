#fastqc_pre.py
from datetime import datetime
import os
import glob
import shutil
from subprocess import check_output
    
startTime = datetime.now()

# Change directory within container
os.chdir('/data')
print('Current working directory is:' +  os.getcwd())
os.mkdir('fastqc_pre_out')
os.chdir('/data/fastq')
print('Current working directory is:' +  os.getcwd())

#Determine if fastq files are present
input_fastq = sorted(glob.glob('*.fastq'))

if len(input_fastq)==0:
	#import input filenames as a list and print out files
	input_files = sorted(glob.glob('*.fastq.gz'))
	if len(input_files)==0:
		print('No files imported. Exiting...')

	print('Input files:')
	print('\n'.join(input_files))
	print('\n')

	print('Decompressing .gz files')
	#create fastq filename from .fastq.gz
	temp_list = [f.replace('.gz', '') for f in input_files]
	#run system command to decompress file with zcat
	for i in range(len(input_files)):
		temp_str = 'zcat ' + input_files[i] + ' > ' + temp_list[i]
		check_output(temp_str, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files = temp_list

else:
	input_files = input_fastq

print('Decompressed files to be analyzed:')
print('\n'.join(input_files))
print('\n')


for i in range(len(input_files)):
	#create string for system command
	temp_str = 'fastqc --outdir=/data/fastqc_pre_out ' + input_files[i]

	print(temp_str)

	check_output(temp_str, shell=True)
	print('\n') 

print('Runtime (hh:mm:ss): ' + str(datetime.now() - startTime))

