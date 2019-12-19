#trim_qc.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output

startTime = datetime.now()

# Change directory within container
os.chdir('/data')
print('Current working directory is:' +  os.getcwd())
print('\n')
os.mkdir('trim_out')
os.mkdir('trim_out/fastqc_post_trim')
os.chdir('/data/fastq')
print('Current working directory is:' +  os.getcwd())
print('\n')

#Determine if fastq files are present
input_fastq = sorted(glob.glob('*.fastq'))
input_fastq_R1 = sorted(glob.glob('*_R1*.fastq'))
input_fastq_R2 = sorted(glob.glob('*_R2*.fastq'))

if len(input_fastq)==0:
	#import input filenames as a list and print out files
	input_files_R1 = sorted(glob.glob('*_R1*.fastq.gz'))
	if len(input_files)==0:
		print('No R1 files imported. Exiting...')
		sys.exit()

	print('Input files R1:')
	print('\n'.join(input_files_R1))
	print('\n')

	print('Decompressing .gz files')
	#create fastq filename from .fastq.gz
	temp_list_R1 = [f.replace('.gz', '') for f in input_files_R1]
	#run system command to decompress file with zcat
	for i in range(len(input_files_R1)):
		temp_str_R1 = 'zcat ' + input_files_R1[i] + ' > ' + temp_list_R1[i]
		check_output(temp_str_R1, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files_R1 = temp_list_R1

	#import input filenames as a list and print out files
	input_files_R2 = sorted(glob.glob('*_R2*.fastq.gz'))
	if len(input_files)==0:
		print('No R2 files imported. Exiting...')
		sys.exit()

	print('Input files R2:')
	print('\n'.join(input_files_R2))
	print('\n')

	print('Decompressing .gz files')
	#create fastq filename from .fastq.gz
	temp_list_R2 = [f.replace('.gz', '') for f in input_files_R2]
	#run system command to decompress file with zcat
	for i in range(len(input_files_R2)):
		temp_str_R2 = 'zcat ' + input_files_R2[i] + ' > ' + temp_list_R2[i]
		check_output(temp_str_R2, shell=True)

	#replace list of input_files with temp_list so future operations are on decompressed files
	input_files_R2 = temp_list_R2

else:
	input_files_R1 = input_fastq_R1
	input_files_R2 = input_fastq_R2

if len(input_files_R1)!=len(input_files_R2):
	print('Unpaired files detected. ...Exiting.')
	sys.exit()

print('Decompressed files to be analyzed:')
print('\n'.join(input_files_R1))
print('\n'.join(input_files_R2))
print('\n')

for i in range(len(input_files_R1)):
#create string for system command
	temp_str = 'trim_galore -o /data/trim_out --paired ' + input_files_R1[i] + ' ' + input_files_R2[i]

	print(temp_str)

	check_output(temp_str, shell=True)
	print('\n') 

os.chdir('/data/trim_out')
print('Current working directory is:' +  os.getcwd())
fastqc_input = sorted(glob.glob('*.fq'))

print('Trimmed files to be analyzed:')
print('\n'.join(fastqc_input))
print('\n')

for i in range(len(fastqc_input)):
#create string for system command
	temp_str = 'fastqc --outdir=fastqc_post_trim ' + fastqc_input[i]

	print(temp_str)

	check_output(temp_str, shell=True)
	print('\n') 


print('Runtime trimming (hh:mm:ss): ' + str(datetime.now() - startTime))
