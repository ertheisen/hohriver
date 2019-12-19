#mapqc_fingerprint.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output

startTime = datetime.now()

os.chdir('/data')
print('Current working directory is:' +  os.getcwd())
bam_files = sorted(glob.glob('*bam'))
if bam_files == 0:
	# Change directory within container
	os.chdir('/data/bams')
	print('Current working directory is:' +  os.getcwd())
	bam_files = sorted(glob.glob('*bam'))
	if bam_files == 0:
		print('No bam files detected. Exiting...')
		sys.exit()

bai_files = sorted(glob.glob('*bai'))
if bai_files == 0:
	print('Bam files must be indexed. Sorting and indexing... Mapping QC output will be in the subdirectory sorted_bams')
	for i in range(len(bam_files)):
		print('\n')
		#create string for system command to sort
		temp_str = 'samtools sort ' + bam_files[i] + ' sorted.' + bam_files[i].split('.')[0]
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')

	sorted_files = sorted(glob.glob('sorted*'))
	os.mkdir('/data/bams/sorted_bams')
	print('Moving sorted bams to new directory')
	output_dir = '/data/bams/sorted_bams'
	for i in range(len(sorted_files)):
		shutil.move(sorted_files[i], output_dir)
	os.chdir('/data/bams/sorted_bams')
	for i in range(len(sorted_files)):
		print('\n')	
		#create string for system command to index
		temp_str = 'samtools index ' + sorted_files[i]
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')
	
#run plotFingerprint
print('Running plotFingerprint for all bam files in folder')
print('\n')
temp_str = 'plotFingerprint -b *.bam --smartLabels --minFragmentLength 20 -T "Mapping Fingerprints" --plotFile mapqc_fingerprint.png --outQualityMetrics mapqc.tab'

print(temp_str)

check_output(temp_str, shell=True)
print('\n')

print('Finished')
print('\n')
print('Runtime mapping QC (hh:mm:ss): ' + str(datetime.now() - startTime))
