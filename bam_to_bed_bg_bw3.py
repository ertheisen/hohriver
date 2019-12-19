#bamtobed_cutnrun.py
from datetime import datetime
import os
import glob
import shutil
import sys
from subprocess import check_output
import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


startTime = datetime.now()
#######################################
# Change directory within container
os.chdir('/data')
print('Current working directory is:' +  os.getcwd())
#######################################
bam_files = sorted(glob.glob('*bam'))
if bam_files == 0:
	# Change directory within container
	os.chdir('/data/bams')
	print('Current working directory is:' +  os.getcwd())
	bam_files = sorted(glob.glob('*bam'))
	if bam_files == 0:
		print('No bam files detected. Exiting...')
		sys.exit()
	else:
		work_dir = 2
else:
	work_dir = 1

bai_files = sorted(glob.glob('*bai'))
if bai_files == 0:
	print('Bam files will be indexed. Sorting and indexing...')
	for i in range(len(bam_files)):
		print('\n')
		#create string for system command to sort
		temp_str = 'samtools sort ' + bam_files[i] + ' sorted.' + bam_files[i].split('.')[0]
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')

	sorted_files = sorted(glob.glob('sorted*'))
	os.mkdir('/data/sorted_bams')
	print('Moving sorted bams to new directory')
	output_dir = '/data/sorted_bams'
	for i in range(len(sorted_files)):
		shutil.move(sorted_files[i], output_dir)
	os.chdir('/data/sorted_bams')
	for i in range(len(sorted_files)):
		print('\n')	
		#create string for system command to index
		temp_str = 'samtools index /data/sorted_bams' + sorted_files[i]
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')

##############################################

###BAM conversion to BED, bedgraph, and BigWig

print('Converting bams to beds, bedgraph, and BigWig without spike-in normalization.')
print('\n')

datafiles = bam_files
##generate bed files from bam files
print('Generating bed files representing whole insert from paired end reads in the data files')
print('\n')

if work_dir == 1:
	os.chdir('/data')

if work_dir == 2:
	os.chdir('/data/bams')

print('Current working directory is:' +  os.getcwd())
print('\n')

#generate bed file names
bed_names = [f.replace('bam', 'bed') for f in datafiles]
size_selected_small = [f.replace('bam', 'small.bed') for f in datafiles]
size_selected_med = [f.replace('bam', 'med.bed') for f in datafiles]
size_selected_big = [f.replace('bam', 'big.bed') for f in datafiles]

#generate file names for length analysis
lengths_names = [f.replace('bam', 'lengths') for f in datafiles]

#generate bed files with bam_to_bed tool (makes bed12 format)
for i in range(len(datafiles)):
	temp_bed = BedTool(datafiles[i]).bam_to_bed(bedpe=True).to_dataframe()

	#need to strip out start and end position of whole insert (bed12 is both reads)
	#column names actually represent <chrom>, <start of insert>, <end of insert>
	temp_bed_stripped = temp_bed.iloc[:,[0,1,5]].sort_values(by = ['chrom', 'start', 'strand'])

	#calculate insert size as column 4 and save file with bed_name
	temp_bed_stripped['length'] = temp_bed_stripped['strand'] - temp_bed_stripped['start']

	temp_bed_stripped.to_csv(bed_names[i], sep="\t", header = False, index = False)

	#analyze lengths of inserts
	temp_lengths = temp_bed_stripped.groupby(by=['length'])['length'].count()

	temp_lengths.to_csv(lengths_names[i], sep="\t", header = [bed_names[i]], index = True, index_label='length')

	#generate size-selected whole insert beds
	subset_small = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=20) & (temp_bed_stripped.iloc[:,3]<=140)]
	subset_small.to_csv(size_selected_small[i], sep="\t", header = False, index = False)

	subset_med = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=150) & (temp_bed_stripped.iloc[:,3]<=300)]
	subset_med.to_csv(size_selected_med[i], sep="\t", header = False, index = False)

	subset_big = temp_bed_stripped[(temp_bed_stripped.iloc[:,3]>=301) & (temp_bed_stripped.iloc[:,3]<=700)]
	subset_big.to_csv(size_selected_big[i], sep="\t", header = False, index = False)



print('Finished generating bed files:')
print('\n')
print('whole insert bed files:' + '\n' + '\n'.join(bed_names))
print('\n')
print('bed files for inserts < 140 bp:' + '\n' + '\n'.join(size_selected_small))
print('\n')
print('bed files for inserts 150 bp - 300 bp:' + '\n' + '\n'.join(size_selected_med))
print('\n')
print('bed files for inserts > 300 bp:' + '\n' + '\n'.join(size_selected_big))


#generate normalized bedgraphs
print('Current working directory is:' +  os.getcwd())
print('\n')

#generate bedgraph names
bg_names = [f.replace('bed', 'bg') for f in bed_names]
size_selected_small_bg = [f.replace('bed', 'bg') for f in size_selected_small]
size_selected_med_bg = [f.replace('bed', 'bg') for f in size_selected_med]
size_selected_big_bg = [f.replace('bed', 'bg') for f in size_selected_big]

print('Generating average normalized bedgraphs for all the bed files')
print('\n')
#count total number of reads in each bed file (before size selection)
read_count = []
for item in bed_names:
	read_count.append(BedTool(item).count())

print(read_count)

#calculate genome size
genome_file = chromsizes('hg19')
DF = pd.DataFrame.from_dict(genome_file, orient='index')
genome_size = DF[1].sum()
	
scaling_factor = []
for i in range(len(read_count)):
	scaling_factor.append(float(genome_size) / read_count[i])

for i in range(len(bed_names)):
	count = str(read_count[i])	
	print('\n')
	print('for ' + bed_names[i] + ' read count:')
	print(count)
	print('\n')
	print('scaling factor for ' + bed_names[i] + ' is:')
	print(scaling_factor[i])


#run bedtools genomecov to generate bedgraph files
for i in range(len(bg_names)):
	BedTool(bed_names[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(bg_names[i])

for i in range(len(size_selected_small_bg)):
	BedTool(size_selected_small[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(size_selected_small_bg[i])

for i in range(len(size_selected_med_bg)):
	BedTool(size_selected_med[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(size_selected_med_bg[i])

for i in range(len(size_selected_big_bg)):
	BedTool(size_selected_big[i]).genome_coverage(bg = True, genome = 'hg19', scale = scaling_factor[i]).moveto(size_selected_big_bg[i])

print('Finished generating bedgraph files:')
print('\n')
print('whole insert bedgraph files:' + '\n' + '\n'.join(bg_names))
print('\n')
print('bedgraph files for inserts < 140 bp:' + '\n' + '\n'.join(size_selected_small_bg))
print('\n')
print('bedgraph files for inserts 150 bp - 300 bp:' + '\n' + '\n'.join(size_selected_med_bg))
print('\n')
print('bedgraph files for inserts > 300 bp:' + '\n' + '\n'.join(size_selected_big_bg))

##make bigwig files
print('Current working directory is:' +  os.getcwd())
print('\n')

print('Generating big_wig files for all the bedgraphs')
print('\n')

#generate bigwig names
bw_names = [f.replace('bg', 'bw') for f in bg_names]
size_selected_small_bw = [f.replace('bg', 'bw') for f in size_selected_small_bg]
size_selected_med_bw = [f.replace('bg', 'bw') for f in size_selected_med_bg]
size_selected_big_bw = [f.replace('bg', 'bw') for f in size_selected_big_bg]

#run bedgraph_to_bigwig tool
for i in range(len(bg_names)):
	bedgraph_to_bigwig(BedTool(bg_names[i]), 'hg19', bw_names[i])
for i in range(len(size_selected_small_bg)):
	bedgraph_to_bigwig(BedTool(size_selected_small_bg[i]), 'hg19', size_selected_small_bw[i])
for i in range(len(size_selected_med_bg)):
	bedgraph_to_bigwig(BedTool(size_selected_med_bg[i]), 'hg19', size_selected_med_bw[i])
for i in range(len(size_selected_big_bg)):
	bedgraph_to_bigwig(BedTool(size_selected_big_bg[i]), 'hg19', size_selected_big_bw[i])


print('Finished generating bigwig files:')
print('\n')
print('whole insert bigwig files:' + '\n' + '\n'.join(bw_names))
print('\n')
print('bigwig files for inserts < 140 bp:' + '\n' + '\n'.join(size_selected_small_bw))
print('\n')
print('bigwig files for inserts 150 bp - 300 bp:' + '\n' + '\n'.join(size_selected_med_bw))
print('\n')
print('bigwig files for inserts > 300 bp:' + '\n' + '\n'.join(size_selected_big_bw))

os.mkdir('/data/beds')
os.mkdir('/data/beds/size_selected_small')
os.mkdir('/data/beds/size_selected_med')
os.mkdir('/data/beds/size_selected_big')
os.mkdir('/data/beds/lengths')
os.mkdir('/data/bedgraphs')
os.mkdir('/data/bedgraphs/size_selected_small')
os.mkdir('/data/bedgraphs/size_selected_med')
os.mkdir('/data/bedgraphs/size_selected_big')
os.mkdir('/data/bigwigs')
os.mkdir('/data/bigwigs/size_selected_small')
os.mkdir('/data/bigwigs/size_selected_med')
os.mkdir('/data/bigwigs/size_selected_big')

output_dir0 = '/data/beds'
output_dir1 = '/data/beds/lengths'
output_dir2 = '/data/beds/size_selected_small'
output_dir3 = '/data/beds/size_selected_med'
output_dir4 = '/data/beds/size_selected_big'
for i in range(len(bed_names)):
	shutil.move(bed_names[i], output_dir0)
for i in range(len(lengths_names)):
	shutil.move(lengths_names[i], output_dir1)
for i in range(len(size_selected_small)):		
	shutil.move(size_selected_small[i], output_dir2)
for i in range(len(size_selected_med)):	
	shutil.move(size_selected_med[i], output_dir3)
for i in range(len(size_selected_big)):	
	shutil.move(size_selected_big[i], output_dir4)

print('Moving bedgraphs to output folder')
print('\n')
output_dir5 = '/data/bedgraphs'
output_dir6 = '/data/bedgraphs/size_selected_small'
output_dir7 = '/data/bedgraphs/size_selected_med'
output_dir8 = '/data/bedgraphs/size_selected_big'
for i in range(len(bg_names)):
	shutil.move(bg_names[i], output_dir5)
for i in range(len(size_selected_small_bg)):		
	shutil.move(size_selected_small_bg[i], output_dir6)
for i in range(len(size_selected_med_bg)):	
	shutil.move(size_selected_med_bg[i], output_dir7)
for i in range(len(size_selected_big_bg)):	
	shutil.move(size_selected_big_bg[i], output_dir8)

print('Moving bigwigs to output folder')
print('\n')
output_dir9 = '/data/bigwigs'
output_dir10 = '/data/bigwigs/size_selected_small'
output_dir11 = '/data/bigwigs/size_selected_med'
output_dir12 = '/data/bigwigs/size_selected_big'
for i in range(len(bw_names)):
	shutil.move(bw_names[i], output_dir9)
for i in range(len(size_selected_small_bw)):		
	shutil.move(size_selected_small_bw[i], output_dir10)
for i in range(len(size_selected_med_bw)):	
	shutil.move(size_selected_med_bw[i], output_dir11)
for i in range(len(size_selected_big_bw)):	
	shutil.move(size_selected_big_bw[i], output_dir12)


print('Finished')
print('\n')
print('Runtime (hh:mm:ss): ' + str(datetime.now() - startTime))

#run plotFingerprint

startTime = datetime.now()

os.chdir('/data/sorted_bams')
print('Current working directory is:' +  os.getcwd())

print('Running plotFingerprint for all bam files in folder')
print('\n')
temp_str = 'plotFingerprint -b *.bam --smartLabels --minFragmentLength 20 -T "Mapping Fingerprints" --plotFile mapqc_fingerprint.png --outQualityMetrics mapqc.tab'

print(temp_str)

check_output(temp_str, shell=True)
print('\n')

print('Finished')
print('\n')
print('Runtime mapping QC (hh:mm:ss): ' + str(datetime.now() - startTime))
