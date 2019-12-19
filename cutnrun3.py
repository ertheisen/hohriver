#fastqc_pre.py
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
		sys.exit()

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

print('Runtime pre-processing fastq (hh:mm:ss): ' + str(datetime.now() - startTime))

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


startTime = datetime.now()

# Change directory within container
os.chdir('/data/trim_out')
print('Current working directory is:' +  os.getcwd())

#prepare input files
input_files = sorted(glob.glob('*val*'))
input_files_R1 = []
input_files_R2 = []

if len(input_files) == 0:
	print('No quality-controlled input files from trim_galore, checking input folder for fastqc output...')
	print('\n')
	input_files = sorted(glob.glob('*.fastq'))
	if len(input_files) == 0:
		print('No uncompressed fastq files detected, looking for fastq.gz...')
		print('\n')
		input_files = sorted(glob.glob('*.fastq.gz*'))
		if len(input_files) == 0:
			print('No valid input files detected. Exiting...')	
			sys.exit()
		else:
			print('Decompressing .gz files')
			#create fastq filename from .fastq.gz
			temp_list = [f.replace('.gz', '') for f in input_files]
			#run system command to decompress file with zcat
			for i in range(len(input_files)):
				temp_str= 'zcat ' + input_files[i] + ' > ' + temp_list[i]
				check_output(temp_str, shell=True)
			input_files = temp_list

for item in input_files:
	if '_R1' in item and '_R2' in item:
		print("Input file: " + item + " contains both strings 'R1' and 'R2'. Not including...")
		sys.exit()
	elif '_R1' in item and '_R2' not in item:
		input_files_R1.append(item)
	elif '_R1' not in item and '_R2' in item:
		input_files_R2.append(item)
	else:
		print("Input file: " + item + "does not contain string 'R1' or 'R2'. Not including...")

if len(input_files_R1) != len(input_files_R2):
	print('Unequal numbers of files assigned as R1 and R2. Check naming convention. Exiting...')
	sys.exit()

if (len(input_files_R1) + len(input_files_R2)) != len(input_files):
	print('Not all of input files assigned as R1 or R2. Exiting...')
	sys.exit()



print('Files assigned as R1:')
print('\n'.join(input_files_R1))
print('\n')

print('Files assigned as R2:')
print('\n'.join(input_files_R2))
print('\n')

#detect spike-in genome preferences
fly_spike = os.path.isfile('/genomes/spike/dm6/Sequence/Bowtie2Index/genome.1.bt2')

yeast_spike = os.path.isfile('/genomes/spike/sacCer3/Sequence/Bowtie2Index/genome.1.bt2')

ecoli_spike = os.path.isfile('/genomes/spike/2008-03-17/Sequence/Bowtie2Index/genome.1.bt2')

#make sam output names
sam_names = []
spike_sam_names = []

for i in range(len(input_files_R1)):
	sam_name = input_files_R1[i].split('_L0')[0] + '.sam'
	sam_names.append(sam_name)

print('Output sam filenames for data species:')
print('\n'.join(sam_names))
print('\n')

if fly_spike == True and yeast_spike+ecoli_spike == 0:
	spike_species = 'dm6'
	print('Spike-in chromatin from Drosophila')
	for i in range(len(sam_names)):
		temp_name = sam_names[i] + '.dm6'
		spike_sam_names.append(temp_name)
	print('Output sam filenames for spike species:')
	print('\n'.join(spike_sam_names))
	print('\n')
	print('Data will be aligned to hg19 and ' + spike_species)

elif yeast_spike == True and fly_spike+ecoli_spike == 0:
	spike_species = 'sacCer3'
	print('Spike-in chromatin from Saccharomyces')
	for i in range(len(sam_names)):
		temp_name = sam_names[i] + '.sacCer3'
		spike_sam_names.append(temp_name)
	print('Output sam filenames for spike species:')
	print('\n'.join(spike_sam_names))
	print('\n')
	print('Data will be aligned to hg19 and ' + spike_species)

elif ecoli_spike == True and fly_spike+yeast_spike == 0:
	spike_species = 'ecoli'
	print('Spike-in chromatin from E. coli')
	for i in range(len(sam_names)):
		temp_name = sam_names[i] + '.ecoli'
		spike_sam_names.append(temp_name)
	print('Output sam filenames for spike species:')
	print('\n'.join(spike_sam_names))
	print('\n')
	print('Data will be aligned to hg19 and ' + spike_species)

else:
	spike_species = 'null'
	print('Data will be aligned to hg19')
	print('No spike-in alignment will be performed.')
	print('\n')
	print('\n')

#define bowtie2 index for hg19
hg19_index = '/genomes/test/Sequence/Bowtie2Index/genome'
print('Bowtie index for data files found at:')
print(hg19_index)
print('\n')

#run bowtie for hg19
print('Running bowtie2 alignment for hg19')
print('\n')
for i in range(len(input_files_R1)):
	print('count = ' +str(i))
	print('\n')
	#create string for system command
	temp_str = 'bowtie2 --no-unal --no-mixed --no-discordant --dovetail --phred33 -q -I 10 -X 700 --threads 16' \
	+ ' -x ' + hg19_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' >> ' + sam_names[i]

	print(temp_str)

	check_output(temp_str, shell=True)
	print('\n')

#define bowtie2 index for spike_species
if spike_species == 'dm6': 
	spike_index = '/genomes/spike/dm6/Sequence/Bowtie2Index/genome'
	print('Bowtie index for spike files found at:')
	print(spike_index)
	print('\n')

elif spike_species == 'sacCer3':
	spike_index = '/genomes/spike/sacCer3/Sequence/Bowtie2Index/genome'
	print('Bowtie index for spike files found at:')
	print(spike_index)
	print('\n')

elif spike_species == 'ecoli':
	spike_index = '/genomes/spike/2008-03-17/Sequence/Bowtie2Index/genome'
	print('Bowtie index for spike files found at:')
	print(spike_index)
	print('\n')

elif spike_species not in ['dm6', 'sacCer3', 'ecoli']:
	print('Bowtie index unavailable for spike species or no spike specied index detected.')
	print('\n')

#run bowtie for spike_species
if spike_species in ['dm6', 'sacCer3', 'ecoli']:
	print('Running bowtie2 alignment for ' + spike_species)
	print('\n')
	for i in range(len(input_files_R1)):
		print('count = ' +str(i))
		print('\n')
		#create string for system command
		temp_str = 'bowtie2 --very-sensitive --no-unal --no-mixed --no-discordant --dovetail --phred33 -q -I 10 -X 700 --threads 16' \
		+ ' -x ' + spike_index + ' -1 ' + input_files_R1[i] + ' -2 ' + input_files_R2[i] + ' >> ' + spike_sam_names[i]

		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')

#make new directory for output
os.chdir('/data')
print('Current working directory is:' +  os.getcwd())
os.mkdir('sams')
os.chdir('/data/trim_out')

#copy files to output folder
output_dir = '/data/sams/'
print('Moving sam files to output folder')
print('\n')
for i in range(len(sam_names)):
	shutil.move(sam_names[i], output_dir)

if spike_species in ['dm6', 'sacCer3', 'ecoli']:
	for i in range(len(spike_sam_names)):
		shutil.move(spike_sam_names[i], output_dir)


print('Alignment Runtime (hh:mm:ss): ' + str(datetime.now() - startTime))
print('\n')

###SAM conversion to bam, bedgraph, and BigWig

os.mkdir('/data/bams')
os.chdir('/data/sams')
print('Current working directory is:' +  os.getcwd())
print('\n')

import pybedtools
from pybedtools import BedTool
import pandas as pd
from pybedtools.helpers import chromsizes
from pybedtools.contrib.bigwig import bedgraph_to_bigwig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

startTime = datetime.now()

#######################################################
## spike in or average normalization with file generation#########
#######################################################

if spike_species in ['dm6', 'sacCer3', 'ecoli']:	
	print('Converting sams to bams, bedgraph, and BigWig with spike-in normalization.')
	print('\n')
	datafiles = sorted(glob.glob('*.sam'))
	spikefiles = sorted(glob.glob('*.' + spike_species))
	
	print('\n')
	print('Data files loaded:')
	print('\n'.join(datafiles))
	print('\n')
	print('Spike files loaded:')
	print('\n'.join(spikefiles))
	print('\n')

	##check equal spike and data files
	if len(datafiles) == 0:
		print('No data files imported. Exiting...')
		sys.exit()
	
	if len(datafiles) != len(spikefiles):
		print('Unequal number of data and spike files. Exiting...')
		sys.exit()

	for i in range(len(datafiles)):
		if datafiles[i] not in spikefiles[i]:
			print('Unmatched pairs and spike files.')
			print(datafiles[i] + ' not paired with ' + spikefiles[i])
			print('Exiting...')
			sys.exit()

else:
	print('Converting sams to bams, bedgraph, and BigWig without spike-in normalization.')
	print('\n')
	datafiles = sorted(glob.glob('*.sam'))

	print('\n')
	print('Data files loaded:')
	print('\n'.join(datafiles))
	print('\n')
	print('No spike-in normalization')
	print('\n')

##convert to bam format
print('Converting to bam format')
print('\n')

bam_names = [f.replace('sam', 'bam') for f in datafiles]

print('\n')
print('\n'.join(bam_names))
print('\n')
print('\n')
print('SAM to BAM')
print('\n')

bam_string = []

for i in range(len(bam_names)):
        bam_string.append('samtools view -b -S ' + datafiles[i] + ' > /data/bams/' + bam_names[i])

for item in bam_string:
        check_output(item, shell = True)

datafiles = bam_names

##Sort and index
print('\n')
print('Sorting bams')
print('\n')

sorted_bam_names = []
for i in range(len(bam_names)):
	sorted_bam_name = 'sorted.' + bam_names[i]
	sorted_bam_names.append(sorted_bam_name)

bam_names_sorted = [f.replace('.bam', '') for f in sorted_bam_names]

sort_string = []

for i in range(len(bam_names)):
		sort_string.append('samtools sort /data/bams/' + bam_names[i] + ' /data/bams/' + bam_names_sorted[i])

for item in sort_string:
		check_output(item, shell = True)

print('\n')
print('Bam files to index:')
print('\n'.join(bam_names_sorted))
print('\n')


print('Indexing bams')
print('\n')             
os.chdir('/data/bams')
print('Current working directory is:' +  os.getcwd())

index_string = []

for i in range(len(bam_names_sorted)):
        index_string.append('samtools index ' + sorted_bam_names[i])

for item in index_string:
        check_output(item, shell = True)

sorted_files = sorted(glob.glob('sorted*'))
os.mkdir('/data/bams/sorted_bams')
output_dir = '/data/bams/sorted_bams'
for i in range(len(sorted_files)):
	shutil.move(sorted_files[i], output_dir)

##generate bed files from bam files
print('Generating bed files representing whole insert from paired end reads in the data files')
print('\n')
print('Also generating size-selected bed files for inserts <= 140bp, those 150-300 bp, and those > 300 bp')


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


if spike_species in ['dm6', 'sacCer3', 'ecoli']:
	print('Generating spike-normalized bedgraphs for all the bed files')
	print('\n')
	#calculate number of reads in spike files
	#spikecount is a list of reads populated by the 'samtools view -c' shell command for each spikefile	
	spikecount = []
	spike_string = []

	for item in spikefiles:
		spike_string.append('samtools view -c -S /data/sams/' + item)

	for item in spike_string:
		spikecount.append((int(check_output(item, shell = True)))/2)

	print('spike counts')
	print(spikecount)
	print('\n')

	#calculate scaling factors
	scaling_factor = []
	for item in spikecount:
		scaling_factor.append(float(10000)/item)

	for i in range(len(spikefiles)):
		count = str(spikecount[i])
		print('\n')
		print('for ' + spikefiles[i] + ' count:') 
		print(count)
		print('\n')
		print('scaling factor for ' + bed_names[i] + ' is:')
		print(scaling_factor[i])

else:
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
print('Runtime file conversion (hh:mm:ss): ' + str(datetime.now() - startTime))
#mapqc_fingerprint.py

startTime = datetime.now()

os.chdir('/data/bams')
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
	print('Bam files must be indexed. Sorting and indexing...')
	for i in range(len(bam_files)):
		print('\n')
		#create string for system command to sort
		temp_str = 'samtools sort ' + bam_files[i] + ' ' + bam_files[i].split('.')[0] + '_sorted'
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')

	sorted_files = sorted(glob.glob('*sorted*'))
	os.mkdir('/data/unsorted_bams')
	print('Moving unsorted bams to new directory')
	output_dir = '/data/unsorted_bams'
	for i in range(len(bam_files)):
		shutil.move(bam_files[i], output_dir)
	
	for i in range(len(sorted_files)):
		print('\n')	
		#create string for system command to index
		temp_str = 'samtools index ' + sorted_files[i]
		print(temp_str)

		check_output(temp_str, shell=True)
		print('\n')
	
#run plotFingerprint
startTime = datetime.now()

os.chdir('/data/bams/sorted_bams')
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
