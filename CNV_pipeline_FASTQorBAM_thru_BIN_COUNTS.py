#Main pipeline script, 2019-09-06
#Contributors: Michael Lindberg, Alex Koeppel, Will Chronister, Inusah Diallo

#Import argparse and sys modules.
import argparse
import sys
import re
import subprocess
from subprocess import *
import os.path

#SETS WORKING PATH (FOR PARALLEL VERSION ONLY)
#workingpath="/net/nas-storage/vol97/mcconnell_lab/cnvpipe/"
#workingpath="/net/midtier8/vol97/mcconnell_lab/cnvpipe/"
workingpath="/nv/vol97/mcconnell_lab/cnvpipe/"

out = subprocess.Popen("echo ~",stdout=subprocess.PIPE, stderr=None, shell=True)
lookout = out.communicate()[0].strip()
m=re.match('/.*/(.*)',lookout)
username=m.group(1)
#Function to locate the directories containing required programs and die if not found.
def find_software(softname):
        findcomm = "find ~/* -type f -name " + softname
        softlook= subprocess.Popen(findcomm, stdout=subprocess.PIPE, stderr=None, shell=True)
        lookout = softlook.communicate()
        alldir = lookout[0].split('\n')
        ct=0
        while True:
                try:
                        sts = call(alldir[ct],shell=True)
                except Exception:
                        ct=ct+1
                        continue
                softdir=alldir[ct]
                break
        if alldir =='':
                sys.exit("No path found for " + softname + ".  If you know that " + softname + " is properly installed, use the appropriate flag to point to the the directory which contains " + softname + ".")
        thepatt=re.match('(.+/)(.*)',softdir)
        softdir=thepatt.group(1)
        return softdir

#Takes arguments from the command line
parser = argparse.ArgumentParser(description='Detection of CNVs')
#parser.add_argument('-d','--input_directory', default='', help='Directory containing all read files in fastq format. Either -d or -f (NOT BOTH) is required.')
parser.add_argument('-f','--input_file', default='', help='A single read file in fastq format. If paired end (-p flag) enter the paths of both files, separated by comma.') #Either -d or -f (NOT BOTH) is required.')
parser.add_argument('-r','--ref_genome', help='Indexed Reference genome fasta file (full file path).')
parser.add_argument('-u','--non_unique', help='Bed file containing the non-Unique regions of the mability track.')
parser.add_argument('-w','--window_size', default='500', help='Size of window to use. Default: ~500kb windows. If window file is not present, it will be created')
parser.add_argument('-p','--paired', help='Use this flag for paired-end data.  Omit for single-end data', action="store_true")
#parser.add_argument('-a','--aligned', help="Use this flag if the sequences are already aligned in bam format (file name is: <Sample Name>.bam)", action="store_true")
parser.add_argument('-b','--bwa_dir', default='', help='Directory in which the BWA aligner can be found')
parser.add_argument('-e','--bed_dir', default='', help='Directory in which bedtools can be found')
parser.add_argument('-s','--sam_dir', default='', help='Directory in which samtools can be found')
parser.add_argument('-c','--picard_dir', default='', help='Directory in which the picard tools can be found')
parser.add_argument('-t','--samstat_dir', default='', help='Directory in which the samstat program can be found')
parser.add_argument('-n','--num_threads', default = '1', help='Number of cpus available')
parser.add_argument('-m','--mem_limit', default = '10000000000', help='Amount of memory available')
parser.add_argument('-o','--bowtie', help='Use this flag to align with bowtie2 instead of the default BWA.', action="store_true")
parser.add_argument('-i','--bow_dir', default='', help='Directory in which the bowtie2 aligner can be found')
parser.add_argument('-sp','--species',default='',help='Species of data being analyzed')
parser.add_argument('-in','--input_type',default='',help='Type of input file being analyzed (\"fastq\" or \"bam\")')
parser.add_argument('-lO','--liftOver',default='',help='Use this flag to indicate that input files are BAMs aligned to mm10 that must be converted to mm9 using liftOver')
parser.add_argument('-en','--output',default='',help='Use this flag to indicate the final output location' )
#Checks that the minimum number of arguments
if len(sys.argv)<5:
        parser.print_help()
        sys.exit("Required arguments are missing!")

#Parse Arguments
args = parser.parse_args()
#Check if liftOver mm10-mm9 indicated
liftOver = args.liftOver
output =args.output
#theinputdir=args.input_directory
theinputfile=args.input_file
if theinputfile == '':
        sys.exit("You must supply input file(s) with -f flag (-f read1.fastq[,read2.fastq]).")
#if theinputdir == '' and theinputfile == '':
#       sys.exit("You must supply either a single fastq file (-f) or a directory containing fastq files (-d).")
#if theinputdir != '' and theinputfile != '':
#       sys.exit("You must supply either a single fastq file (-f) or a directory containing fastq files (-d), BUT NOT BOTH!")
species=args.species
if species == 'human' and species == 'mouse':
        sys.exit("You must indicate the species of the data, either \"human\" or \"mouse\".")
if species != 'human' and species != 'mouse':
        sys.exit("You must indicate the species of the data, either \"human\" or \"mouse\".")
input_type=args.input_type
if input_type == 'fastq' and input_type == 'bam':
        sys.exit("You must indicate one type of input data, either \"fastq\" or \"bam\".")
if input_type != 'fastq' and input_type != 'bam':
        sys.exit("You must indicate one type of input data, either \"fastq\" or \"bam\".")
if input_type == 'fastq':
        align = False
else:
        align = True
reffile=args.ref_genome
paired = args.paired
bowtie = args.bowtie
if os.path.isfile(workingpath + username+"_software.txt"):
        softfile = open(workingpath + username+"_software.txt",'r')
        softList = softfile.readlines()
        if bowtie == False:
                if len(softList) == 5:
                        bwadir = softList[0].strip()
                        beddir = softList[1].strip()
                        samdir = softList[2].strip()
                        statdir = softList[3].strip()
                        picdir = softList[4].strip()
                elif len(softList) == 6:
                        bwadir = softList[0].strip()
                        beddir = softList[1].strip()
                        samdir = softList[2].strip()
                        statdir = softList[3].strip()
                        picdir = softList[4].strip()
                        bowdir = softList[5].strip()
                else:
                        sys.exit("Software list has incorrect number of file paths.  Edit mysoftware.txt, or delete it from the working directory to allow the pipeline to automatically detect dependencies.")
        if bowtie == True:
                if len(softList) == 6:
                        bwadir = softList[0].strip()
                        beddir = softList[1].strip()
                        samdir = softList[2].strip()
                        statdir = softList[3].strip()
                        picdir = softList[4].strip()
                        bowdir = softList[5].strip()
                else:
                        sys.exit("Software list has incorrect number of file paths.  Edit mysoftware.txt, or delete it from the working directory to allow the pipeline to automatically detect dependencies.")

else:
        softfile = open(workingpath + username+"_software.txt",'w')
        bwadir = args.bwa_dir
        if bwadir=='':
                bwadir = find_software("bwa")
        beddir = args.bed_dir
        if beddir=='':
                beddir = find_software("bedtools")
        samdir = args.sam_dir
        if samdir=='':
                samdir = find_software("samtools")
        statdir = args.samstat_dir
        if statdir=='':
                statdir = find_software("samstat")
        picdir = args.picard_dir
        if picdir=='':
                picdir = find_software("MarkDuplicates.jar")
        if bowtie == False:
                softfile.write(bwadir +'\n' + beddir + '\n' + samdir + '\n' + statdir + '\n' + picdir)
                softfile.close()
        if bowtie == True:
                bowdir = args.bow_dir
                if bowdir=='':
                        bowdir = find_software("bowtie2")
                softfile.write(bwadir +'\n' + beddir + '\n' + samdir + '\n' + statdir + '\n' + picdir + '\n' + bowdir)
                softfile.close()
        sts = call("chmod 666 " + username+"_software.txt", shell=True)

numthreads = args.num_threads
#print numthreads
maxmem = args.mem_limit
nubed = args.non_unique
winsize = args.window_size
fullwinsize = int(winsize) * 1000
if species == 'human':
        winbed = "all.window.40mer.b37.map." + winsize + "kb.bed"
if species == 'mouse':
        winbed = "MM9.all.window.40mer.b37.map." + winsize + "kb.bed"
#Check if a window file of the appropriate size exists.  If so use it.  Otherwise create one.
if os.path.isfile(workingpath + winbed):
        print("Using window file: " + winbed)
else:
        if species == 'human':
                sts=call("python " + workingpath + "windowMaker.py -c " + workingpath + "human_g1k_v37.chrlist -n " + str(fullwinsize) + " -f " + workingpath + "b37_full.maskedMap.40mer.fa > " + workingpath + winbed,shell=True)
        if species == 'mouse':
                sts=call("python " + workingpath + "windowMaker.py -c " + workingpath + "mouse_mm9.chrlist -n " + str(fullwinsize) + " -f " + workingpath + "MM9.b37_full.maskedMap.40mer.fa > " + workingpath + winbed,shell=True)
        sts=call("chmod 666 " + workingpath + winbed, shell=True)
        print("Creating window file: " + winbed)
#align = args.aligned

#Generate a list file with full file path of all inputs.
#if theinputdir != '':
#       print "Input is a directory."
#       print theinputdir
#       sts=call("ls -d " + theinputdir + "/* | grep " + input_type + " > mylist.txt", shell=True)
#       inputfiles = open(workingpath + "mylist.txt",'r')
#       #Read in the fastq file list to a python list.
#       fastqList = inputfiles.readlines()
#       inputfiles.close()
if theinputfile != '':
        print("Input is a file.")
        if input_type == 'fastq':
                inputList=[]
                if paired == False:
                        inputList.append(theinputfile)
                if paired == True:
                        splitFile = theinputfile.split(',')
                        leftfile = splitFile[0]
                        rightfile = splitFile[1]
                        inputList.append(leftfile)
                        inputList.append(rightfile)
        if input_type == 'bam':
                inputList=[]
                inputList.append(theinputfile)
DEBUGFILE = open(workingpath + "debugfile.txt",'w')
#for testdebug in fastqList:
#       DEBUGFILE.write(testdebug + '\n')
#print winbed

#Get the file path
m=re.match('(.+/)(.*)',inputList[0])
path=m.group(1)
filename=m.group(2)
#print path

#Loop through the list of files, and perform the first step of a BWA alignment (bwa aln) on each.
if input_type == "fastq":
        for afastq in inputList:
                m=re.match('(.+/)(.*)',afastq)
#               m=re.match('(.*)(?=\.)',afastq)
                path=m.group(1)
                filename=m.group(2)
                breakName = filename.rstrip().split('.fastq')
                if len(breakName) ==2:
                        head=breakName[0]
                        tail='fastq' + breakName[1]
#               elif len(breakName) > 2:
#                       head=breakName[0]
#                       for anamepiece in breakName[1:-1]:
#                               head = head + "." + anamepiece
#                       tail = breakName[-1]
                else:
                        sys.exit("Fastq file names must contain one \".fastq\"")
                DEBUGFILE.write(path + '\n')
                DEBUGFILE.write(head + '\n')
                DEBUGFILE.write(tail + '\n')
#               if align ==False:
                if bowtie == False:
                        #print(bwadir + "/bwa aln -t " + numthreads + " " + reffile + " " + path + head + "." + tail + " -f " + path + head + ".sai")
                        sts=call(bwadir + "/bwa aln -t " + numthreads + " " + reffile + " " + path + head + "." + tail + " -f " + path + head + ".sai" , shell=True)
        DEBUGFILE.close()

if input_type == 'fastq':
        #If single end FASTQs, run bwa samse on each file.
        if paired ==False:
                bamSEFileList=[]
#               if align == False:
                print("Aligning single-end reads...")
                for afastq in inputList:
                        m=re.match('(.+/)(.*)',afastq)
                        path=m.group(1)
                        filename=m.group(2)
                        breakName = filename.rstrip().split('.fastq')
                        if len(breakName) == 2:
                                head=breakName[0]
                                tail='fastq' + breakName[1]
#               elif len(breakName) > 2:
#                       head=breakName[0]
#                       for anamepiece in breakName[1:-1]:
#                               head = head + "." + anamepiece
#                       tail = breakName[-1]
                        else:
                                sys.exit("Fastq file names must contain one \".fastq\"")
                        #print bwadir + "/bwa samse " + reffile + " " + path + head + ".sai " + path + head + "." + tail + " | " + samdir + " samtools view -Sb - > " + path + head + ".bam"
#                       if align == False:
                        if bowtie == False:
                                sts=call(bwadir + "/bwa samse " + reffile + " " + path + head + ".sai " + path + head + "." + tail + " | " + samdir + "/samtools view -Sb - > " + path + head + ".bam", shell=True)
                        if bowtie == True:
                                sts=call(bowdir + "/bowtie2 -x " + reffile + " -U " + path + head + tail + " -S "  + path + head + ".sam",shell=True)
                                sts=call(samdir + "/samtools view -Sb " + path + head + ".sam -o " + path + head + ".bam", shell=True)
                        thebamfile = path + head + ".bam"
                        bamSEFileList.append(thebamfile)
                theBamFileList=bamSEFileList
        #If paired-end reads, run bwa sampe.
        if paired ==True:
                bamPEFileList=[]
#               if align ==False:
                print("Aligning paired-end reads...")
                for j in range(0,len(inputList),2):
                        readpair1 = inputList[j].rstrip()
                        readpair2 = inputList[j+1].rstrip()
                        m1=re.match('(.+/)(.*)',readpair1)
                        m2=re.match('(.+/)(.*)',readpair2)
                        path1=m1.group(1)
                        filename1=m1.group(2)
                        path2=m2.group(1)
                        filename2=m2.group(2)
                        pieces1 =filename1.split(".fastq")
                        pieces2 =filename2.split(".fastq")
                        if len(pieces1) ==2:
                                head1=pieces1[0]
                                tail1="fastq" + pieces1[1]
#                       elif len(pieces1) > 2:
#                               head1=pieces1[0]
#                               for anamepiece in pieces1[1:-1]:
#                                       head1 = head1 + "." + anamepiece
#                               tail1 = pieces1[-1]
                        else:
                                sys.exit("Fastq file names must contain at least one \".fastq\"")
                        if len(pieces2) ==2:
                                head2=pieces2[0]
                                tail2="fastq" + pieces2[1]
#                       elif len(pieces2) > 2:
#                               head2=pieces2[0]
#                               for anamepiece in pieces2[1:-1]:
#                                       head2 = head2 + "." + anamepiece
#                               tail2 = pieces2[-1]
                        else:
                                sys.exit("Fastq file names must contain at least one \".fastq\"")
                        #print bwadir + "/bwa sampe " + reffile + " " + path1 + head1 + ".sai " + path2 + head2 + ".sai " + path1 + head1 + "." + tail1 + " " + path2 + head2 + "." + tail2 + " | " + samdir + " samtools view -Sb - > " + path1 + head1 + "_" + head2 + ".bam"
#                       if align == False:
                        if bowtie == False:
                                sts=call(bwadir + "/bwa sampe " + reffile + " " + path1 + head1 + ".sai " + path2 + head2 + ".sai " + path1 + head1 + "." + tail1 + " " + path2 + head2 + "." + tail2 + " | " + samdir + "/samtools view -Sb - > " + path1 + head1 + "_" + head2 + ".bam", shell=True)
                        if bowtie == True:
                                sts=call(bowdir + "/bowtie2 -x " + reffile + " -1 " + path1 + head1 + tail + " -2 " + path2 + head2 + tail + " -S "  + path1 + head1 + "_" + head2 + ".sam",shell=True)
                                sts=call(samdir + "/samtools view -Sb " + path1 + head1 + "_" + head2 + ".sam -o " + path1 + head1 + "_" + head2 + ".bam", shell=True)
                        thebamfile =  path1 + head1 + "_" + head2 + ".bam"
                        bamPEFileList.append(thebamfile)
                theBamFileList=bamPEFileList
if input_type == 'bam':
        theBamFileList=[]
        theBamFileList.append(theinputfile)
#print(theinputfile)
#print(input_type)
#print(theBamFileList)

for abam in theBamFileList:
        print(abam)
#       print(theinputfile)
        m=re.match('(.+/)(.*)',abam)
        path=m.group(1)
        filename=m.group(2)
        breakName = filename.rstrip().split('.bam')
        if len(breakName) ==2:
                bamhead=breakName[0]
        elif len(breakName) > 2:
                print("Something is off with this bam name. Too many \".bam\"s in filename? :")
                print(abam)
                bamhead=breakName[0]
                for anamepiece in breakName[1:-1]:
                        bamhead = bamhead + "." + anamepiece
        else:
                sys.exit("Bam file names must contain at least one \".\"")
        #Sort alignment files.
        print("Sorting bam file (samtools sort)...")
#       print samdir + "/samtools sort -m " + maxmem + " " + path + bamhead + ".bam " + path + bamhead + "_sorted"
        sts = call(samdir + "/samtools sort -m " + maxmem + " " + path + bamhead + ".bam " + path + bamhead + "_sorted", shell=True)
#       print "java -jar " + picdir + "/MarkDuplicates.jar INPUT=" + path + bamhead + "_sorted.bam OUTPUT=" + path + bamhead + "_sorted_rmdup.bam METRICS_FILE=" + path + bamhead + "_rmdup.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true"
        #Remove optical/pcr duplicates
        print("Removing duplicates with picard-tools...")
        sts = call("java -jar " + picdir + "/MarkDuplicates.jar INPUT=" + path + bamhead + "_sorted.bam OUTPUT=" + path + bamhead + "_sorted_rmdup.bam METRICS_FILE=" + path + bamhead + "_rmdup.metrics ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true", shell=True)
        print("Generating alignment statistics with flagstat")
        #Generate Alignment stats (flagstat and samstat)
        sts = call(samdir + "/samtools flagstat " + path + bamhead + "_sorted_rmdup.bam > " + path + bamhead + "_sorted_rmdup.flagstat", shell=True)
        sts = call(samdir + "/samtools flagstat " + path + bamhead + "_sorted.bam > " + path + bamhead + "_sorted.flagstat", shell=True)
#       sts = call(statdir + "/samstat " + path + bamhead + "_sorted_rmdup.bam",shell=True)
#       sts = call(statdir + "/samstat " + path + bamhead + "_sorted.bam",shell=True)
        #Parse all flagstats into one file
        #print "dupeStats_MOD.py -i " + path + bamhead + "_sorted.flagstat -d" +  path + bamhead + "_sorted_rmdup.flagstat >> " + path + "All_FlagStats.txt"
        sts = call("python " + workingpath + "dupeStats_MOD.py -i " + path + bamhead + "_sorted.flagstat -d" +  path + bamhead + "_sorted_rmdup.flagstat >> " + output +"/sample_stats/" + "All_FlagStats.txt", shell=True)
        #print samdir + "/samtools view -bu " + path + bamhead + "_sorted_rmdup.bam | " + beddir + "/bamToBed -i stdin | awk '{print $1\"\\t\"$2\"\\t\"$2+1}' | " + beddir + "/intersectBed -a stdin -b " + nubed + " -v | " + beddir + "/coverageBed -a stdin -b " + winbed + " -counts > " + path + bamhead + "_" + winsize + "kb.cov"
        #Index the sorted no-dupe bam file
        print("Indexing sorted, duplicate-free bam file...")
        sts = call(samdir + "/samtools index " + path + bamhead + "_sorted_rmdup.bam", shell=True)
        #Get read counts per window (coverageBed)
        #sts = call(samdir + "/samtools view -bu " + path + bamhead + "_sorted_rmdup.bam | " + beddir + "/bamToBed -i stdin | awk '{print $1\"\\t\"$2\"\\t\"$2+1}' | " + beddir + "/intersectBed -a stdin -b " + nubed + " -v | " + beddir + "/coverageBed -a stdin -b " + winbed + " -counts > " + path + bamhead + "_" + winsize + "kb.cov", shell=True)
        sts = call(samdir + "/samtools view -bu " + path + bamhead + "_sorted_rmdup.bam | " + beddir + "/bamToBed -i stdin | sed -e \"s/\\-1/0/g\" > " + path + bamhead + "_fullread.bed", shell=True)
        sts = call(samdir + "/samtools view -bu " + path + bamhead + "_sorted_rmdup.bam | " + beddir + "/bamToBed -i stdin | sed -e \"s/\\-1/0/g\" | awk '{print $1\"\\t\"$2\"\\t\"$2+1}' > " + path + bamhead + ".bed",shell=True)
        if liftOver == 'yes':
                #Adding lines to change mm10 to mm9 coordinates in BED file
                sts = call("mv " + path + bamhead + ".bed " + path + bamhead + "_mm10.bed", shell=True)
                sts = call("mv " + path + bamhead + "_fullread.bed " + path + bamhead + "_mm10_fullread.bed", shell=True)
                #Usage for liftOver: liftOver oldFile map.chain newFile unMapped
                sts = call("/nv/vol97/mcconnell_lab/cnvpipe/software/liftover/liftOver " + path + bamhead + "_mm10.bed /nv/vol97/mcconnell_lab/cnvpipe/software/liftover/mm10ToMm9.over.chain " + path + bamhead + ".bed " + path + bamhead + "_mm10to9unmapped.bed",shell=True)
                sts = call("/nv/vol97/mcconnell_lab/cnvpipe/software/liftover/liftOver " + path + bamhead + "_mm10_fullread.bed /nv/vol97/mcconnell_lab/cnvpipe/software/liftover/mm10ToMm9.over.chain " + path + bamhead + "_mm9_fullread.bed " + path + bamhead + "_mm10to9unmapped_fullread.bed",shell=True)
        sts=call(beddir + "/intersectBed -a "+ path + bamhead + ".bed"+ " -b " + nubed + " -v > " + path + bamhead + "_winIntersect.bed",shell=True)
        sts=call(beddir + "/coverageBed -a " + path + bamhead + "_winIntersect.bed" + " -b " + workingpath + winbed + " -counts > " + path + bamhead + "_" + winsize + "kb.cov", shell=True)
        #sts = call(beddir + "/bamToBed -i " + path + bamhead + "_sorted_rmdup.bam | awk '{print $1\"\\t\"$2\"\\t\"$2+1}' | sed -e \"s/\\-1//g\" | " + beddir + "/intersectBed -a stdin -b " + nubed + " -v | " + beddir + "/coverageBed -a stdin -b " + winbed + " -counts > " + path + bamhead + "_" + winsize + "kb.cov", shell=True)
        #New normalizing script.  Normalizes then segments with CBS.
        #print "Rscript CNV_countNorm.R " + path + " " + bamhead + "_" + winszie + "kb.cov " + bamhead + "_norm.num"
        if species == 'human':
                sts = call("cat "  + path + bamhead + "_" + winsize + "kb.cov | sed -e \"s/^chr//g\" | sed -e \"s/X/23/g\" | sed -e \"s/Y/24/g\" | sort -k 1 -n -k 2 -n > " + output +"/" +bamhead + "_" + winsize + "kb_SORT.cov", shell=True)
        if species == 'mouse':
                sts = call("cat "  + path + bamhead + "_" + winsize + "kb.cov | sed -e \"s/^chr//g\" | sed -e \"s/X/20/g\" | sed -e \"s/Y/21/g\" | sort -k 1 -n -k 2 -n > " + output+ "/"+ bamhead + "_" + winsize + "kb_SORT.cov", shell=True)
        sts = call("Rscript " + "/nv/vol97/mcconnell_lab/cnvpipe/TEST_TEST/" + "Normalize_read_counts_and_segment.R " + output + " " + path + " " + bamhead + "_" + winsize + "kb_SORT.cov" , shell=True)
        #Calculate confidence score.
        #sts = call("python " + workingpath + "residualsEdit.py -l " + output + "/"+ bamhead + "_" + winsize + "kb_SORT_norm.num.long -s 0.5 > " + path + "tempres.txt" , shell=True)
        #sts = call("cat " + path + "tempres.txt | grep " + str("\"#\"") + " >> " + output + "/sample_stats/" + "All_Conf_Scores.txt", shell=True)

#Change permissions for output to be accessible to all users
#sts = call("chmod 666 " + path + "*", shell=True)

#Move files


sts = call("mv " + path + bamhead + "_sorted_rmdup.bam " + output + "/bams", shell=True)
sts = call("mv " + path + bamhead + "_sorted_rmdup.bam.bai " + output + "/bams", shell=True)
if liftOver == "yes":
        sts = call("mv " + path + bamhead + "_mm10_fullread.bed " + output + "/beds", shell=True)
        sts = call("mv " + path + bamhead + "_mm9_fullread.bed " + output + "/beds", shell=True)
else:
        sts = call("mv " + path + bamhead + "_fullread.bed " + output + "/beds", shell=True)
sts = call("mv " + output + "/"+ bamhead + "_500kb_SORT.cov " + output + "/cov_files_badbinsIN", shell=True)
sts = call("mv " + output + "/" + bamhead + "_500kb_SORT_PLOT.png " + output + "/pngs_badbinsIN", shell=True)
sts = call("mv " + path + bamhead +"_sorted_rmdup.flagstat " + output + "/sample_stats", shell=True)
sts = call("mv " + path + bamhead + "_rmdup.metrics " + output + "/sample_stats", shell=True)
#sts = call("mv " + path + "All_Conf_Scores.txt " + path + "../sample_stats", shell=True)
#sts = call("mv " + path + "All_FlagStats.txt " + path + "../sample_stats", shell=True)
sts = call("mv " + output + "/" + bamhead +"_500kb_SORT_SoS.png " + output + "/sos_CN_mult_badbinsIN", shell=True)
sts = call("mv " + output + "/" + bamhead +"_500kb_SORT_norm.num.AllBinInfo.long " + output + "/full_longs_badbinsIN", shell=True)
#sts = call("mv " + path + "bin_CN_values.txt" + path + "../other_outputs", shell=True)
#sts = call("mv " + path + "bin_raw_counts.txt" + path + "../other_outputs", shell=True)
#sts = call("mv " + path + "BICs_and_other_info_badbinsIN.txt" + path + "../other_outputs", shell=True)
#sts = call("mv " + path + "SoS_CN_multiplier_info.txt" + path + "../sos_error_CN_multiplier", shell=True)


#Delete unneeded files as necessary
if input_type=="fastq":
        sts = call("rm " + path + bamhead + ".bam", shell=True)
        if paired == True:
                sts = call("rm " + path1 + head1 + ".sai", shell=True)
                sts = call("rm " + path2 + head2 + ".sai", shell=True)
        else:
                sts = call("rm " + path + head + ".sai", shell=True)
sts = call("rm " + path + bamhead + "_500kb.cov", shell=True)
sts = call("rm " + output +"/" + bamhead + "_500kb_SORT_norm.num", shell=True)
sts = call("rm " + output + "/" +bamhead + "_500kb_SORT_norm.num.long", shell=True)
sts = call("rm " + output + "/" +bamhead + "_500kb_SORT_norm.num.short", shell=True)
#Remove mm10 bed file in case of liftOver
if liftOver == 'yes':
        sts = call("rm " + path + bamhead + "_mm10.bed", shell=True)
        sts = call("rm " + path + bamhead + "_mm10to9unmapped.bed", shell=True)
        sts = call("rm " + path + bamhead + "_mm10to9unmapped_fullread.bed", shell=True)
#Remove main bed file (will contain mm9 coords. in case of liftOver)
sts = call("rm " + path + bamhead + ".bed", shell=True)
sts = call("rm " + path + bamhead + "_sorted.bam", shell=True)
#sts = call("rm " + path + bamhead + "_sorted.bam.samstat.html", shell=True)
sts = call("rm " + path + bamhead + "_sorted.flagstat", shell=True)
#sts = call("rm " + path + bamhead + "_sorted_rmdup.bam.samstat.html", shell=True)
sts = call("rm " + path + bamhead + "_winIntersect.bed", shell=True)
##sts = call("rm " + path + "tempres.txt", shell=True)
