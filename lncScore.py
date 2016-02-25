#!/usr/bin/env python

'''determine whether the version of user's python comply with the requirements of  this procedure'''
import sys
if sys.version_info[0] != 2 or sys.version_info[1] != 7:
	print >>sys.stderr, "\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " CPAT needs python2.7!\n"
	sys.exit()

import os
import re
import math
import string
import optparse
import time
#import pysam
from multiprocessing import Process
import shutil

from cpmodule import fickett
from cpmodule import orf
from cpmodule import ireader
import numpy as np
from sklearn.linear_model import LogisticRegressionCV
 
#==============================================================================
# functions definition
#==============================================================================
def bed_or_fasta(infile):
	'''determine if the input file is bed or fasta format'''
	format = "UNKNOWN"
	for line in ireader.reader(infile):
		#line = line.strip()
		if line.startswith('#'):
			continue
		if line.startswith('>'):
			format="FASTA"
			return format
		elif len(line.split())==12:
			format='BED'
			return format
	return format
#==============================================================================
# checkout the availability of the file
#==============================================================================
def index_fasta(infile):
	if os.path.isfile(infile):
		pass
	else:
		print >>sys.stderr, "Indexing " + infile + ' ...',
		pysam.faidx(infile)
		print >>sys.stderr, "Done!"
#==============================================================================
# transfer the .bed file into a .fasta file
#==============================================================================
def bed_to_fasta(inbed,refgenome):
	'''extract features of sequence from bed line'''
		
	transtab = string.maketrans("ACGTNX","TGCANX")
	mRNA_seq = ''
	if inbed.strip():
		try:
			fields = inbed.split()
			chrom = fields[0]
			tx_start = int( fields[1] )
			geneName = fields[3]
			strand = fields[5].replace(" ","_")			
			exon_starts = map(int, fields[11].rstrip( ',\n' ).split( ',' ) )
			exon_starts = map((lambda x: x + tx_start ), exon_starts)
			exon_ends = map( int, fields[10].rstrip( ',\n' ).split( ',' ) )
			exon_ends = map((lambda x, y: x + y ), exon_starts, exon_ends);   
		except:
			print >>sys.stderr,"Wrong format!" + inbed 
			return None
		for st,end in zip(exon_starts, exon_ends):
			exon_coord = chrom + ':' + str(st +1) + '-' + str(end)
			tmp = pysam.faidx(refgenome,exon_coord)
			mRNA_seq += ''.join([i.rstrip('\n\r') for i in tmp[1:]])
		if strand =='-':
			mRNA_seq = mRNA_seq.upper().translate(transtab)[::-1]				
		return (geneName, mRNA_seq)
#==============================================================================
# transfer multiple-line fasta format into two-line fasta format
#==============================================================================
def TwoLineFasta (Seq_Array):
    Tmp_sequence_Arr = []
    Tmp_trans_str = ''
    for i in range(len(Seq_Array)):
        if '>' in Seq_Array[i]:
            if i == 0:
                Tmp_sequence_Arr.append(Seq_Array[i])
            else:
                Tmp_sequence_Arr.append(Tmp_trans_str)
                Tmp_sequence_Arr.append(Seq_Array[i])
                Tmp_trans_str = ''
        else:
            if i == len(Seq_Array) - 1:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
                Tmp_sequence_Arr.append(Tmp_trans_str)
            else:
                Tmp_trans_str = Tmp_trans_str + str(Seq_Array[i])
    return Tmp_sequence_Arr
#==============================================================================
# Remove transcripts shorter than 200nt
#==============================================================================
   
def Check_length (input_arr,temp_log):
    LogResult = temp_log + '.log'
    LogFile = open(LogResult,'w')
    LogFile.write('Transcripts shorter than 200nt:\n')
    Tmp_Arr = []
    for n in range(len(input_arr)):
        if n == 0 or n % 2 == 0:
            label = input_arr[n]
        else :
            seq = input_arr[n]
            if len(seq.lower().replace('n', ''))> 200:
                Tmp_Arr.append(label)
                Tmp_Arr.append(seq)
            else:
                TempLabel = label.split(' ')
                LogString = TempLabel[0] +'\n'
                LogFile.write(LogString)
                LogString = seq +'\n'
                LogFile.write(LogString)
    LogFile.close()
    return Tmp_Arr
#==============================================================================
# construct id array and sequence array
#==============================================================================
   
def Tran_Seq (input_arr):
    label_Arr = []
    FastA_seq_Arr = []
    for n in range(len(input_arr)):
        if n == 0 or n % 2 == 0:
            label = input_arr[n]
            label_Arr.append(label)
        else :
            seq = input_arr[n]
            FastA_seq_Arr.append(seq)
    return (label_Arr,FastA_seq_Arr)
#==============================================================================
# To split a file into multiple files, used for parallel computing 
#==============================================================================
def split(files,totallen,number,out,filename):
    seq_num = len(files)/2
    split_step = int(int(totallen) / int(number))
    title = ''+out+ '/'+ filename
    [labels,Sequences] = Tran_Seq(files)    
    start = 0
    end = seq_num
    length = 0
    for i in range(1,int(number)):
        temp_title = title + str(i)
        TEMP_FILE = open(temp_title,'w')
        for j in range(start,end):
            Tmp = labels[j]
            Tmp = str(Tmp) + '\n'
            TEMP_FILE.write(Tmp)
            Tmp = Sequences[j]
            length = length + len(Tmp)
            Tmp = str(Tmp) + '\n'
            TEMP_FILE.write(Tmp)
            if (length > split_step) :
                break
        TEMP_FILE.close()
        start = j + 1
        length = 0
    temp_title = title + str(number)
    TEMP_FILE = open(temp_title,'w')
    for j in range(start,end):
        Tmp = labels[j]
        Tmp = str(Tmp) + '\n'
        TEMP_FILE.write(Tmp)
        Tmp = Sequences[j]
        Tmp = str(Tmp) + '\n'
        TEMP_FILE.write(Tmp)
    TEMP_FILE.close()
#==============================================================================
# transfer the transcript sequence into codons 
#==============================================================================
def InitCodonSeq(num,length,step,Arr):
    TempStrPar = ''
    for w in range(num,length,step):
        index = w
        code1 = Arr[index]
        index += 1
        code2 = Arr[index]
        index += 1
        code3 = Arr[index]
        Temp = code1+code2+code3
        TempStrPar = TempStrPar+Temp+' '
    return TempStrPar
#==============================================================================
# Searching for the maximum coding subsequence (MCSS)    
#==============================================================================
  
def MaxInterval(TempArray,seqLength,coding,noncoding):
    maxStart = 0
    maxEnd = 0
    Max_cur = 0
    Max_so = 0
    MaxStemp = maxStart
    for n in range(0,seqLength-1):
        k = TempArray[n] + TempArray[n+1]
        if (not coding.has_key(k)) or (not noncoding.has_key(k)):
            continue
        if coding[k] >0 and noncoding[k] > 0:
            codingscore = math.log(coding[k]/noncoding[k])
            Max_cur = Max_cur + codingscore
            if codingscore >= Max_cur:
                Max_cur = codingscore
                MaxStemp = n
        elif coding[k] > 0 and noncoding[k] == 0:
            Max_cur = Max_cur + 10
        elif coding[k] == 0 :
            Max_cur = 0
            MaxStemp = n
        else:
            continue                
        
        if Max_cur > Max_so:
            Max_so = Max_cur
            maxEnd = n
            maxStart = MaxStemp

    return(Max_so,maxStart,maxEnd)
#==============================================================================
# The process for the calculation of MCSS features
#==============================================================================
def mcssProcess(tran_sec_seq,coding,noncoding):
    sequence_process_Arr = list(tran_sec_seq)
    Seq_len = len(sequence_process_Arr) - 1
    max_coding_Value = []

    max_coding_String = []
    coding_length_store_array = []

    for o in range(0,3): #three kinds of open reading frame of each sequence
        TempStr = ''
        TempStr = InitCodonSeq(o,Seq_len-1,3,sequence_process_Arr)
        TempArray = TempStr.split(' ') # construct codon array
        TempArray.pop()
        seqLength = len(TempArray)
        (Max_so,maxStart,maxEnd) = MaxInterval(TempArray,seqLength,coding,noncoding)
        max_coding_Value.append(Max_so)
        OutStr_coding = ''
        for out in range(maxStart,maxEnd+1):
            OutStr_coding = OutStr_coding+TempArray[out]+' '
        
        max_coding_String.append(OutStr_coding)
        coding_length_store_array.append(maxEnd+1-maxStart)
  
    TotalCodingScore = sum(max_coding_Value)
    MaxCodingScore = max(max_coding_Value)
    indexCoding = max_coding_Value.index(MaxCodingScore)
    CodingSequenceLen = coding_length_store_array[indexCoding]
    CodingPercent = MaxCodingScore/float(TotalCodingScore)

    return(MaxCodingScore,CodingSequenceLen,CodingPercent)    
  
#==============================================================================
# The process for the calculation of hexamer score and distance
#==============================================================================
  
def HexamerFeatures(seq,hash_matrix):
    if len(seq) < 6:
		return(0,0)
    frame_sequence = list(seq)
    frame_seq_length = len(frame_sequence)
    CDS_array = []
    for o in range(0,3):
        frame_TempStr = ''
        frame_TempStr = InitCodonSeq(o,frame_seq_length-2,3,frame_sequence)
        frame_array = frame_TempStr.split(' ') ## codon array
        frame_array.pop()
        other_num = 0
        frame_array_Len = len(frame_array) - 1
        for j in range(frame_array_Len):
            temp2 = frame_array[j]+frame_array[j+1]
            temple4 = re.compile('[atcg]{6}')
            if temple4.match(temp2):
                other_num = string.atof(other_num) + string.atof(hash_matrix[temp2])
        frame_array_Len = frame_array_Len + 2
        other_num = other_num / frame_array_Len
        CDS_array.append(other_num)      
    Mscore = max(CDS_array)
    score_distance = 0
    for m in range(0,3): #problem location
        score_distance += Mscore - CDS_array[m]
    score_distance = score_distance/float(2)   

    return(Mscore,score_distance)

#==============================================================================
#   Print the final result 
#==============================================================================
def PrintResult(ids,labels,probability,outputfile):

        
    Tabel = 'Transcript_id' + '\t' + 'Index' + '\t' + 'Coding_score' + '\n'
    outputfile.write(Tabel)
    for i in range(len(ids)):
        transcriptid = ids[i]
        Arr_label = transcriptid.split('>')
        line = Arr_label[1]
        if labels[i] == 1:
            line = line + '\t' + 'coding'
        else:
            line = line + '\t' + 'noncoding'
        line = line + '\t' + str(probability[i]) +'\n'
        outputfile.write(line)
        
   
#==============================================================================
# The main process for the feature calculation
#==============================================================================
def mainProcess(input,output,number,c_tab,g_tab,codonArr,hash_matrix,classifier):

    if number > 1:
        Temp_Dir = output + '_Tmp_Dir'
        temp_score = ''+Temp_Dir+'/'+ output + str(number)
#        temp_feature = ''+Temp_Dir+'/temp_feature' + str(number)
        SCORE = open(temp_score,'w')
#        DATA = open(temp_feature,'w')
        sequence_Arr = input.split('\n')
        sLen = len(sequence_Arr) - 1
        del sequence_Arr[sLen]
    if number == 1:
        SCORE = open(output,'w')
        sequence_Arr = input        
    
    label_Arr_tmp = []
    FastA_seq_Arr_tmp = []
    for n in range(len(sequence_Arr)):
        if n == 0 or n % 2 == 0:
            label = sequence_Arr[n]
            label_Arr_tmp.append(label)
        else :
            seq = sequence_Arr[n]
            FastA_seq_Arr_tmp.append(seq)
    data = []
    ids = []
    for i in range(len(label_Arr_tmp)):
        Seq = FastA_seq_Arr_tmp[i]
        tran_fir_seq = Seq.lower()
        tran_sec_seq_one = tran_fir_seq.replace('u','t')
        strinfo = re.compile('[^agctn]')                   
        tran_sec_seq = strinfo.sub('n',tran_sec_seq_one)                
        tran_sec_seq2 = tran_sec_seq.upper()
        tmp = orf.ORFFinder(tran_sec_seq2)
        (CDS_start, CDS_stop, CDS_size, CDS_frame, CDS_seq) = tmp.longest_orf(direction="+")
        (MCS,CSL,CP) = mcssProcess(tran_sec_seq2,c_tab,g_tab)
        fickett_score = fickett.fickett_value(CDS_seq)
        (orfscore,orfdistance) = HexamerFeatures(CDS_seq.lower(),hash_matrix)
        labels_Arr = label_Arr_tmp[i].split()
        ids.append(labels_Arr[0])
        Exons_mscore = []
        Exons_distance =[]
        Exons_GC = []

        Site_start = 0
        for j in range(1,len(labels_Arr)):
            
            seq = tran_sec_seq[Site_start:Site_start+int(labels_Arr[j])]
            if (len(seq) > 0):
                GCnum = seq.count('c') + seq.count('g')
                GCratio = GCnum/float(len(seq))          
                Exons_GC.append(GCratio)
                (mscore,distance) = HexamerFeatures(seq,hash_matrix)
                Exons_mscore.append(mscore)
                Exons_distance.append(distance)
                Site_start = Site_start + int(labels_Arr[j])
            else:
                continue
        Max_Mscore_exon = max(Exons_mscore)
        Max_distance = max(Exons_distance)
        Max_GCcontent = max(Exons_GC)
       
        full_len = len(tran_sec_seq)
        orf_ratio = CDS_size/float(full_len)
        
        transcript_features = [CDS_size,orf_ratio,fickett_score,orfscore,orfdistance,Max_Mscore_exon,Max_distance,Max_GCcontent,MCS,CSL,CP]        
        data.append(transcript_features)
#        PROPERTY_STR = labels_Arr[0]  + ' ' + str(CDS_size) + ' '+ str(orf_ratio) + ' ' + str(fickett_score) + ' '+ str(orfscore) + ' '+ str(orfdistance)+' '+ str(Max_Mscore_exon)+ ' ' + str(Max_distance)+ ' ' + str(Max_GCcontent)+ ' ' +str(MCS) +' '+str(CSL)+' '+str(CP)+'\n'
#        DATA.write(PROPERTY_STR)
    testing_data = np.array(data)
    del data
    testing_data = testing_data.reshape(len(label_Arr_tmp),11)
    prob = classifier.predict_proba(testing_data)
    labels = classifier.predict(testing_data)
    PrintResult(ids,labels,prob[:,1],SCORE) 
    SCORE.close()
#    return(PROPERTY_ARR)
#==============================================================================
# Main program
#==============================================================================

parse=optparse.OptionParser()
parse.add_option('-f','--file',dest='file',action='store',metavar='input files',help="enter transcripts in .bed or .fasta format: if this is a .bed format file, '-r' must be specified; if this is a .fasta format file, ignore the '-r'.")
parse.add_option('-g','--gtf',dest='gtf',action='store',metavar='gtf file name',help='please enter your gtf files')
parse.add_option('-o','--out',dest='outfile',action='store',metavar='output files',help='assign your output file')
parse.add_option('-p','--parallel',dest='parallel',action='store',metavar='prallel numbers',default=1,help='please enter your specified speed ratio')
parse.add_option('-x','--hex',dest='hexamer_dat',action="store",metavar='hexamer matrix',help="Prebuilt hexamer frequency table (Human, Mouse, Fly, Zebrafish, C. elegans, Sheep and Rat). Run 'make_hexamer_tab.py' to make this table out of your own training dataset.")
parse.add_option('-t','--train',dest='training_dat',action="store",metavar='training dataset',default="dat/Human_training",help="Please enter your specified training dataset")
parse.add_option("-r","--ref",dest="ref_genome",action="store",metavar='reference genome files',help="Reference genome sequences in FASTA format. Ignore this option if sequences file was provided to '-f'. Reference genome file will be indexed automatically (produce *.fa file along with the original *.bed file within the same directory) if hasn't been done.")

(options,args) = parse.parse_args()

#check input and output files
for file in ([options.file,options.gtf,options.outfile,options.hexamer_dat,options.training_dat]):
	if not (file):
		parse.print_help()
		sys.exit(0)

file_format = bed_or_fasta(options.file)
if file_format == 'UNKNOWN':
	print >>sys.stderr, "\nError: unknown file format of '-g'\n"
	parse.print_help()
	sys.exit(0)		
elif file_format == 'BED':
    import pysam
    print >>sys.stderr, "Input gene file is in BED format"
    if not options.ref_genome:
		print >>sys.stderr, "\nError: Reference genome file must be provided\n"
		parse.print_help()
		sys.exit(0)
    index_fasta(options.ref_genome)
    filearray = options.file.split('.')
    inPutFileName = filearray[0] + '.fasta'
    TMP = open(inPutFileName,'w')
    for line in ireader.reader(options.file):
		if line.startswith('track'):continue
		if line.startswith('#'):continue
		if line.startswith('browser'):continue
		#if not line.strip(): continue
		(gene_id, sequence)=bed_to_fasta(line, options.ref_genome)   
		print >>TMP, '\n'.join([str(i) for i in [gene_id, sequence]])

		
elif file_format == 'FASTA':
    inPutFileName = options.file

LNCSCOREPATH = os.path.split(os.path.realpath(__file__))[0]

temp_inPutFileName = 'temp_inputfile.fasta'
os.system('perl '+ LNCSCOREPATH + '/cpmodule/shortID.pl '+ inPutFileName + ' '+ temp_inPutFileName)
exon_inPutFileName = 'inputfile.fasta'
os.system('perl '+ LNCSCOREPATH + '/cpmodule/exon_extraction.pl '+ temp_inPutFileName +' ' + options.gtf + ' '+ exon_inPutFileName)
os.remove(temp_inPutFileName)
inPutFileName = exon_inPutFileName

outPutFileName = options.outfile
Parallel = options.parallel


#==============================================================================

MatrixPath = LNCSCOREPATH + "/dat/Matrix"
inMatrix = open(MatrixPath)
Matrix = inMatrix.read()
inMatrix.close()

Alphabet = ['ttt','ttc','tta','ttg','tct','tcc','tca','tcg','tat','tac','tgt','tgc','tgg','ctt','ctc','cta','ctg','cct','ccc','cca','ccg','cat','cac','caa','cag','cgt','cgc','cga','cgg','att','atc','ata','atg','act','acc','aca','acg','aat','aac','aaa','aag','agt','agc','aga','agg','gtt','gtc','gta','gtg','gct','gcc','gca','gcg','gat','gac','gaa','gag','ggt','ggc','gga','ggg']
Matrix_hash = {}
Matrix_Arr=Matrix.split('\n')
length = len(Matrix_Arr) - 1
del Matrix_Arr[length]
for line in Matrix_Arr :
    each = line.split('\t')
    key = each[0]
    value = each[1]
    Matrix_hash[key] = value
#==============================================================================
coding={}
noncoding={}	
for line in open(options.hexamer_dat):
	line = line.strip()
	fields = line.split()
	if fields[0] == 'hexamer':continue
	coding[fields[0]] = float(fields[1])
	noncoding[fields[0]] =  float(fields[2])
#==============================================================================
#  To built a logit model with the training dataset
#==============================================================================
f = open(options.training_dat)
data = np.loadtxt(f)
training_data = data[:,1:]
labels = data[:,0]
classifier = LogisticRegressionCV().fit(training_data,labels)   
#################################### 64 alphabet and hash dictionary ############################################
inFiles = open(inPutFileName)
inFilesArr = inFiles.read()
inFiles.close()

Compute_time = time.time()
#==============================================================================


if int(Parallel) == 1:
    sequence_Arr = inFilesArr.split('\n')
    del inFilesArr
    sLen = len(sequence_Arr) - 1
    del sequence_Arr[sLen]
    ARRAY_temp =  TwoLineFasta(sequence_Arr)
    ARRAY = Check_length(ARRAY_temp,outPutFileName)
    del ARRAY_temp
    inFileLength = len(ARRAY)/2
    del sequence_Arr
    mainProcess(ARRAY,outPutFileName,1,coding,noncoding,Alphabet,Matrix_hash,classifier)
    print('lncScore: The calculation of cpat features wws completely done!')
    print("%f second for" % (time.time() - Compute_time) + ' ' + str(inFileLength) + ' ' + "transcript's computation.")
#    shutil.rmtree(Temp_Dir,True)   


if int(Parallel) > 1:
    sequence_Arr = inFilesArr.split('\n')
    del inFilesArr
    sLen = len(sequence_Arr) - 1
    del sequence_Arr[sLen]
    ARRAY_temp =  TwoLineFasta(sequence_Arr)
    del sequence_Arr
    ARRAY = Check_length(ARRAY_temp,outPutFileName)
    del ARRAY_temp
    Label_Array,FastA_Seq_Array = Tran_Seq(ARRAY)
    inFileLength = len(Label_Array)
    TOT_STRING = []
    totallen = 0
    for i in range(len(Label_Array)):
        tmp_label_one = Label_Array[i]
        tmp_label = tmp_label_one.replace('\r','')
        tmp_seq = FastA_Seq_Array[i]
        Temp_Seq = tmp_seq.replace('\r','')
        TOT_STRING.append(tmp_label)
        TOT_STRING.append(Temp_Seq)
        totallen = totallen + len(Temp_Seq)
    del Label_Array
    del FastA_Seq_Array  
    Proc_Thread = []
    Temp_Dir = outPutFileName + '_Tmp_Dir'
    os.mkdir(Temp_Dir)
    split(TOT_STRING,totallen,Parallel,Temp_Dir,'sequence_file')
    del TOT_STRING
    for i in range(1,int(Parallel)+1):
        temp_inPutFileName = ''+Temp_Dir+'/sequence_file' + str(i)
        temp_inFiles = open(temp_inPutFileName)
        temp_inFilesArr = temp_inFiles.read()
        Proc_Thread.append(Process(target=mainProcess, args=(temp_inFilesArr,outPutFileName,str(i),coding,noncoding,Alphabet,Matrix_hash,classifier)))
    for p in Proc_Thread:
        p.start()
    for i in Proc_Thread:
        p.join()
    feature_string = ''
    data_string = ''

    Files = open(outPutFileName,'w')
    features = ''
    i = 1
    while i < int(Parallel)+1:
        feature_string = ''+Temp_Dir+'/'+outPutFileName + str(i)
        tempfile = open(feature_string)
        tempfeatures = tempfile.read()
        if len(tempfeatures) > 0:
            features = features+tempfeatures
            i = i+1
        tempfile.close()
    Files.write(features)
    Files.close()
   
    print('lncScore: The classification of candidate transcripts was completely done!')
    print("%f second for" % (time.time() - Compute_time) + ' ' + str(inFileLength) + ' ' + "transcript's computation.")
    shutil.rmtree(Temp_Dir,True)   

os.remove(inPutFileName)


