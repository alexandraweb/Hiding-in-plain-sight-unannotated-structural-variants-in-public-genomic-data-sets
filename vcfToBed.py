import sys
import vcf
import time
import logging
import argparse
import os
import re
import subprocess


class MakeFastaFile:

    def __init__(self,inVCFFile,outFastaFile,chr,minInsertLength,logger):
        self.inVCFFile = inVCFFile
        self.outFile = outFastaFile
        self.chr = chr
        self.minInsertLength = minInsertLength
        self.logger = logger
        self.convertToFasta()

    def convertToFasta(self):
        start = time.time()
        if self.inVCFFile.endswith('.gz'):

            progressCounter = 0
            # vcf_reader = vcf.Reader(open(vcf_file, 'r'))
            vcf_reader = vcf.Reader(filename=self.inVCFFile)
            record = next(vcf_reader)
            all_coord = {}
            for record in vcf_reader.fetch(self.chr):

                # for record in vcf_reader:
                if (record.var_type != "snp") and (record.var_subtype == "ins"):
                    alt_seq = record.ALT
                    ref_seq = record.REF
                    chrom = record.CHROM
                    pos = record.POS
                    ref_seq = str(ref_seq).replace('[', '').replace(']', '')
                    alt_seq = str(alt_seq).replace('[', '').replace(']', '')
                    alt_seq_length = len(alt_seq) - len(ref_seq)
                    if alt_seq_length > self.minInsertLength:
                        coord = (">chr" + str(chrom) + ":" + str(pos))
                        if coord in all_coord:
                            all_coord[coord] += 1
                        else:
                            all_coord[coord] = 1
                        text = (coord + "_" + str(all_coord[coord]) + "\n" + alt_seq)
                        # print(text)
                        # add status print statements every 1e7 records
                        if progressCounter == 10000:
                            print ("10000 records reached")
                            progressCounter = 0
                        progressCounter += 1
                        with open(self.outFile, 'a') as fasta:
                            fasta.write(text + '\n')

            # print("done")
            end = time.time()
            program_time = (end - start)
            self.logger.info("The VCF file for chr" + self.chr + " has been converted to a fastA file")
            print (program_time)


class MakeBlastFile:

    def __init__(self,inFastFile,outBlastFile,blastDir,outMeiFile,chr,meiOnly,logger):
        self.inFastFile = inFastFile
        self.outBlastFile = outBlastFile
        self.outMeiFile = outMeiFile
        self.blastDir = blastDir
        self.chr = chr
        self.blastDb = self.blastDir + "/chr" + str(chr) + ".fa"
        self.blastMEI = self.blastDir + "/MEI_refs.fa"
        self.meiOnly = meiOnly
        self.logger = logger
        self.createBlastFile()

    def createBlastFile(self):
        # TODO: make blast dictionary configurable
        if self.meiOnly is False:    
            print("Using genome BLAST database: " + self.blastDb)

        print("Using MEI BLAST database: " + self.blastMEI)

        if self.meiOnly is False:
            blastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                            + ' -db ' + self.blastDb + ' -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -perc_identity 90 -out ' + self.outBlastFile \
                            + ' -outfmt "6 qseqid sseqid qcovs stitle qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand" '
            blast = os.system(blastCommand)

        meiBlastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                        + ' -db ' + self.blastMEI + ' -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -perc_identity 90 -out ' + self.outMeiFile \
                        + ' -outfmt "6 qseqid sseqid qcovs stitle qlen slen qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand" '

        # print(blastCommand)
        meiBlast = os.system(meiBlastCommand)
        self.logger.info("The fastA file for chr" + self.chr + " has been BLASTed")
        # print (program_time)
# 

class ApplyFilters:

    def __init__(self,inBlastFile,inMeiFile,outBEDFile,meiOnly,logger):

        self.inBlastFile = inBlastFile
        self.inMeiFile = inMeiFile
        self.outBEDFile = outBEDFile
        self.meiOnly = meiOnly
        self.logger = logger

        self.filter()
        # filter(self)
# 


    def filter(self):

        for i in range(2):
            if i==0:
                if self.meiOnly is True:
                    continue
                else:
                    blast_file = open(self.inBlastFile, 'r')
                    mei = False
            elif i==1:
                blast_file = open(self.inMeiFile, 'r')
                mei = True 
            else:
                break

            # print("here")
            blast_file_review = {}
            count = 0

            # print ("qseqid", "\t", "bitscore", "\t", "(length/qlen)", "\t", "pident", "\t", "gaps", "\t", "sstrand", "\t", "send", "\t", "sstart\n")
            # for i in self.inBlastFile:
            for i in blast_file:
                if not i.startswith("#"):
                    count +=1
                    items = i.split("\t")
                    qseqid = items[0]
                    sseqid = items[1]
                    qcovs = items[2]
                    stitle = items[3]
                    qlen = float(items[4])
                    slen = items[5]
                    qstart = items[6]
                    qend = items[7]
                    sstart = items[8]
                    send = items[9]
                    qseq = items[10]
                    sseq = items[11]
                    evalue = items[12]
                    bitscore = float(items[13])
                    score = items[14]
                    length = float(items[15])
                    pident = float(items[16])
                    nident = items[17]
                    mismatch = items[18]
                    positive = items[19]
                    gapopen = items[20]
                    gaps = float(items[21])
                    ppos = items[22]
                    sstrand = items[23].replace("\n", "")
                    chr = re.findall('^chr[0-9]+',qseqid)
                    chr = chr[0]
                    # if (sstart > send) and (length):
                    #	blast_file_review[qseqid] = []

                    # chr21:47836788_
                    processed_qseqid = int(qseqid.replace(chr + ":", "").split("_")[0])

                    bp_length = int(sstart) - int(processed_qseqid)
                    blast_length = (int(length) / int(qlen))
                    mei_start = int(processed_qseqid) - 1
                    # print(qseqid, processed_qseqid, sstart, bp_length)

                    Deletionfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
                    Inversionfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
                    TandomDupfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)

                    if (int(bitscore) > 50) and (blast_length > 0.9) and (float(pident) > 90) and (int(gaps) < 3):
                        if mei is False:
                            self.deletions(chr,Deletionfamily_strand_uniqueNumber,sstrand,send,sstart,
                                                   processed_qseqid,bitscore,bp_length)

                            self.inversions(chr,bitscore,send,sstart,Inversionfamily_strand_uniqueNumber,sstrand,bp_length,processed_qseqid)

                            self.tandumDup(chr,bitscore,send,sstart,sstrand,qseqid,TandomDupfamily_strand_uniqueNumber,bp_length,processed_qseqid)

                        else:
                            self.meis(chr,mei_start,processed_qseqid,sstart,bitscore,sseqid,sstrand,count,send)
                        #print qseqid, "\t", bitscore, "\t", blast_length, "\t", pident, "\t", gaps, "\t", sstrand, "\t", send, "\t", sstart
                if count % 1000 == 0:
                    print('Done with BLAST result ' + str(count) )
        self.logger.info("The BLAST results have been filtered out output as BED")
            # print (program_time)


    def deletions(self,chr,Deletionfamily_strand_uniqueNumber,sstrand,send,sstart,processed_qseqid,bitscore,bp_length):
        #require that alignment is on plus strand nad downstream from org site
        deletionsFile = open(str(chr) + '_deletions.bed','a+')
        # print(sstrand, bp_length)
        if (sstrand.strip().rstrip() == 'plus') and (int(bp_length) > 50) and (int(bp_length) < 1e5):
            deletionsFile.write(chr + "\t" +  str(sstart) +  "\t" +  str(send)  +  "\t" + Deletionfamily_strand_uniqueNumber + "\t" + str(bitscore) + "\t" + "deletion" + "\n")
        deletionsFile.close()

    def meis(self,chr,mei_start,processed_qseqid,sstart,bitscore,sseqid,sstrand,count,send):
        #stop is org position of inserted sequence that is in query id
        #start is stop -1

        meisFile = open(str(chr) + '_meis.bed', 'a+')
        MEIfamily_strand_uniqueNumber = sseqid + "_" + sstrand + "_" + str(count)
        meisFile.write(chr + "\t" + str(mei_start) + "\t" + str(processed_qseqid) + "\t" + MEIfamily_strand_uniqueNumber + "\t" + str(bitscore) + "\t" + "MEI" + "\n")
        meisFile.close()

    def inversions(self,chr,bitscore,send,sstart,sstrand,Inversionfamily_strand_uniqueNumber,bp_length,processed_qseqid):
        inversionsFile = open(str(chr) +  '_inversions.bed','a+')
        #strand is minus
        if ( (sstrand.strip().rstrip() == 'minus') and (abs(bp_length) > 50 ) and (abs(bp_length) < 1e5) ):
            start = min([sstart, processed_qseqid])
            end = max([sstart, processed_qseqid])
            inversionsFile.write(chr +  "\t" +  str(start) +  "\t" +  str(end) +  "\t" +  Inversionfamily_strand_uniqueNumber +  "\t" +  str(bitscore) +  "\t" +  "inversion" + "\n")
        inversionsFile.close()
    def tandumDup(self,chr,bitscore,send,sstart,sstrand,qseqid,TandomDupfamily_strand_uniqueNumber,bp_length,processed_qseqid):
        # sstart < qseqid
        # plus strand
        tandumDupFile = open(str(chr) + '_tandumDup.bed', 'a+')
        if (sstrand.strip().rstrip() == 'plus') and (int(sstart) < processed_qseqid) and (abs(bp_length) > 50) and (abs(bp_length) < 1e5):
            tandumDupFile.write(chr +  "\t" +  str(sstart) +  "\t" +  str(processed_qseqid) +  "\t" +  TandomDupfamily_strand_uniqueNumber +  "\t" +  str(bitscore) +  "\t" +  "tandomDup" + "\n")
        tandumDupFile.close()



def main():

    #inputfile = ''
    #outputfile = ''
    #minInsertLength = 25
    # chromosome='21'


    parser = argparse.ArgumentParser(prog='vcfToFasta.py',
                                     description='''Convert a vcf file to a fasta file''')

    parser.add_argument('-inFile', '--inputFile', type=str, required=False,
                        help='''VCF file to be parsed''')

    parser.add_argument('-outFile', '--outputFile', type=str, required=False,
                        help='''Name of fasta file that gets outputted''')

    parser.add_argument('-chr', '--chromosome', type=str, required=False,
                        help='''Chromosome''')

    parser.add_argument('-dir', '--blastDir', type=str, required=False,
                        help='''BLAST database directory''')

    parser.add_argument('-meiOnly', '--meiOnly', type=bool, required=False,
                        help='''Analyze MEIs only''')

    parser.add_argument('-v', '--version', action='version', version='0.0.1 ')

    args = parser.parse_args()

    meiOnly = False
    inVCFFile = args.inputFile
    outFastaFile = args.outputFile
    chr = args.chromosome
    minInsertLength = 25
    blastDir = args.blastDir
    outBlastFile = chr + 'Blast.txt'
    outMeiFile = chr + 'BlastMEI.txt'
    outBEDFile = chr + 'BedFile.bed'
    meiOnly = args.meiOnly


    # create logger
    logger = logging.getLogger('vcf_to_fasta')
    logger.setLevel(logging.DEBUG)

    # create console handler
    vcf2Fasta = logging.StreamHandler()
    vcf2Fasta.setLevel(logging.DEBUG)

    # create formatter and add it to the handler
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    vcf2Fasta.setFormatter(formatter)

    # add the handler to the logger
    logger.addHandler(vcf2Fasta)
    # print ("here")
    try:
        outputfile = chr + "FastA_AS.fa"
        # simple test to try and open file
        f = open(inVCFFile, 'r')
        print ('Input file is ', inVCFFile)
        print ('Output file is ', outFastaFile)
        print ('Minimum insert length is ', minInsertLength)
        print ('Chromosome is ', str(chr))

    except IOError as e:
        logger.critical('Cannot open file.\n%s' % e)
        exit(1)

    #convert VCF to Fasta file
    MakeFastaFile(inVCFFile,outFastaFile,chr,minInsertLength,logger)

    # take newly created fasta file and create a blast file and a mei blast file
    inFastFile = outFastaFile
    MakeBlastFile(inFastFile,outBlastFile,blastDir,outMeiFile,chr,meiOnly,logger)

    #take newly created blast file and mei blast file and apply filters to create bed file
    inBlastFile = outBlastFile
    inMeiFile = outMeiFile

    # blast_file = open('/Volumes/bioinfo/users/Alexandrea/hack-a-thon/blastchr21_output6.txt', 'r')
    ApplyFilters(inBlastFile,inMeiFile,outBEDFile,meiOnly,logger)
    # ApplyFilters(blast_file) #inBlastFile, inMeiFile, outBEDFile)

if __name__ == "__main__":
    main()