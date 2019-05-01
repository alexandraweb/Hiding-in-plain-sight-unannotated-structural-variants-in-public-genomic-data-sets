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
                            print "10000 records reached"
                            progressCounter = 0
                        progressCounter += 1
                        with open(self.outFile, 'a') as fasta:
                            fasta.write(text + '\n')

            # print("done")
            end = time.time()
            program_time = (end - start)
            self.logger.info("The VCF file for chr" + self.chr + " has been converted to a fastA file")
            print program_time


class MakeBlastFile:

    def __init__(self,inFastFile,outBlastFile,outMeiFile):
        self.inFastFile = inFastFile
        self.outBlastFile = outBlastFile
        self.outMeiFile = outMeiFile
        self.createBlastFile()

    def createBlastFile(self):

        blastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                        + ' -db chr21db  -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -out ' + self.outBlastFile \
                        + ' -outfmt "6 qseqid sseqid qcovs stitle qlen qstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand'

        meiBlastCommand = 'blastn -task megablast -query ' + self.inFastFile \
                        + ' -db MEIdb  -best_hit_score_edge 0.1 -gapopen 5 -gapextend 2 -reward 1 -penalty -1 -word_size 11 -out ' + self.outMeiFile \
                        + ' -outfmt "6 qseqid sseqid qcovs stitle qlen qstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos sstrand'


        blast = os.system(blastCommand)
        meiBlast = os.system(meiBlastCommand)


class ApplyFilters:

    def __init__(self,inBlastFile): #inBlastFile,inMeiFile,outBEDFile):

        self.inBlastFile = inBlastFile
        #self.inMeiFile = inMeiFile
        #self.outBEDFile = outBEDFile
        self.filter()



    def filter(self):

        #blast_file = open(self.inBlastFile, 'r')

        blast_file_review = {}
        count = 0

        print "qseqid", "\t", "bitscore", "\t", "(length/qlen)", "\t", "pident", "\t", "gaps", "\t", "sstrand", "\t", "send", "\t", "sstart\n"
        for i in self.inBlastFile:
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
                processed_qseqid = int(qseqid.replace(chr + ":", "").replace("_", ""))

                bp_length = int(sstart) - int(processed_qseqid)
                blast_length = (length / qlen)
                mei_start = int(processed_qseqid) - 1

                Deletionfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
                Inversionfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
                TandomDupfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)

                if (bitscore > 50) and ((length / qlen) > 0.9) and (pident > 90) and (gaps < 3) and (bp_length > 50 and bp_length < 100000):

                    self.deletions(chr,Deletionfamily_strand_uniqueNumber,sstrand,send,sstart,
                                           processed_qseqid,bitscore)
                    self.meis(chr,mei_start,processed_qseqid,sstart,bitscore,qseqid,sstrand,count,send)
                    #print qseqid, "\t", bitscore, "\t", blast_length, "\t", pident, "\t", gaps, "\t", sstrand, "\t", send, "\t", sstart

                # inversions - less then 1kb in size
                if (bitscore > 50) and ((length / qlen) > 0.9) and (pident > 90) and (gaps < 3) and (bp_length < 100000):
                    self.inversions(chr, bitscore, send, sstart, Inversionfamily_strand_uniqueNumber, sstrand)

                if (bitscore > 50) and ((length / qlen) > 0.9) and (pident > 90) and (gaps < 3) and (bp_length > 50):
                    self.tandumDup(chr,bitscore,send,sstart,sstrand,qseqid,TandomDupfamily_strand_uniqueNumber)

    def deletions(self,chr,Deletionfamily_strand_uniqueNumber,sstrand,send,sstart,processed_qseqid,bitscore):
        #require that alignment is on plus strand nad downstream from org site
        #if send > sstart then plus, else minus
        #if sstart > qseqid - sequence aligns downstream
        deletionsFile = open('deletionsBED.txt','a+')
        if (sstrand.strip().rstrip() == 'plus') and (int(send) > int(sstart)) and (sstart > processed_qseqid):
            deletionsFile.write(chr + "\t" +  str(sstart) +  "\t" +  str(send)  +  "\t" + Deletionfamily_strand_uniqueNumber + "\t" + str(bitscore) + "\t" + "deletion" + "\n")
        deletionsFile.close()

    def meis(self,chr,mei_start,processed_qseqid,sstart,bitscore,qseqid,sstrand,count,send):
        #stop is org position of inserted sequence that is in query id
        #start is stop -1
        #if strand send > sstart then plus else minus

        meisFile = open('meisBED.txt', 'a+')
        if (send > sstart):
            sstrand="plus"
            MEIfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
            meisFile.write(chr + "\t" + str(mei_start) + "\t" + str(processed_qseqid) + "\t" + MEIfamily_strand_uniqueNumber + "\t" + str(bitscore) + "\t" + "MEI" + "\n")
        else:
            sstrand = "minus"
            MEIfamily_strand_uniqueNumber = qseqid + "_" + sstrand + "_" + str(count)
            meisFile.write(chr + "\t" + str(mei_start) + "\t" + str(processed_qseqid) + "\t" + MEIfamily_strand_uniqueNumber + "\t" + str(bitscore) + "\t" + "MEI" + "\n")
        meisFile.close()

    def inversions(self,chr,bitscore,send,sstart,Inversionfamily_strand_uniqueNumber,sstrand):
        inversionsFile = open('inversionsBED.txt','a+')
        #strand is minus
        if (sstrand.strip().rstrip() == 'minus'):
            inversionsFile.write(chr +  "\t" +  str(sstart) +  "\t" +  str(send) +  "\t" +  Inversionfamily_strand_uniqueNumber +  "\t" +  str(bitscore) +  "\t" +  "Inversions" + "\n")
        inversionsFile.close()
    def tandumDup(self,chr,bitscore,send,sstart,sstrand,qseqid,TandomDupfamily_strand_uniqueNumber):
        # sstart < qseqid
        # plus strand
        tandumDupFile = open('tandumDupBED.txt', 'a+')
        if (sstrand.strip().rstrip() == 'plus') and (sstart < qseqid):
            tandumDupFile.write(chr +  "\t" +  str(sstart) +  "\t" +  str(send) +  "\t" +  TandomDupfamily_strand_uniqueNumber +  "\t" +  str(bitscore) +  "\t" +  "Tandom Dups" + "\n")
        tandumDupFile.close()






def main():

    #inputfile = ''
    #outputfile = ''
    #minInsertLength = 25
    # chromosome='21'


    parser = argparse.ArgumentParser(prog='vcfToFasta.py',
                                     description='''Convert a vcf file to a fasta file''')

    parser.add_argument('-inFile', '--inputFile', type=str, required=True,
                        help='''VCF file to be parsed''')

    parser.add_argument('-outFile', '--outputFile', type=str, required=False,
                        help='''Name of fasta file that gets outputted''')

    parser.add_argument('-chr', '--chromosome', type=str, required=True,
                        help='''Chromosome''')

    parser.add_argument('-v', '--version', action='version', version='0.0.1 ')

    args = parser.parse_args()

    inVCFFile = args.inputFile
    outFastaFile = args.outputFile
    chr = args.chromosome
    minInsertLength = 25
    outBlastFile = chr + 'Blast.txt'
    outMeiFile = chr + 'BlastMEI.txt'
    outBEDFile = chr + 'BedFile.bed'


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
    print "here"
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
    #MakeFastaFile(inVCFFile,outFastaFile,chr,minInsertLength,logger)

    # take newly created fasta file and create a blast file and a mei blast file
    inFastFile = outFastaFile
    #MakeBlastFile(inFastFile,outBlastFile,outMeiFile)

    #take newly created blast file and mei blast file and apply filters to create bed file
    inBlastFile = outBlastFile
    inMeiFile = outMeiFile

    blast_file = open('blastchr21_output6.txt', 'r')
    #ApplyFilters(inBlastFile,inMeiFile,outBEDFile)
    ApplyFilters(blast_file) #inBlastFile, inMeiFile, outBEDFile)

if __name__ == "__main__":
    main()