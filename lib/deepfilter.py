import HTSeq
import stream
#from config import QUALITY_ENCODING
import sys

def dsopen(name,quality_encoding):
    return HTSeq.FastqReader(name,quality_encoding)

def mean(sequence):
    return(sum(sequence)/float(len(sequence)))

def trim(read):
    return HTSeq.SequenceWithQualities(read.seq[:-1],read.name,read.qualstr[:-1])


class QualityFilter():
    def __init__(self,treshold):
        self.treshold=treshold
        self.message="quality < %d" % treshold
        self.discarded=0

    def __call__(self,read):
        if all(read.qual >= self.treshold):
            return True
        self.discarded+=1
        return False
            

class MaxLenFilter():
    def __init__(self, maxlen, filehandler = None):
        self.maxlen=maxlen
        self.message="len > %d" % maxlen
        self.discarded=0
        self.long_reads_file = filehandler

    def __call__(self,read):
        if len(read) <= self.maxlen:
            return True
        if self.long_reads_file:
            read.write_to_fastq_file(self.long_reads_file)
        self.discarded+=1
        return False

class MeanQualityFilter():
    def __init__(self,treshold):
        self.treshold=treshold
        self.message="mean quality < %d" % treshold
        self.discarded=0

    def __call__(self,read):
        if mean(read.qual) >= self.treshold:
            return True
        self.discarded+=1
        return False
            
class MorePositionsLessThan():
    def __init__(self,n,treshold):
        self.n=n
        self.treshold=treshold
        self.message="more than %d bases with quality < %d" % (n,treshold)
        self.discarded=0

    def __call__(self,read):
        if sum(read.qual < self.treshold) < self.n:
            return True
        self.discarded+=1
        return False

class Counter(stream.Stream):
    def __init__(self):
        self.n = 0

    def __call__(self,iterator):
        def counter():
            for el in iterator:
                self.n += 1
                yield el
        return counter()

    def __pipe__(self,inpipe):
        self.iterator = self.__call__(iter(inpipe))
        return self


class FastqWriter(stream.Stream):
    def __init__(self,path):
        if path == "-":
            self.out = sys.stdout
        else:
            self.out = open(path,"w")

    def __pipe__(self,inpipe):
        iterator = iter(inpipe)
        for read in iterator:
            read.write_to_fastq_file(self.out)
        self.out.close()
