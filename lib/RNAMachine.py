#!/usr/bin/env python

import os
import rna
import pdb
#from config import THRESHOLD, END_THRESHOLD, MAX_RNA_LEN, MIN_RNA_LEN, NEW_RNA_OBJECT_THRESHOLD

 
class ReadsWrapper:
    def __init__(self,data):
        self.data=data
        
    def __iter__(self):
        return self

    def push(self,element):
        self.data.insert(0,element)

    def has_data(self):
        return self.data

    def next(self):
        if not self.data:
            raise StopIteration
        element = self.data[0]
        self.data = self.data[1:]
        return element



class RNA_from_reads:
    WORKING = 0
    DONE    = 1
    #def __init__(self,reads,threshold=NEW_RNA_OBJECT_THRESHOLD,name=""):
    def __init__(self,reads,NEW_RNA_OBJECT_THRESHOLD,name=""):
        self.reads=reads
        self.status = RNA_from_reads.WORKING
        self.running_start = None
        self.start = None
        self.running_end = None
        self.end   = None
        self.count = 0
        self.new   = True
        self.threshold = NEW_RNA_OBJECT_THRESHOLD #threshold
        self.assigned_reads = []
        self.seq = ""
        self.consume()
        self.name=name

    def __str__(self):
        return "name='%s', start=%d, end=%d, seq=%s, count=%d" % (self.name,self.start,self.end,\
        self.seq, self.count)

    def __repr__(self):
        return "RNA(%s)" % self.__str__()


    def __contains__(self,item):
        return  item >= self.start and item <= self.end


    def consume(self):
        while self.status != RNA_from_reads.DONE:
            read = self.reads.next()
            seq , start, end = read[0]
            count = read[1]
            
            ## when the read starts and/or ends too far from the other reads
            ## processed, do not add it to the read pack, but push it
            ## back into the stack so it will be re-evaluated for other packs
            if self.start and abs(start - self.running_start) > self.threshold \
               and self.end and end > (self.running_end + self.threshold):
               #?and self.end and abs(end - self.running_end) > self.threshold:?
                self.status = RNA_from_reads.DONE
                self.reads.push(read)
                break
            
            ## start-a-new-pack signal
            if not self.start or start < self.start:
                self.start = start
                self.running_start = self.start
                
            if not self.end or end > self.end:
                self.end   = end
                self.running_end = self.end

            ## assign pack
            self.count += count
            self.running_start = start
            self.running_end = end
            self.assigned_reads.append(read)
            self.seq = seq
            if not self.reads.has_data():
                self.status = RNA_from_reads.DONE
                break
            

    
def rna_generator(blob, NEW_RNA_OBJECT_THRESHOLD):
    reads = blob.data.items()
    ## sort reads by start and then by end position
    reads = ReadsWrapper(sorted(reads, key=lambda x:(x[0][1], x[0][2])) )
    while reads.has_data():
        new_rna = RNA_from_reads(reads, NEW_RNA_OBJECT_THRESHOLD)
        yield new_rna
        
def rna_objects_from_blob(blob, THRESHOLD, NEW_RNA_OBJECT_THRESHOLD):
    return [rna_ for rna_ in rna_generator(blob, NEW_RNA_OBJECT_THRESHOLD) if rna_.count > THRESHOLD]
