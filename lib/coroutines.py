#!/usr/bin/env python

from __future__ import print_function
import pdb
import sys
import HTSeq
from time import time, sleep
from collections import defaultdict
from itertools import groupby
#import cPickle as pickle
import pickle




def broadcast(*targets):
    targets = [callable(target) and target() or target for target in targets]
    while True:
        item = (yield)
        for target in targets:
            target.send(item)

class Pump(object):
    pass

class Transformer(object):
    def set_target(self,target):
        self.target = target
    
    def __call__(self):
        while True:
            item = (yield)
            self.target.send(item)


class Filter(object):
    def __init__(self):
        self.queue = []

    def set_target(self,target):
        self.target = target

    def enqueue(self,*args):
        self.queue.append(args)

    def run_queue(self):
        for job in self.queue:
            name = job[0]
            args = job[1:]
            method = getattr(self,name)
            method(*args)


class Writer(object):
    def __init__(self,wid):
        self.n = 0
        self.wid = wid
        
    def __call__(self):
        # pdb.set_trace()
        while True:
            item = yield
            self.n += 1
            print("writer {}: {}".format(self.wid,item))
        

    
class Init(Pump):
    def run(self):
        for i in [1,1,2,2,3,3,3,3,3,3,3,1,1,4,4,1]:
            self.send(i)

class Tail(Transformer):
    def __call__(self):
        while True:
            item = ( yield )
            print(item)
    

class Meter(object):
    def __init__(self):
        self.n = 0
        self.target = None

    def __call__(self):
        start = time()
        old = 0
        frmstr = '\r{:<'+str(getTerminalSize()[0])+'}'
        item = yield
        while True:
            self.n += 1
            now = time()
            if  now - old > 1:
                delta = int(round(now-start))
                if not delta:
                    delta = 1
                speed = self.n / delta
                message =  "mean speed: %d reads/s" % speed
                sys.stderr.write(frmstr.format(message))
                sys.stderr.flush()
                old = now
            item = yield item


def getTerminalSize():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl, termios, struct, os
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (env['LINES'], env['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])

class FastqPump(Pump):
    def __init__(self,filename):
        self.fastq = iter(HTSeq.FastqReader(filename,"solexa")) 
    
    def __call__(self):
        for read in self.fastq:
            yield read

class SamPump(Pump):

    def __init__(self, filename, filter_unmapped = True):
        self.filter_unmapped = filter_unmapped
        sam_file = HTSeq.SAM_Reader(filename)
        self.sam = groupby(iter(sam_file), lambda x: x.read.name)

    def __call__(self):
        for aln in self.sam:
            if self.filter_unmapped:
                filt_aln = (aln[0], list(filter(lambda x: x.aligned, aln[1])))
                if len(filt_aln[1]) > 0:
                    yield filt_aln
            else:
                yield aln

class CountFilter(Filter):
    def __init__(self,threshold):
        Filter.__init__(self)
        self.threshold = threshold
        self.counter   = defaultdict(int)
        self.on_hold   = defaultdict(list)
        self._keeped = 0

    def preprocess(self,item):
        return item

    def key(self,item):
        return item

    def keeped(self):
        return self._keeped

    def rejected(self):
        return sum(map(lambda x: len(x[1]),self.on_hold.items()))

    def __call__(self):
        item = yield
        while True:
            item = self.preprocess(item)
            k    = self.key(item)
            if self.counter[ k ] >= self.threshold:
                self._keeped += 1
                item = yield item
            else:
                self.counter[k] += 1
                self.on_hold[k].append(item)
                item = yield

    def sync(self):
        for k, v in self.on_hold.items():
            if self.counter[k] >= self.threshold and v:
                for e in v:
                    self._keeped += 1
                    yield e
                self.on_hold[k] = []
 
class SamCountFilter(CountFilter):
    def __init__(self,threshold):
        CountFilter.__init__(self,threshold)

    def preprocess(self,item):
        name, alnmts = item
        alnmts  = list(alnmts)
        return (name, alnmts)

    def key(self,item):
        name, alnmts = item
        return alnmts[0].read.seq

class MultipleHitsGenomicFilter(Filter):
    def __init__(self,threshold,filename):
        Filter.__init__(self)
        self.threshold = threshold
        self.discarded = []
        with open(filename, 'rb') as f:
            self.genomic_hits = pickle.load(f)
        self.dist = defaultdict(int)

    def __call__(self):
        item = yield
        while True:
            name, alnmts = item
            alnmts = list(alnmts)
            n_in_hairpins = len([ alnmt for alnmt in alnmts if alnmt.iv != None and alnmt.iv.strand == "+" ])
            first_alnmt = alnmts[0]
            if first_alnmt.aligned: 
                seq = first_alnmt.read.seq
                n_in_genomic  = self.genomic_hits.get(seq,0)
                delta = n_in_genomic - n_in_hairpins
                if delta <= self.threshold:
                    self.dist[delta]+=1
                    item = yield item
                        #print item
                else:
                    self.discarded.append(delta)
                    item = yield
            else:
                item = yield

    def write_summary(self,filename):
        class Summary():
            def __init__(self,lenght,d_,dist):
                self.lenght = lenght
                self.d_ = d_
                self.dist = dist
                
            def write_to_file(self,filename):
                with open(filename,"w") as f:
                    f.write("total discarded because of multiple genomic mapping: %d\n" % self.lenght)
                    f.write("distribution:\n")
                    keys = sorted(self.d_.keys())
                    for key in keys:
                        f.write("%d:%d\n" % (key,self.d_[key]))
                    f.write("not discarded distribution:\n")
                    keys = sorted(self.dist.keys())
                    for key in keys:
                        f.write("%d:%d\n" % (key,self.dist[key]))
        d_ = {}
        discarded = sorted(self.discarded)
        for delta, l in groupby(discarded):
            d_[delta] = len(list(l))
        summary = Summary(len(self.discarded),d_,self.dist)
        summary.write_to_file(filename)
        #return Summary(len(self.discarded),d_,self.dist)
 
 
class MultipleHitsGenomicFilterVerbose(MultipleHitsGenomicFilter):
    def __init__(self,*args,**kwargs):
        super(MultipleHitsGenomicFilterVerbose,self).__init__(*args,**kwargs)
        self.aligned_to_genome = 0
        self.aligned_to_hairpins_exact = 0

    def __call__(self):
        while True:
            item = yield
            name, alnmts = item
            alnmts = list(alnmts)
            try:
                aligned_hairpins = [ alnmt for alnmt in alnmts if alnmt.iv != None and alnmt.iv.strand == "+" ]
                if [alnmt.optional_field("NM") == 0 for alnmt in aligned_hairpins]:
                    self.aligned_to_hairpins_exact += 1
                n_in_hairpins = len(aligned_hairpins)
            except AttributeError:
                pdb.set_trace()
            first_alnmt = alnmts[0]
            if first_alnmt.aligned: 
                seq = first_alnmt.read.seq
                n_in_genomic  = self.genomic_hits.get(seq,0)
                if n_in_genomic:
                    self.aligned_to_genome += 1
                delta = n_in_genomic - n_in_hairpins
                if delta <= self.threshold:
                    self.dist[delta]+=1
                else:
                    self.discarded.append(delta)

    def write_summary(self,filename):
        with open(filename,"w") as f:
            n_discarded = len(self.discarded)
            fraction = float(n_discarded) / self.aligned_to_hairpins_exact * 100
            f.write("%d reads were mapped to the human genome\n" % self.aligned_to_genome)
            f.write("%d were mapped to extended hairpins\n" % self.aligned_to_hairpins_exact)
            f.write("%d (%.1f%%) were discarded since mapping in more than %d genomic loci out of miRNA hairpins.\n" % (n_discarded,fraction,self.threshold) )
            f.write("... obtaining a final set of %d reads WOW\n" % (self.aligned_to_hairpins_exact - n_discarded))

