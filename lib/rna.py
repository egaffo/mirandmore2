#!/usr/bin/env python
from collections import defaultdict
import pdb
import os
import cPickle as pickle
from Bio import pairwise2
#import stream
from functools import partial
import re
import sys
#try:
#    from config import EXACT, SHORTER_OR_LONGER, MIS_1, MIS_2, MIN_COUNT,\
#    SISTER_OVERHANG_LEN, SISTER_MM_PROP_THRESHOLD, MIN_MORNA_LEN,\
#    SISTER_MATCHES_THRESHOLD, SPECIES, ALLOWED_OVERHANG
#except ImportError:
#    from deep_config import EXACT, SHORTER_OR_LONGER, MIS_1, MIS_2, MIN_COUNT,\
#    SISTER_OVERHANG_LEN, SISTER_MM_PROP_THRESHOLD, MIN_MORNA_LEN,\
#    SISTER_MATCHES_THRESHOLD, SPECIES, ALLOWED_OVERHANG


#MIRANDMORE_HOME = os.environ["MIRANDMORE_HOME"]



def give_me_0():
    return 0

def build_pre_to_mature_table(mature_table):
    table = defaultdict(list)
    #mirna_annotations = open(os.path.join(MIRANDMORE_HOME,"annotations","mirbase",mature_table),"r")
    mirna_annotations = open(mature_table,"r")
    for line in mirna_annotations:
        _,pre,name,_,start,end,pre_len = line.strip().split("\t")
        if int(pre_len) - int(end) < int(start):
            order = "3p"
        else:
            order = "5p"

        mature=Mature(pre,name,start,end,order)
        table[pre].append(mature)
    return table

def build_pre_to_mature_map(mature_table):
    table = defaultdict(dict)
    #mirna_annotations = open(os.path.join(MIRANDMORE_HOME,"annotations","mirbase",mature_table),"r")
    mirna_annotations = open(mature_table,"r")
    for line in mirna_annotations:
        _,pre,name,_,start,end,pre_len = line.strip().split("\t")
        if int(pre_len) - int(end) < int(start):
            order = "3p"
        else:
            order = "5p"

        mature=Mature(pre,name,start,end,order)
        table[pre][mature.name]=mature
    return table


def build_mature_to_pre_table(pre_to_mature_table):
    d={}
    for key,values in pre_to_mature_table.items():
        for value in values:
            d[value.name]=key
    return d

class Mature:
    def __init__(self,pre,name,start,end,order):
        self.name = name
        self.start = int(start)
        self.end = int(end)
        self.pre = pre
        self.order = order

    def __repr__(self):
        return "Mature %s: from %d to %d in pre %s (%s) " % (self.name,self.start,self.end,self.pre,self.order)



class MatureResult:
    def __init__(self):
        self.exact = defaultdict(int)
        self.shorter_or_longer = defaultdict(int)
        self.mismatch_1 = defaultdict(int)
        self.mismatch_2 = defaultdict(int)
        self.categories = [self.exact,
                           self.shorter_or_longer,
                           self.mismatch_1,
                           self.mismatch_2]
        
    def add(self,entry):
        category, mature, sequence , pre, start, end, variant_end, weight = entry
        self.categories[category][sequence] += weight

class MatureResultInPre(MatureResult):
    def __init__(self):
        MatureResult.__init__(self)

    def add(self,entry):
        category, mature, sequence , pre, start, end, variant_end, weight = entry
        self.categories[category][(sequence,start,end,variant_end)] += weight

        
class MatureResultSet:
    def __init__(self):
        self.data = defaultdict(MatureResult)

    def __setitem__(self,key,value):
        self.data[key]=value

    def __getitem__(self,key):
        return self.data[key]

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def add(self,entry):
        # pdb.set_trace()
        mature = entry[1]
        self.data[mature].add(entry)



class PreResult:
    def __init__(self):
        self.data = defaultdict(MatureResultInPre)

    def add(self, entry):
        category, mature, sequence , pre, start, end, variant_end, weight = entry
        self.data[mature].add(entry)

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def __getitem__(self,key):
        return self.data[key]

    def __contains__(self,name):
        return name in self.data


    
class PreResultSet:
    def __init__(self):
        self.data = defaultdict(PreResult)

    def __setitem__(self,key,value):
        self.data[key]=value

    def __getitem__(self,key):
        return self.data[key]

    def keys(self):
        return self.data.keys()

    def items(self):
        return self.data.items()

    def add(self,entry):
        # pdb.set_trace()
        pre = entry[3]
        self.data[pre].add(entry)
    

class PreExactResult:
    def __init__(self):
        self.data = defaultdict(int)

    def add(self,entry):
        sequence, pre_name, start, end, weight = entry
        self.data[(sequence,start,end)] += weight
        
    def items(self):
        return self.data.items()

    def __setitem__(self,key):
        self.data[key]=value

    def __getitem__(self,key):
        return self.data[key]
    
class PreExactResultSet:
    def __init__(self):
        self.data = defaultdict(PreExactResult)
    
    def __setitem__(self,key,value):
        self.data[key]=value

    def __getitem__(self,key):
        return self.data[key]

    def keys(self):
        return self.data.keys()

    def values(self):
        return self.data.values()

    def items(self):
        return self.data.items()

    def add(self,entry):
        pre = entry[1]
        self.data[pre].add(entry)

class PreAnnotation(object):
    def __init__(self,name,seq,structure,energy):
        self.name = name
        self.seq  = seq
        self.structure = structure
        self.energy = energy

class HairpinProxy(object):
    #def __init__(self,blob="hairpin." + SPECIES + ".extended.blob"):
    def __init__(self, HAIRPIN_EXTENDED_BLOB):
        #in_file = os.path.join(MIRANDMORE_HOME,"annotations","mirbase",blob)
        in_file = HAIRPIN_EXTENDED_BLOB
        with open(in_file,"r") as f:
            self.data = pickle.load(f)

    def getDNA(self,hairpin,start,end):
        seq = self.data[hairpin]
        return seq[start -1: end]

    def lookupRNA(self,hairpin,rna_):
        seq = self.data[hairpin]
        return seq[rna_.start - 1: rna_.end]

WT = {'A':'T',
      'T':'A',
      'C':'G',
      'G':'C'}

def reverse_complement(seq):
    seq = seq[::-1]
    seq = map(lambda x: WT[x.upper()],seq)
    return "".join(seq)


def rnafold_to_dict(structure):
    stack = []
    tuples = []
    i = 1
    for element in structure:
        if element == '(':
            stack.append(i)
        elif element == '.':
            pass
        elif element == ')':
            tuples.append((stack.pop(),i))
        i+=1
    d_ = dict(tuples)
    tuples_reverse = [(b,a) for (a,b) in tuples]
    d_.update(tuples_reverse)
    return d_


def mature_couple_score(m1,m2,match_dict, SISTER_OVERHANG_LEN):
    score = 0
    overhang_score = 0
    mismatches = 0
    m1_len = m1.end - m1.start + 1
    for i in xrange(m1.start,m1.end+1):
        match_pos =  match_dict.get(i)
        if match_pos and match_pos >= m2.start and match_pos <= m2.end:   #and in_m2(match_pos):
            score += 1
        elif i > m1.end - SISTER_OVERHANG_LEN and (not match_pos or (match_pos < m2.start or match_pos > m2.end)):
            overhang_score += 1
        else:
            mismatches += 1
            
    return score, mismatches, overhang_score,  m1_len,  m1.order


    
def hairpin_summary(name, MATURE_TABLE, HAIRPIN_ANNOTATION_BLOB):
#    table = build_pre_to_mature_table("mature-table.txt")
#    with open(os.path.join(MIRANDMORE_HOME,"annotations","mirbase","hairpin." + SPECIES + ".annotations.blob")) as f:
#       annotations_table = pickle.load(f)
    table = build_pre_to_mature_table(MATURE_TABLE)
    with open(HAIRPIN_ANNOTATION_BLOB) as f:
       annotations_table = pickle.load(f)
    m1, m2  = table[name]
    match_dict = rnafold_to_dict(annotations_table[name].structure)
    m1_scores = mature_couple_score(m1,m2,match_dict)
    m2_scores = mature_couple_score(m2,m1,match_dict)
    return {m1.name:m1_scores,m2.name:m2_scores}
    
class Oracle:
    def __init__(self,mature_table,annotations):
        self.table = build_pre_to_mature_table(mature_table)
        #with open(os.path.join(MIRANDMORE_HOME,"annotations","mirbase",annotations)) as f:
        with open(annotations) as f:
            self.annotations_table = pickle.load(f)

    def is_mature_sister(self,rna_,mature, SISTER_OVERHANG_LEN, SISTER_MATCHES_THRESHOLD):
        score = 0
        overhang_score = 0
        mismatches = 0
        rna_len = rna_.end - rna_.start + 1
        match_dict = rnafold_to_dict(self.annotations_table[mature.pre].structure)
        for i in xrange(rna_.start,rna_.end+1):
            match_pos =  match_dict.get(i)
            if match_pos and match_pos >= mature.start and match_pos <= mature.end: #and in_mature(match_pos):
                score += 1
            elif i > rna_.end - SISTER_OVERHANG_LEN and (not match_pos or (match_pos < mature.start or match_pos > mature.end)):
                overhang_score += 1
            else:
                mismatches += 1
        #mm_prop = float(mismatches)/rna_len #never used here
        return (score >= SISTER_MATCHES_THRESHOLD)


    def score_sisters(self,rna_,mature, SISTER_OVERHANG_LEN):
        score = 0
        overhang_score = 0
        mismatches = 0
        rna_len = rna_.end - rna_.start + 1
        match_dict = rnafold_to_dict(self.annotations_table[mature.pre].structure)
        for i in xrange(rna_.start,rna_.end+1):
            match_pos =  match_dict.get(i)
            if match_pos and match_pos >= mature.start and match_pos <= mature.end:   #and in_mature(match_pos):
                score += 1
            elif i > rna_.end - SISTER_OVERHANG_LEN and (not match_pos or (match_pos < mature.start or match_pos > mature.end)):
                overhang_score += 1
            else:
                mismatches += 1
        mm_prop = float(mismatches)/rna_len
        return {"mismatches":mismatches,
                "matches":score,
                "overhang_score":overhang_score,
                "mm_prop":mm_prop}

    def get_shadow(self,mature):
        match_dict = rnafold_to_dict(self.annotations_table[mature.pre].structure)
        matching_positions =  map(lambda x: match_dict.get(x),xrange(mature.start,mature.end + 1))
        matching_positions =  filter(lambda x: not x == None,matching_positions)

        ## WARNING: it might happens that there is no matching positions in the folding
        ## TODO: raise some an exception here/handle the case.
        ## The current workaround is to give fake-nonsense negative matching positions
        if not matching_positions:
            matching_positions = [-1]
            ## TODO: set a proper logging system
            sys.stderr.write('WARNING: no folding matches for mature ' + str(mature))
        return (min(matching_positions),max(matching_positions))
            
    
def rna_len(rna_):
    return rna_.end - rna_.start + 1

def midpoint(obj):
    return (obj.start + obj.end)/2

class TwoMatureAssignmentError(BaseException):
    pass

class AssignmentError(BaseException):
    pass

class NamingError(BaseException):
    pass

def mir_suffix(name, SPECIES):
    name = re.sub(SPECIES + "-(mir-)?", "", name)
    return name

def normalize_mir_name(name):
    return name.replace("-ext", "")


def name_of_sister(pre_name, mature):
    """Format and assign the name of the sister miRNA.
    
    Parameters
    ----------

    pre_name : string
        The name of the precursor.
    mature: an object generated by the build_pre_to_mature_table() function
        The mature miRNA on the opposite arm of the hairpin.

    Returns
    -------
    sister_miRNA_name : string
        The formatted name for the sister miRNA.

    Raises
    ------
    NamingError
        In case it is not possible to infer the 3p/5p suffix for the sister miRNA.
        
    """

    pre_name = pre_name.replace("r","R").replace("-ext","")
    if mature.name.endswith("-5p") or mature.order == "5p":
        return pre_name + "-3p"
    elif mature.name.endswith("-3p") or mature.order == "3p":
        return pre_name + "-5p"
    else:
        raise NamingError("Cannot assign a suffix to {} sister miRNA".format(mature.name))


def name_of_mor(pre_name, suffix):
    """Format and assign the name of the moRNA.

    Parameters
    ----------
    pre_name: string
        The name of the precursor.
    suffix: string
        Either '-5p' or '-3p'. This must be inferred prior to this function

    Returns
    -------
    moRNA_name : string
        The formatted name for the moRNA.
    """

    sp_pre_name = pre_name.split("-")
    if "mir" in sp_pre_name:
        pre_name = pre_name.replace("mir","moR").replace("-ext","")
    elif "let" in sp_pre_name:
        pre_name = pre_name.replace("let","moR-let").replace("-ext","")
    else:
        pre_name = pre_name.replace("-ext","") + "-moR"
    return pre_name + suffix


def name_of_loop(pre_name):
    """Format the name of the loop.

    Parameters
    ----------
    pre_name: string
        The name of the precursor

    Returns
    -------
    loop_name : string
        The formatted name of the loop
    """

    sp_pre_name = pre_name.split("-")
    if "mir" in sp_pre_name:
        pre_name = pre_name.replace("mir","loop").replace("-ext","")
    elif "let" in sp_pre_name:
        pre_name = pre_name.replace("let","loop-let").replace("-ext","")
    else:
        pre_name = pre_name.replace("-ext","") + "-loop"
    return pre_name


class MockRna:
    def __init__(self,name="",count=0,start=0,end=0,status=0):
        self.name = name
        self.count = count
        self.start = start
        self.end  = end
        self.status = status


class PreSummary:
    ass_methods = ["try_as_known_mir", "try_as_sister_mir", "try_as_loop", "try_as_mor"]
    
    def __init__(self, name, mature_table, oracle, MIN_MORNA_LEN, SPECIES,
                 SISTER_OVERHANG_LEN, SISTER_MATCHES_THRESHOLD):
        self.name = name
        self.stdName = normalize_mir_name(name)
        self.suffix = mir_suffix(self.stdName, SPECIES)
        self.matures = mature_table[name]
        self.n_matures = len(self.matures)
        self.five_prime_mor = None
        self.five_prime_mir = None
        self.loop = None
        self.three_prime_mir = None
        self.three_prime_mor = None
        self.oracle = oracle 
        self.stack = []
        self.data = {}
        self.MIN_MORNA_LEN = MIN_MORNA_LEN
        self.SPECIES = SPECIES
        self.SISTER_OVERHANG_LEN = SISTER_OVERHANG_LEN
        self.SISTER_MATCHES_THRESHOLD = SISTER_MATCHES_THRESHOLD


    def __repr__(self):
        return str(self)


    def __str__(self):
        return """
        PreSummary(
        mor-5p='%s',
        mir-5p='%s',
        loop='%s',
        mir-3p='%s',
        mor-3p='%s') 
        """ % (str(self.five_prime_mor), str(self.five_prime_mir), str(self.loop), \
               str(self.three_prime_mir), str(self.three_prime_mor))


    def any_mature_assigned(self):
        return self.five_prime_mir or self.three_prime_mir


    def all_matures_assigned(self):
        return self.five_prime_mir and self.three_prime_mir


    def try_as_known_mir(self,rna_):
        for mature in self.matures:
            mature_mid_point = midpoint(mature)
            if mature_mid_point > rna_.start and mature_mid_point < rna_.end: #\
            #and mature.start-ALLOWED_OVERHANG <= rna_.start and \
            #mature.end+ALLOWED_OVERHANG >= rna_.end:
                if mature.order == "5p":
                    ## prevent the read pack trying to get an already assigned slot
                    if not self.five_prime_mir:
                        self.five_prime_mir = rna_
                    else:
                        continue
                else:
                    if not self.three_prime_mir:
                        self.three_prime_mir = rna_
                    else:
                        continue
                rna_.name = mature.name.replace("-5p", "").replace("-3p", "") + '-' + mature.order
                rna_.new = False
                return True
        return False

    def try_as_sister_mir(self,rna_):
        if self.n_matures == 1:
            mature = self.matures[0]
            # obtain an oracle instance somehow somewhere
            if self.oracle.is_mature_sister(rna_, mature, 
                                            self.SISTER_OVERHANG_LEN, 
                                            self.SISTER_MATCHES_THRESHOLD):
                if mature.order == "5p":
                    self.three_prime_mir = rna_
                else:
                    self.five_prime_mir  = rna_
                try:
                    rna_.name = name_of_sister(self.name, mature)
                except NamingError, ne:
                    rna_.name = mature.name + "*"
                    ## TODO: set a proper logging system
                    sys.stderr.write("WARNING: {} may have a non standard name."\
                                     "Please check name of {}".format(rna_.name, mature.name))
                return True
        return False

    def try_as_loop(self,rna_):
		rna_midpoint = midpoint(rna_)
		if self.n_matures == 2:
	   	    if self.matures[0].order == "5p":	
			    end = self.matures[0].end
			    start = self.matures[1].start
		    else:
				end = self.matures[1].end
				start = self.matures[0].start
		elif self.n_matures == 1:
			if self.five_prime_mir:
				end = self.five_prime_mir.end
				if self.matures[0].order == "5p":
					return False		
				if self.matures[0].order == "3p":
				   start = self.matures[0].start
			elif self.three_prime_mir:
				start = self.three_prime_mir.start
				if self.matures[0].order == "3p":
					return False		
				if self.matures[0].order == "5p":
				   end = self.matures[0].end
			else:
				return False	
		else:
			raise AssignmentError("Can't assign '%s' in pre '%s'. "\
			"Trying sequence as loop failed because none mature miRNAs is known. "\
			"This should not happen!" % (rna_, self.name))
		if rna_midpoint > end and rna_midpoint < start:
		        self.loop = rna_
		        rna_.name = name_of_loop(self.name)
		        return True
		return False


    def try_as_mor(self,rna_):
        ## first condition: appropriate length
        if rna_len(rna_) > self.MIN_MORNA_LEN:

            ## check first with known miRs, even when they have no ERE assigned
            for mature in self.matures:
                mature_mid_point = midpoint(mature)
                if mature_mid_point > rna_.end and mature.order == "5p":
                    self.five_prime_mor = rna_
                    rna_.name = name_of_mor(self.name, '-5p')
                    return True
                if mature_mid_point < rna_.start and mature.order == "3p":
                    self.three_prime_mor = rna_
                    rna_.name = name_of_mor(self.name, '-3p')
                    return True
            
            ## check now the case of moRs adjacent to new sister miRs
            if self.any_mature_assigned():
                if self.five_prime_mir and rna_.end < midpoint(self.five_prime_mir):
                    self.five_prime_mor = rna_
                    rna_.name = name_of_mor(self.name, '-5p')
                    return True                    
                elif self.three_prime_mir and rna_.start > midpoint(self.three_prime_mir):
                    self.three_prime_mor = rna_
                    rna_.name = name_of_mor(self.name, '-3p')
                    return True
        return False


    def project_shadow(self):
        if self.n_matures == 2:
            raise TwoMatureAssignmentError("Two matures and still can't assign! WTF!")
        mature = self.matures[0]
        return self.oracle.get_shadow(mature)


    def inject_mock(self):
        start, end =  self.project_shadow()
        mock = MockRna(start=start,end=end)
        mature = self.matures[0]
        if mature.order == "5p":
            self.three_prime_mir = mock
        else:
            self.five_prime_mir  = mock


    def _assign(self,rna_):
                
        ## This code will apply every method defined in the PreSummary class
        ## to the rna_ object (that is a pack of reads/alignments) and 
        ## save the results to a list of booleans telling which methods
        ## were able to assign the read pack. In other words, this code
        ## will check if the read pack can be assigned to any element 
        ## on the miRNA precursor (i.e. 5' moR, 5' miR, loop, 3' miR and 
        ## 3' moR), but since it tries all precursor elements, it could 
        ## be assigned to more than one element (then the result list 
        ## will have > 1 True values).
        #assigned = []
        #for method in PreSummary.ass_methods:
        #    #assigned.append(getattr(self, method)(rna_))?
        #    if getattr(self, method)(rna_):
        #       assigned.append(rna_)
        ## check for ambiguous assignments
        ##if sum(assigned) > 1:
        ##    log.debug.write('''WARNING: read assigned to multiple elements '''
        ##              '''in precursor''')
        #return any(assigned)

        ## Different behaviour would be 'stop checking as soon as the 
        ## read pack is assigned' if using the inline any(), since the 
        ## any() function stops iterating at the first True value
        ## Doing so, (i) it is not possible to check possible multiple 
        ## assignments, and (ii) a possible assignment bias is introduced 
        ## in favour of 'left-sided' elements in the precursor 
        return any(getattr(self, method)(rna_) for method in PreSummary.ass_methods)
        
        ## old code using the stream package (deprecated)
        #methods = map(partial(getattr,self),PreSummary.ass_methods)
        #methods = map(lambda x:  partial(x,rna_=rna_),methods)
        #assigned = iter(methods) >> stream.dropwhile(lambda x: not x()) >> stream.take(1) >> list
        #return assigned
      

    def assign(self,rna_):
        assigned = self._assign(rna_)
        if not assigned:
            self.stack.append(rna_)

    def assign_from_stack(self,rna_):
        assigned = self._assign(rna_)
        if not assigned:
            if not (self.five_prime_mir and self.three_prime_mir):
                try:
		    self.inject_mock()
		except TwoMatureAssignmentError, tmae:
		    raise TwoMatureAssignmentError("Can't assign rna '%s' in pre '%s'. %s" % (rna_, self.name, tmae.message))
            if not(self._assign(rna_)):
                # pdb.set_trace()
                raise AssignmentError("Can't assign rna '%s' in pre '%s' " % (rna_,self.name))

    def clean(self):
        for name in ["five_prime_mir","three_prime_mir"]:
            mir = getattr(self,name)
            if isinstance(mir,MockRna):
                setattr(self,name,None)
                break

    #def as_dict(self):
    #    _dict = {}
    #    for attr in ("five_prime_mir", "five_prime_mor", "loop",
    #                 "three_prime_mir", "three_prime_mor"):
    #        obj = getattr(self,attr) 
    #        if obj  is not None: 
    #            _dict[obj.name] = obj
    #    return _dict

    def populate(self, rnas):
        for rna_ in rnas:
            self.assign(rna_)
        while self.stack:
            rna_ = self.stack.pop()
            self.assign_from_stack(rna_)
        
        self.clean()
        #self.data = self.as_dict()
