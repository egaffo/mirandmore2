from datetime import datetime
import urllib2
import time
import re
import SCons.Node
import os.path
from ftplib import FTP
from urlparse import urlparse


MDTM_OK = "213"
MDTM_ERROR = "550"


def basename(name):
    name = os.path.basename(name)
    return os.path.splitext(name)[0]

def mdtm_to_timestamp(mdtm):
    status, mdtm = mdtm.split()
    if status == MDTM_ERROR:
        raise ValueError("file not found")
    mdtm = re.sub("\.\d+","",mdtm)
    year =   int(mdtm[0:4])
    month =  int(mdtm[4:6])
    day   =  int(mdtm[6:8])
    hour  =  int(mdtm[8:10])
    minute = int(mdtm[10:12])
    second = int(mdtm[12:14])
    cal = datetime(year,month,day,hour,minute,second)
    return int(time.mktime(cal.timetuple()))



class DepValue(SCons.Node.Python.Value):
    def __init__(self,name):
        value = globals().get(name)
        SCons.Node.Python.Value.__init__(self,value)
        self.name = name

    def __str__(self):
        return self.name

    def changed_since_last_build(self,target,prev_ni):
        if not hasattr(prev_ni,"csig"):
            return True
        if prev_ni.csig == str(self.value):
            return False
        else:
            return True
        return True


class FakeUri(SCons.Node.Python.Value):
    def __init__(self,value):
        SCons.Node.Python.Value.__init__(self,value)
        self.value=value

    def __str__(self):
        return self.value
    
    def changed_since_last_build(self,target,prev_ni):
        return False

    def get_csig(self,calc=None):
        return "OxDEADBEEF"

class Uri(SCons.Node.Python.Value):
    '''
    A dependency type based on a remote URI. It supports HTTP and FTP.
    To check freshness uses etag header for HTTP and MDTM command for FTP
    '''
    def __init__(self,value):
       SCons.Node.Python.Value.__init__(self,value)
       self.value = value

    def str_for_display(self):
        return self.value

    def __str__(self):
        return self.value

    def changed_since_last_build(self,target,prev_ni):
        return getattr(self,self.get_protocol()+"_changed_since_last_build")(target,prev_ni)

    def ftp_changed_since_last_build(self,target,prev_ni):
        if not hasattr(prev_ni,"csig"):
            return True
        if int(prev_ni.csig) >= int(self.get_csig()):
            return False
        else:
            return True
        return True
        
    def http_changed_since_last_build(self,target,prev_ni):
        if not hasattr(prev_ni,"csig"):
            return True
        if prev_ni.csig == self.get_csig():
            return False
        else:
            return True
        return True

    def https_changed_since_last_build(self,target,prev_ni):
        return self.http_changed_since_last_build(target,prev_ni)

    def get_protocol(self):
        return self.value.split("://")[0]

    def get_http_sig(self):
        request = urllib2.Request(self.value)
        opener  = urllib2.build_opener()
        response = opener.open(request)
        try:
            return response.headers.get("etag").replace('"','')
        except AttributeError:
            pass

    def get_https_sig(self):
        return self.get_http_sig()

    def get_ftp_sig(self):
        parsed_url = urlparse(self.value)
        ftp=FTP(parsed_url.netloc)
        ftp.login()
        ftp.cwd(os.path.dirname(parsed_url.path))
        mdtm = ftp.sendcmd("MDTM %s" % os.path.basename(parsed_url.path))
        return str(mdtm_to_timestamp(mdtm))

    def get_csig(self,calc=None):
        try:
            return self.ninfo.csig
        except AttributeError:
            pass
        csig = getattr(self,"get_"+self.get_protocol()+"_sig")()
        self.get_ninfo().csig = csig
        return csig


def extension(filename):
    return filename.split(".")[-1]

def smart_decider(dependency,target,prev_ni):
    '''
    smart_decider looks at content to evaluate dependecy freshness if dependency is a python file.
    Otherwise it looks at timestamp.
    In such way we prevent costly rebuilds every time we switch git branch
    '''
    if extension(str(dependency)).lower() in ('fq','tfq'):
        return False

    if extension(str(dependency)).lower() in ('py','r'):
        cur_csig = dependency.get_csig()
        try:
            return cur_csig != prev_ni.csig
        except AttributeError:
            return True
    else:
        try:
            return dependency.get_timestamp() > target.get_timestamp() 
        except AttributeError:
            return True
