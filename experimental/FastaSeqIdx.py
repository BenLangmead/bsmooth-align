import os
import csv

class FastaSeqIdx(object):
    
    ''' Object encapsulating a FASTA file that has been indexed by sequence.
        That is, the index makes it easy to seek directly to the beginning of
        a specific sequence (but not to specific points within it). '''
    
    def __init__(self, fa_fns, truncate=True, idx_fn=None):
        self.fa_fns = fa_fns
        self.idx_fn = idx_fn
        self.idx = dict()
        self.lens = dict()
        self.truncate = truncate
        self.fh = None
        self._lastc = ''
        if idx_fn is not None and os.path.exists(idx_fn):
            self.loadIdx(idx_fn)
        else:
            self.buildIdx()
        if idx_fn is not None and not os.path.exists(idx_fn):
            self.saveIdx(idx_fn)
    
    def __del__(self):
        if self.fh is not None: self.fh.close()
    
    def close(self):
        if self.fh is not None: self.fh.close()
    
    def nextc(self, n=1, skip=False):
        ''' Return the next string of n characters. '''
        self._lastc = ''
        if n == 1:
            c = self.fh.read(1)
            while c.isspace():
                c = self.fh.read(1)
            self._lastc = c
            return c # Empty string if EOF
        else:
            s = bytearray()
            slen = 0
            lastc = None
            while slen < n:
                c = self.fh.read(1)
                while c.isspace():
                    c = self.fh.read(1)
                if not skip:
                    s += c
                lastc = c
                slen += 1
            if slen > 0 and lastc is not None:
                self._lastc = lastc
            return str(s)
    
    def get(self, ref, off=0, ln=None):
        ''' Get the string at 0-based offset off, with length ln '''
        if off == 0 and ln is None:
            return self.getWhole(ref)
        if ln is None: ln = self.lens[ref]
        self.scanTo(ref)
        self.nextc(n=off, skip=True)
        return self.nextc(n=ln, skip=False)
    
    def getpad(self, ref, off=0, ln=None, pad='x'):
        ''' Get the string at 0-based offset off, with length ln.  If this
            string falls off either end of the reference, fill those slots with
            the "pad" string. '''
        if off == 0 and ln is None:
            return self.getWhole(ref)
        if ln is None: ln = self.lens[ref]
        lf = ''
        rt = '' 
        if off < 0: lf = pad * -off
        if off+ln > self.lens[ref]: rt = pad * (off+ln - self.lens[ref])
        self.scanTo(ref)
        self.nextc(n=off, skip=True)
        return lf + self.nextc(n=ln, skip=False) + rt
    
    def getWhole(self, ref):
        self.scanTo(ref)
        lns = []
        while True:
            ln = self.fh.readline().rstrip()
            if len(ln) == 0 or ln[0] == '>': break
            lns.append(ln)
        return "".join(lns)
    
    def peek(self):
        ''' Take a peek at the next character and return it '''
        off = self.fh.tell()
        c = self.fh.read(1)
        while c.isspace():
            c = self.fh.read(1)
        self.fh.seek(off)
        return c
    
    def lastc(self):
        ''' Return last character read '''
        return self._lastc
    
    def scanTo(self, ref):
        ''' Scan to the beginning of the FASTA sequence of the given name '''
        if self.fh is not None:
            self.fh.close()
            self.fh = None
        if self.truncate:
            idx = ref.find(' ')
            if idx >= 0: ref = ref[:idx]
        assert ref in self.idx
        fn, off = self.idx[ref]
        self.fh = open(fn, 'rb')
        self.fh.seek(off, 0)
    
    def loadIdx(self, idx_fn):
        ''' Load an already-calculated index from a file. '''
        csvrd = csv.reader(open(idx_fn, 'rb'), delimiter='\t')
        self.idx = dict()
        self.lens = dict()
        for row in csvrd:
            ref, fn, off, ln = row
            self.idx[ref] = (fn, int(off))
            self.lens[ref] = int(ln)
    
    def saveIdx(self, idx_fn):
        ''' Save an index to a file. '''
        csvw = csv.writer(open(idx_fn, 'wb'), delimiter='\t')
        for ref, v in self.idx.iteritems():
            fn, off = v
            csvw.writerow([ref, fn, off, self.lens[ref]])
    
    def length(self, nm):
        assert nm in self.lens
        return self.lens[nm]
    
    def ref_key_iter(self):
        return self.lens.iterkeys()
    
    def ref_item_iter(self):
        return self.lens.iteritems()
    
    def buildIdx(self):
        ''' Build an index dictionary over for all the reference sequences in
            all the FASTA files. '''
        for fn in self.fa_fns:
            fh = open(fn, 'rb')
            cur = None
            while 1:
                line = fh.readline()
                if not line:
                    break
                line = line.rstrip()
                if len(line) == 0: continue
                off = fh.tell()
                if line[0] == '>': # Header
                    ref = line[1:]
                    if self.truncate:
                        idx = ref.find(' ')
                        if idx >= 0: ref = ref[:idx]
                    assert ref not in self.idx, "Already saw ref '%s' from line '%s'" % (ref, line)
                    self.idx[ref] = (fn, off)
                    cur = ref
                    self.lens[cur] = 0
                elif cur is not None:
                    self.lens[cur] += len(line)
            fh.close()

if __name__ == '__main__':
    import unittest
    unittest.main()
