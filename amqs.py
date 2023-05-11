import numpy as np
from bloom_filter import BloomFilter
import sys
import time
import bbhash
import hashlib
import wyhash
from pympler import asizeof
import xxhash
import mmh3

class Bloom:
    def __init__(self, keys, error_rate):
        self.keys = keys
        self.error_rate = error_rate
        self.bloom = None
    def build_bloom_filter(self):
        bloom = BloomFilter(max_elements = len(self.keys), error_rate = self.error_rate)
        for k in self.keys:
            bloom.add(k)
        self.bloom = bloom
    def query_bf(self, query_keys): 
        bf_quer = []
        for k in query_keys:
            bf_quer.append(k in self.bloom)
        return bf_quer
    def compute_fp_rate(self, guess, query_keys, tn):
        fp = 0
        fn = 0
        for i in range(len(query_keys)):
            in_actual = query_keys[i] in self.keys
            in_bf = guess[i]
            if in_actual and not in_bf:
                fn+=1
            elif in_bf and not in_actual:
                fp+=1
            else:
                continue
        fp_rate = fp/(fp+tn)
        return fp_rate
    def compute_fn_rate(self, guess, query_keys, tp):
        fp = 0
        fn = 0
        for i in range(len(query_keys)):
            in_actual = query_keys[i] in self.keys
            in_bb = guess[i]
            if in_actual and in_bb is None:
                fn+=1
            elif in_bb is not None and not in_actual:
                fp+=1
            else:
                continue
        fn_rate = fn/(fn+tp)
        return fn_rate
    
class MPH:
    def __init__(self, keys):
        self.keys = keys
        self.seed = 1
        self.sec = wyhash.make_secret(self.seed)
        self.hashkeys = [wyhash.hash(k.encode("UTF-8"), self.seed, self.sec) for k in keys]
        # [mmh3.hash(x, signed =False) for x in keys]
        # [wyhash.hash(k.encode("UTF-8"), self.seed, self.sec) for k in keys]
        # [xxhash.xxh64_intdigest(k) for k in keys]                
        self.fingerprints = None
        self.mphf = None
    def build_mph(self):
        self.mphf = bbhash.PyMPHF(self.hashkeys, len(self.hashkeys), 1, 1.0)
       
    def query_bb(self, query_keys):
        bb_quer = []
        for k in query_keys:
            hv = wyhash.hash(k.encode("UTF-8"), self.seed, self.sec)
            bb_quer.append(self.mphf.lookup(hv))
        return bb_quer
    def build_fingerprints(self, nbits):
        fingerprints = np.empty(len(self.keys)).astype(int)
        for i in range(len(self.keys)):
            mp = self.mphf.lookup(self.hashkeys[i])
            fprint = extract_bbits(mp, nbits, nbits)
            fingerprints[mp] = fprint
        self.fingerprints = fingerprints
        self.hashkeys = None
    def query_bb_fingerprint(self, query_keys, nbits):
        bb_quer = []
        totg = 0
        for k in query_keys:
            hv = wyhash.hash(k.encode("UTF-8"), self.seed, self.sec)
            # mmh3.hash(k, signed =False)
            # hv = xxhash.xxh64(k.encode("UTF-8"))
            # hv = xxhash.xxh64_intdigest(k)
            mphval = self.mphf.lookup(hv)
            if mphval:
                ev = extract_bbits(hv, nbits, nbits)
                if(ev>len(self.fingerprints)):
                    mphval = None
                    totg+=1
                    print("MORE")
                elif (ev != self.fingerprints[mphval]):
                    mphval = None
            bb_quer.append(mphval)
        return bb_quer
    
    
    def compute_fp_rate(self, guess, query_keys, tn):
        fp = 0
        fn = 0
        for i in range(len(query_keys)):
            in_actual = query_keys[i] in self.keys
            in_bb = guess[i]
            if in_actual and in_bb is None:
                fn+=1
            elif in_bb is not None and not in_actual:
                fp+=1
            else:
                continue
        fp_rate = fp/(fp+tn)
        return fp_rate
    def compute_fn_rate(self, guess, query_keys, tp):
        fp = 0
        fn = 0
        count_none = 0
        for i in range(len(query_keys)):
            in_actual = query_keys[i] in self.keys
            in_bb = guess[i]
            if(in_bb is None):
                count_none+=1
            if in_actual and in_bb is None:
                fn+=1
            elif in_bb is not None and not in_actual:
                fp+=1
            else:
                continue
       
        fn_rate = fn/(fn+tp)
        return fn_rate



def extract_bbits(number, k, p):
    shifted_number = int(bin((number >> (p-1)) & ((1 << k) - 1))[2:],2)
    return shifted_number

def main():
    print(sys.argv[2])
    kmerTruef = open(sys.argv[1], "r")
    # line = kmerTruef.readline()
    kmersT = []
    while True:
        line = kmerTruef.readline()
        line = line.strip("\n")
        if line=='':
            break
        kmersT.append(line)
    kmerTruef.close()
    kmerFpf = open(sys.argv[2], "r")
    kmersFp = []
    while True:
        line = kmerFpf.readline()
        line = line.strip("\n")
        if line=='':
            break
        kmersFp.append(line)
    kmerFpf.close()
    error = float(sys.argv[3])
    nbits = int(sys.argv[4])
    true_num_pos = int(sys.argv[5])
    true_num_neg = int(sys.argv[6])

    ### TASK 1 ###
    bf = Bloom(kmersT, error)
    bf.build_bloom_filter()

    #Timing bf query...
    timeS = time.time_ns()

    guesses = bf.query_bf(kmersFp)
    timeelapsedquery = time.time_ns()-timeS
    print("Time elapsed for kmer query: "+ str(timeelapsedquery))
    print("Size of bloom filter: ", asizeof.asizeof(bf.bloom))

    #Building fp rate
    fp_rate = bf.compute_fp_rate(guesses, kmersFp, true_num_neg)
    print("False Positive Rate for error TASK1:"+str(error)+" is " + str(fp_rate) +" On query kfp of size: "+str(len(kmersFp)))
    fn_rate = bf.compute_fn_rate(guesses, kmersFp, true_num_pos)
    print("False Negative Rate for error TASK1:"+str(error)+" is " + str(fn_rate) +" On query kfp of size: "+str(len(kmersFp)))


    
    
    
    ## TASK 2 ####
    mph = MPH(kmersT)
    mph.build_mph()


    #Timing bb query...
    timeS = time.time_ns()
    guesses = mph.query_bb(kmersFp)
    timeelapsedquery = time.time_ns()-timeS
    print("Time elapsed for mph kmer query: "+ str(timeelapsedquery))
    print("Size of MPHF: ", asizeof.asizeof(mph.mphf))
    fp_rate = mph.compute_fp_rate(guesses, kmersFp, true_num_neg)
    print("False Positive Rate for error TASK2:"+str(error)+" is " + str(fp_rate) +" On query kfp of size: "+str(len(kmersFp)))
    fn_rate = mph.compute_fp_rate(guesses, kmersFp, true_num_pos)
    print("False Negative Rate for error TASK2:"+str(error)+" is " + str(fn_rate) +" On query kfp of size: "+str(len(kmersFp)))



    #### TASK 3 ####
    mph.build_fingerprints(nbits)

    #Timing bb+fp query...
    timeS = time.time_ns()
    guesses = mph.query_bb_fingerprint(kmersFp, nbits)
    timeelapsedquery = time.time_ns()-timeS
    print("Time elapsed for mph+finger kmer query: "+ str(timeelapsedquery))
    mph.mphf.save("size_"+str(len(kmersT))+".out")
    print("Size of MPHF+FP: ", asizeof.asizeof(mph.mphf)+ asizeof.asizeof(mph.fingerprints))
    print("Size of MPHF+FPPACKED: ", asizeof.asizeof(mph.mphf)+ asizeof.asizeof(np.packbits(mph.fingerprints)))
    fp_rate = mph.compute_fp_rate(guesses, kmersFp, true_num_neg)
    print("False Positive Rate for error TASK3:"+str(error)+" is " + str(fp_rate) +" On query kfp of size: "+str(len(kmersFp)))
    fn_rate = mph.compute_fn_rate(guesses, kmersFp, true_num_pos)
    print("False Negative Rate for error TASK3:"+str(error)+" is " + str(fn_rate) +" On query kfp of size: "+str(len(kmersFp)))




if __name__ =='__main__':
    main()