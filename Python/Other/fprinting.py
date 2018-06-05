from scipy.io import wavfile
import numpy as np
import itertools
import pickle
import requests

import time
import os
#print(os.path.dirname(os.path.realpath(__file__)))
BASE_DIR = os.path.dirname(os.path.realpath(__file__))

# dead simple user yes/no promt (for adding prints to a library)
def yes_or_no(question):
    while "the answer is invalid":
        reply = str(raw_input(question+' (y/n): ')).lower().strip()
        if reply[:1] == 'y':
            return True
        if reply[:1] == 'n':
            return False


# function for taking spectogram (short ffts) of the signal
#IN:
#  signal - output of wavfile.read(mono.wav), for e.g. fs1, signal = wavfile.read(Donot)
#  note, wav file should be mono, or else take one channel: signal = signal[:,0]
#  chunk_size = 400, the length of fft
#  overlap = 30, its really a step size(in time), chunk_size - overlap = true overlap
#OUT:
#  numpy 2d array: specgram[chunk_size, total_ffts]
# 
def specgram(signal, chunk_size, overlap):
  total_ffts = ( len(signal) - chunk_size ) / overlap
  specgram = np.empty([chunk_size, total_ffts])
  
  for i in range(total_ffts):
    segment = signal[overlap*i:overlap*i+chunk_size] 
    pad = np.zeros(chunk_size) #pad segment with zeros (same length)
    fft = np.fft.fft( np.append(segment,pad)  ) / float(chunk_size) # get fft of padded segment 
    fft = fft[0:chunk_size] # take only chunk_size of fft coefs
    fft = np.absolute(fft)
    specgram[:,i] = fft
   
  return(specgram)

#fingerprinting
#IN:
#  spec - spectogram of track, i.e. output of specgram() function
#  bins - list of bins where to look for frequency peaks, e.g. bins = [10, 30, 50, 70, 100] (4 bins)
#OUT:
#  dictionary of fingerprints: {fprint:time}, fprint is a tuple of location of frequency peaks
def take_fprints(spec, bins):

  total_ffts = spec.shape[1] 
#  print "in take_fprints, total_ffts= ", total_ffts
  fprints = {}
  
  #loop through all ffts
  for t in range(total_ffts): 
      
     # fingerprint - find max frequencies in bins
#     print 'in take_fprints, t = ', t
     instance = spec[:,t]
     fprint = [0, 0, 0, 0]
     for k in range(len(bins)-1):
         fprint[k] = np.argmax( instance[ bins[k]:bins[k+1] ] ) + bins[k] 
     # make fprint immutable for use as key in dict
     fprint = tuple(fprint)
     #put fprint into dictionary as key and t as its value
     #skip same fprints
     if fprint not in fprints:
       fprints[fprint] = t
  return(fprints)  





# fprint matching
#IN:
#  fprint dictionaries i.e. outputs from take_fprint()
#OUT:
#  number of matches between dict1 and dict2
def match_fprints(dict1,dict2):
  matches = 0
  for k1,k2 in itertools.combinations(dict1,2):
    if k1 in dict2 and k2 in dict2:
      if abs(dict1[k1]-dict1[k2]) == abs(dict2[k1]-dict2[k2]):
        matches = matches + 1
  return(matches)


#def process_recording(fprints_file_path, recording_file_path,k):  # testing
def process_recording(fprints_file_path, recording_file_path):
    fs, rec = wavfile.read(recording_file_path)
    # take last k-seconds of the recording
    k = 5
    rec = rec[-k*8000:] #last k seconds
    #rec = rec[-k:] # testing
    # take specgram of rec
    chunk_size = 400
    overlap = 30
    rec_spec = specgram(rec, chunk_size, overlap)
    # fprint res_spec
    bins = [10, 30, 50, 70, 100]
    rec_fprints = take_fprints(rec_spec, bins)

    #fprints_lib = load_obj(library)
    # load pickled fprints file 
    with open(fprints_file_path, 'rb') as f:
         fprints_lib =  pickle.load(f)

    #t0 = time.time()   
    # match fprints
    result = 'noresult'
    for k in fprints_lib:
        matches = match_fprints(rec_fprints, fprints_lib[k])
        #print 'in fprinting, process_recording: ', k,matches
        #optionally print out number of matches for each track in the lib
        #print k, matches # testing

        if matches > 300:
            result = k
            #result = (k,matches) # testing
    #t1 = time.time()
    #print 'time matching =', t1-t0
    return(result)







