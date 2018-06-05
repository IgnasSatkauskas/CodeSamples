#!/usr/bin/python
#
#------------------------------------------------------------------------------------------------
# old school prototype code 
# "skype-dialer"
#------------------------------------------------------------------------------------------------

import sys
sys.path.append("/home/ignas/Notify4/")
import autocallbot as acb

import pyaudio
import wave
import datetime
import time
import subprocess
import numpy as np
from scipy.io import wavfile
import MySQLdb
import logging
import audioop


logging.basicConfig(filename='/home/ignas/Notify4/test_run3.log',level=logging.DEBUG)
#db name to use and total number of users +1
TABLE_NAME='MyUsers'
#TABLE_NAME = 'TestUsers3'
#TABLE_NAME = 'TestUsers8pm'
USER_TOTAL = acb.fetch_db_total(TABLE_NAME) + 1

#parameters for recording
CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 2
RATE = 8000
RECORD_SECONDS = 14
#WAVE_OUTPUT_FILENAME = "/home/ignas/Python/Sound/Tracks/8k/test.wav"


#FOR PROCESSING (RECOGNITION)
#WORD1 = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/do_not.wav'
#WORD2 = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/you_are_req.wav'
#FULL = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/do_not_test_full_12s.wav'
#FULL = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/you_are_required_full_12s.wav'
#FULL = '/home/ignas/Notify2/Tracks/test.wav'


#detect sound (needs audioop)
def detect_sound():
   np_data = np.zeros(1)
   cutoff = 800
   chunk = 1024 
   rate_listen = 8000
   seconds_listen=20

   p = pyaudio.PyAudio()

   stream = p.open(format=pyaudio.paInt16,
                   channels=2,
                   rate=rate_listen,
                   input=True,
                   frames_per_buffer=chunk)

#needs timeout for listening?

#   while True:
   for i in range(0, int(rate_listen / chunk * seconds_listen)):
      data = stream.read(chunk)
      rms = audioop.rms(data, 2)  #width=2 for format=paInt16
      if rms > cutoff:
         sound_detected = 1
         break
      sound_detected = 0 # no sound in seconds_listen
#      np_data = np.append(np_data,rms)


   stream.stop_stream()
   stream.close()

   p.terminate()
   return(sound_detected)
#----------------------------------------------------------------------------------


#-----------------------------------------------------------------
def record(IdNumber,PhoneNumber):
   
   wav_output_filename = "/home/ignas/Notify4/Tracks/test" + str(IdNumber) + ".wav" 

   p = pyaudio.PyAudio()
   stream = p.open(format=FORMAT,
                   channels=CHANNELS,
                   rate=RATE,
                   input=True,
                   frames_per_buffer=CHUNK)
   		
   
   frames = []
   
   for i in range(0, int(RATE / CHUNK * RECORD_SECONDS)):
       data = stream.read(CHUNK)
       frames.append(data)
   
   
   stream.stop_stream()
   stream.close()
   
   wf = wave.open(wav_output_filename, 'wb')
   wf.setnchannels(CHANNELS)
   wf.setsampwidth(p.get_sample_size(FORMAT))
   wf.setframerate(RATE)
   wf.writeframes(b''.join(frames))
   wf.close()
  
   p.terminate() 
   time.sleep(1)
   
#----------------------------------------------------------------------------------

def normalize(track):
        factor = max(track)-min(track)
        track = track / factor
        return(track)
#----------------------------------------------------------------------------------


def process(IdNumber):

   WORD1 = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/do_not.wav'
   WORD2 = '/home/ignas/Python/Sound/Tracks/8kgreen2blue/you_are_req.wav'

   FULL = "/home/ignas/Notify4/Tracks/test" + str(IdNumber) + ".wav" 
   
   fs1, donot = wavfile.read(WORD1)
   donot = donot[:,0]
   
   fs2, req = wavfile.read(WORD2)
   req = req[:,0]

   try:
      fs3, full = wavfile.read(FULL)
      full = full[:,0]
      full = full / (2.**15)
      full = normalize(full) + 1.
   except:
      print "error opening file"
      logging.debug("error opening file FULL") 
      return(2) 

   donot= donot[0:1000]
   req= req[0:1000]

   donot = donot / (2.**15)
   req = req / (2.**15)

   donot = normalize(donot) + 1.
   req = normalize(req) + 1.



   norm2_donot = np.zeros([1,1])
   k=2 #shift every k samples
   #for index in range(len(full)-2*len(donot)-1):
   hits_donot = 0
   cutoff = 6
   for index in range(0,len(full)-len(donot)-k,k):
        piece = full[index:index+len(donot)]
        norm2 = np.linalg.norm(piece-donot,2)
        #norm2_donot = np.append(norm2_donot,norm2) #for printing
        if norm2 < cutoff:
                hits_donot = hits_donot+1;

   hits_req = 0
   cutoff = 6
   norm2_req = np.zeros([1,1])
   for index in range(0,len(full)-len(req)-k,k):
        piece = full[index:index+len(req)]
        norm2 = np.linalg.norm(piece-req,2)
        #norm2_req = np.append(norm2_req,norm2) #for printing
        if norm2 < cutoff:
                hits_req = hits_req+1;

    
   print "below 6 donot:", hits_donot
   print "below 6 req:", hits_req
   logging.debug("hits_donot= " + str(hits_donot) + "  hits_req= " + str(hits_req))
   if hits_req > 80:
      return(1)
   elif hits_donot > 80: 
      return(0)
   else:
      return(2)

#---------------------------------------------------------------------------------
#new processing (no cable and gap measuring bw instances)

def process2(IdNumber):

   WORD1 = '/home/ignas/Notify4/Tracks/NoCable/donot.wav'
   WORD2 = '/home/ignas/Notify4/Tracks/NoCable/req.wav'

   FULL = "/home/ignas/Notify4/Tracks/test" + str(IdNumber) + ".wav"

   fs1, donot = wavfile.read(WORD1)
   donot = donot[:,0]

   fs2, req = wavfile.read(WORD2)
   req = req[:,0]

   try:
      fs3, full = wavfile.read(FULL)
      full = full[:,0]
      full = full / (2.**15)
      full = normalize(full) + 1.
   except:
      print "error opening file"
      return(2)

   donot= donot[0:1000]
   req= req[0:1000]

   donot = donot / (2.**15)
   req = req / (2.**15)

   donot = normalize(donot) + 1.
   req = normalize(req) + 1.


#----------against donot------------------------------------------
   norm2_donot = np.zeros([1,1])
   # indeces of hits_donot (for measuring the distance bw instances)
   hits_donot_ind = np.zeros([1,1])
   k=2 #shift every k samples
   #for index in range(len(full)-2*len(donot)-1):
   hits_donot = 0
   cutoff = 6
   for index in range(0,len(full)-len(donot)-k,k):
        piece = full[index:index+len(donot)]
        norm2 = np.linalg.norm(piece-donot,2)
        norm2_donot = np.append(norm2_donot,norm2) #for printing
        if norm2 < cutoff:
                hits_donot = hits_donot+1;
                hits_donot_ind = np.append(hits_donot_ind,index)

#----------against required------------------------------------------
   norm2_req = np.zeros([1,1])
   # indeces of hits_donot (for measuring the distance bw instances)
   hits_req_ind = np.zeros([1,1])
   hits_req = 0
   cutoff = 6
   for index in range(0,len(full)-len(req)-k,k):
        piece = full[index:index+len(req)]
        norm2 = np.linalg.norm(piece-req,2)
        norm2_req = np.append(norm2_req,norm2) #for printing
        if norm2 < cutoff:
                hits_req = hits_req+1;
                hits_req_ind = np.append(hits_req_ind,index)


   print "below 6 donot:", hits_donot
   print "below 6 req:", hits_req

   logging.debug("hits_donot= " + str(hits_donot) + "  hits_req= " + str(hits_req))
# print indeces of donot and req instances   
   #print "hits_donot_ind", hits_donot_ind
   #print "hits_req_ind", hits_req_ind


   # calculate distances bw hits in DONOT case ---------------------------------
   if (hits_donot > 40) & (hits_req == 0):
      #take hits_donot_ind and find two max distances bw them
      x = hits_donot_ind
      y = x[1:len(x)-1]
      z = x[2:len(x)]
      d = z-y
      d = np.sort(d)
      max1 = d[len(d)-1]
      max2 = d[len(d)-2]
      max3 = d[len(d)-3] # just to see that 3rd max is much smaller

#      print "d sorted ", d
      print "max1, max2, max3 ", max1, max2, max3
      logging.debug("max1, max2, max3 " + str(max1) + " " + str(max2) + " " +  str(max3) )
      if (max1 > 15000) & (max2 > 15000):
         return(0)
      else:
         logging.debug("ERROR PROCESSING: maxes")
         return(2)

   # calculate distances bw hits in REQUIRED case ---------------------------------

   elif (hits_req > 80) & (hits_donot == 0):
      #take hits_req_ind and find two max distances bw them
      x = hits_req_ind
      y = x[1:len(x)-1]
      z = x[2:len(x)]
      d = z-y
      d = np.sort(d)
      max1 = d[len(d)-1]
      max2 = d[len(d)-2]
      max3 = d[len(d)-3] # just to see that 3rd max is much smaller

#      print "d sorted ", d
      print "max1, max2, max3 ", max1, max2, max3
      logging.debug("max1, max2, max3 " + str(max1) + " " + str(max2) + " " +  str(max3) )
      if (max1 > 15000) & (max2 > 15000):
         return(1)
      else:
         logging.debug("ERROR PROCESSING: maxes")
         return(2)
   else:
      logging.debug("ERROR PROCESSING: hits")
      return(2)


#----------------------------------------------------------------------------------
#connect to database and retrieve IdNuber given Id(number in the database)
def fetch_IdNumber(Id):
   # open database connection
   #db = mysqldb.connect("localhost","user","pass","database")
   db = MySQLdb.connect("localhost","root","Apiat&20piat&","autocallbot")
   
   #prepare a cursor object 
   cursor = db.cursor()
   
   #prepare SQL query 
   #sql = "SELECT email FROM MyUsers WHERE id=2;"
   sql = "SELECT IdNumber FROM " + TABLE_NAME + " WHERE " + "Id=" + str(Id) + ";"
   
   try:
      #execute SQL command
      cursor.execute(sql)
      #fetch all the rows in a list of lists
      output = cursor.fetchone()
      #result = cursor.fetchall()
   
      #print "result type", type(result)
      return(output[0])
      #print type(email)
      print "sql output", output 
      
   except:
      print "Error: unable to fetch data"
      logging.debug("Unable to fetch data")
   
   #disconnect from database
   db.close()   
   
#----------------------------------------------------------------------------------
#connect to database and retrieve PhoneNumber given Id(number in the database)
def fetch_PhoneNumber(Id):
   # open database connection
   #db = mysqldb.connect("localhost","user","pass","database")
   db = MySQLdb.connect("localhost","root","Apiat&20piat&","autocallbot")
   
   #prepare a cursor object 
   cursor = db.cursor()
   
   #prepare SQL query 
   #sql = "SELECT email FROM MyUsers WHERE id=2;"
   sql = "SELECT PhoneNumber FROM " + TABLE_NAME + " WHERE " + "Id=" + str(Id) + ";"
   
   try:
      #execute SQL command
      cursor.execute(sql)
      #fetch all the rows in a list of lists
      output = cursor.fetchone()
      #result = cursor.fetchall()
   
      #print "result type", type(result)
      return(output[0])
      #print type(email)
      print "sql output", output 
      
   except:
      print "Error: unable to fetch data"
      logging.debug("Unable to fetch data")
   
   #disconnect from database
   db.close()   
   
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
#connect to DB and INSERT result 
def insert_Result(IdNumber, Result):
   # open database connection
   #db = mysqldb.connect("localhost","user","pass","database")
   db = MySQLdb.connect("localhost","root","Apiat&20piat&","autocallbot")
   
   #prepare a cursor object 
   cursor = db.cursor()
   
   #prepare SQL query 
   #sql = "SELECT id8dig FROM TestUsers WHERE " + "id=" + str(number) + ";"
   sql = "UPDATE " + TABLE_NAME + " SET Result=" + str(Result) + " WHERE IdNumber=" + str(IdNumber) + ";" 
   #print "sql insert statement", sql 
   
   try:
      #execute SQL command
      cursor.execute(sql)
      #fetch all the rows in a list of lists
      #outcome = cursor.fetchone()
      #result = cursor.fetchall()
   
      #print "sql outcome", outcome 
      print "Result inserted" 
   except:
      print "Error: unable to insert data"
      logging.debug("Unable to insert data")
  
   #commit changes
   db.commit() 
   #disconnect from database
   db.close() 


#----------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
#connect to DB and fetch entire Result column 
def fetch_all_Result():
   # open database connection
   #db = mysqldb.connect("localhost","user","pass","database")
   db = MySQLdb.connect("localhost","root","Apiat&20piat&","autocallbot")

   #prepare a cursor object 
   cursor = db.cursor()

   #prepare SQL query 
   sql = "SELECT Result FROM " + TABLE_NAME + ";"
   #sql = "UPDATE TestUsers SET Result=" + str(Result) + " WHERE IdNumber=" + str(IdNumber) + ";"
   #print "sql insert statement", sql 

   try:
      #execute SQL command
      cursor.execute(sql)
      #fetch all the rows in a list of lists
      #outcome = cursor.fetchone()
      result = cursor.fetchall()

      #print "sql outcome", outcome 
      print "Results fetched=", result
      logging.debug("RESULT COLUMN=" + str(result) ) 
   except:
      print "Error: unable to fetch Results"
      logging.debug("Unable to fetch Results")

   #disconnect from database
   db.close()
#--------------------------------------------------------------------------------------------------------------

#connect to DB and INSERT results as 2 (error) 
def insert_all_Result_as2():
   # open database connection
   #db = mysqldb.connect("localhost","user","pass","database")
   db = MySQLdb.connect("localhost","root","Apiat&20piat&","autocallbot")

   #prepare a cursor object 
   cursor = db.cursor()

   #prepare SQL query 
   #sql = "SELECT id8dig FROM TestUsers WHERE " + "id=" + str(number) + ";"
   sql = "UPDATE " + TABLE_NAME + " SET Result=2;" 
   #print "sql insert statement", sql 

   try:
      #execute SQL command
      cursor.execute(sql)
      #fetch all the rows in a list of lists
      #outcome = cursor.fetchone()
      #result = cursor.fetchall()

      #print "sql outcome", outcome 
      print "All Result inserted as 2"
   except:
      print "Error: unable to insert all results as 2"
      logging.debug("Unable to insert all results as 2")

   #commit changes
   db.commit()
   #disconnect from database
   db.close()

#-----------------------------------------------------------------------------------


logging.debug("--------------------------------------------------------")
logging.debug(time.asctime()) 
logging.debug("--------------------------------------------------------")
# make few tries to test how robust it is 
for i in range(1,2):
   logging.debug("TRY--------------------------------------------------------------" + str(i))

   for Id in range(1,USER_TOTAL):
      start_time = time.time()
      print "WORKING ON Id",Id
      logging.debug("WORKING ON Id " + str(Id) + "----------------------")
      print 'check if call is allowed at this time and is needed'
      logging.debug("check if call is allowed at this time and is needed")
      # find current hour
      now = datetime.datetime.now()
      nowhour = now.hour
      # find call gap (allowed calling times)
      tstart = acb.fetch_db(TABLE_NAME,Id,'tstart')
      tend = acb.fetch_db(TABLE_NAME,Id,'tend')
      # fetch current Result 
      current_Result = acb.fetch_db(TABLE_NAME,Id,'Result')
      # if call is allowed, and is needed, proceed... 
      if (nowhour >= tstart) & (nowhour < tend) & (current_Result==2):
         print 'call allowed at this time and is needed, proceeding...'
         logging.debug('call allowed at this time and is needed, proceeding...')

         #fetch IdNumber and PhoneNumber
         IdNumber = fetch_IdNumber(Id) 
         print "IdNumber= ", IdNumber 
         PhoneNumber = fetch_PhoneNumber(Id) 
         print "PhoneNumber= ", PhoneNumber 
         logging.debug("IdNumber= " + str(IdNumber))
         logging.debug("PhoneNumber= " + str(PhoneNumber))
      
         #make 3 attempts if necessary
         for attempt in range(1,4):
            time.sleep(1)
            print "attempt=", attempt
            logging.debug("attempt= " + str(attempt))
            #make call
            logging.debug("Calling")
            #login to Skype
            subprocess.call(["/home/ignas/Notify2/skype_login.sh"])
            time.sleep(10)
            #make Skype call
            subprocess.call(["skype", "--callto", "+1" + PhoneNumber])
            #detect sound
            print "Detecting sound..."
            sound_detected = detect_sound()
            #print "sound_detected=", sound_detected
            if sound_detected:
               print "Sound detected"
            else: 
               print "No sound detected"
               #kill Skype
               subprocess.call(["killall", "-9", "skype"])
               time.sleep(1)
               continue
            time.sleep(5)
            #run xdotool sequence
            print "running xdotool"
            subprocess.call(["/home/ignas/Notify4/xdo.sh", PhoneNumber, IdNumber]) 
            print "done running xdotool"
            print "recording"
            record(IdNumber,PhoneNumber)
            print "done recording" 
            #kill Skype
            subprocess.call(["killall", "-9", "skype"])
            time.sleep(1)
            #process
            print "processing..."
            Result = process(IdNumber)
            print "done processing"
      
            print "Result=", Result 
            logging.debug("Result= " + str(Result))
            if Result == 0:
               print "DO NOT TEST"
               insert_Result(IdNumber,Result)
               logging.debug("Inserted Result =" + str(Result))
               break
            elif Result == 1:
               print "REQUIRED TO TEST"
               insert_Result(IdNumber,Result)
               logging.debug("Inserted Result =" + str(Result))
               break
            else:
               print "ERROR in PROCESSING"
               insert_Result(IdNumber,Result)
               logging.debug("Inserted Result =" + str(Result))
      #when call is not allowed at this time
      else:
         print 'call is not allowed at this time or not needed' 
         logging.debug('call is not allowed at this time or not needed')
 
      end_time = time.time()
      ellapsed_time = end_time - start_time
      print "Ellapsed time for one Id is", ellapsed_time 
      logging.debug("Ellapsed time for one Id is " + str(ellapsed_time))
   
   print "DONE with all Ids"
   logging.debug("DONE with all Ids")
   time.sleep(1)
      
            
