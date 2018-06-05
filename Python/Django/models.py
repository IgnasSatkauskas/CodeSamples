from __future__ import unicode_literals
from django.utils.encoding import python_2_unicode_compatible

from django.db import models
from django.utils import timezone

from accounts.models import User

#for date fields
from datetime import date
import datetime

# Create your models here.

@python_2_unicode_compatible  # only if you need to support Python 2 (for __str__ method)
class Autocall(models.Model):

# max call attempts, no need to makemigrations/migrate upon changing
    MAX_ATTEMPTS = 3

# one-to-one relationship with user
# places and restaurants
    user = models.OneToOneField(
      User,
      on_delete = models.CASCADE,
      primary_key = True,
    )

 
    testline_phone = models.CharField(
        max_length=12,
        help_text = 'Enter the phone number you have to call daily'
    )

    @property
    def testline_phone_friendly(self):
        area = self.testline_phone[2:5]
        prefix = self.testline_phone[5:8]
        number = self.testline_phone[8:12]
        return '(' + area + ') ' + prefix + '-' + number 
    # can use it as: if instance.testline_number_friendly 


    testline_id = models.CharField(
        max_length=20,
        verbose_name = 'Testline ID',
        help_text = 'Enter your ID number'
    )
    last_name = models.CharField(
        max_length=3,
        verbose_name = 'Lastname',
        help_text = 'Enter the first 3 letters of your lastname you have to confirm during the call'
    )


# twilio call sid
    twilio_call_sid = models.CharField(
                max_length = 255,
                blank = True
    )

# twilio recording sid
    twilio_rec_sid = models.CharField(
                max_length = 255,
                blank = True
    )

# call time (better way)
    AM3 = datetime.time(hour=3) 
    AM4 = datetime.time(hour=4)
    AM5 = datetime.time(hour=5)
    AM6 = datetime.time(hour=6)
    AM7 = datetime.time(hour=7)
    AM8 = datetime.time(hour=8) 
    AM9 = datetime.time(hour=9) 
    AM10 = datetime.time(hour=10) 
    AM11 = datetime.time(hour=11) 
    PM12 = datetime.time(hour=12) 
    PM1 = datetime.time(hour=13) 
    PM2 = datetime.time(hour=14) 
    PM3 = datetime.time(hour=15) 
    PM4 = datetime.time(hour=16) 
    PM5 = datetime.time(hour=17) 
    PM6 = datetime.time(hour=18) 
    PM7 = datetime.time(hour=19) 
    PM8 = datetime.time(hour=20) 
    PM9 = datetime.time(hour=21) 
    TIME_CHOICES = (
      (AM3,'between 3AM and 4AM'), 
      (AM4,'between 4AM and 5AM'), 
      (AM5,'between 5AM and 6AM'), 
      (AM6,'between 6AM and 7AM'), 
      (AM7,'between 7AM and 8AM'), 
      (AM8,'between 8AM and 9AM'), 
      (AM9,'between 9AM and 10AM'), 
      (AM10,'between 10AM and 11AM'), 
      (AM11,'between 11AM and 12PM'), 
      (PM12,'between 12PM and 1PM'), 
      (PM1,'between 1PM and 2PM'), 
      (PM2,'between 2PM and 3PM'), 
      (PM3,'between 3PM and 4PM'), 
      (PM4,'between 4PM and 5PM'), 
      (PM5,'between 5PM and 6PM'), 
      (PM6,'between 6PM and 7PM'), 
      (PM7,'between 7PM and 8PM'), 
      (PM8,'between 8PM and 9PM'), 
      (PM9,'between 9PM and 10PM'), 
    )

    time = models.TimeField(
      choices = TIME_CHOICES,
      default = AM5,
      verbose_name = 'call time', 
      help_text = 'Select a time window for your automatic call'
    )

# for agg calls: different time for different days of the week

# MON
# make call? 
    mon_call = models.BooleanField(
      default = True 
    )
# calling time
    mon_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'mon call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    mon_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'mon call time default', 
    )

# TUE
# make call? 
    tue_call = models.BooleanField(
      default = True 
    )
# calling time
    tue_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'tue call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    tue_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'tue call time default', 
    )

# WED 
# make call? 
    wed_call = models.BooleanField(
      default = True 
    )
# calling time
    wed_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'wed call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    wed_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'wed call time default', 
    )

# THU 
# make call? 
    thu_call = models.BooleanField(
      default = True 
    )
# calling time
    thu_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'thu call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    thu_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'thu call time default', 
    )

# FRI 
# make call? 
    fri_call = models.BooleanField(
      default = True 
    )
# calling time
    fri_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'fri call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    fri_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'fri call time default', 
    )


# SAT 
# make call? 
    sat_call = models.BooleanField(
      default = True 
    )
# calling time
    sat_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'sat call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    sat_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'sat call time default', 
    )


# SUN 
# make call? 
    sun_call = models.BooleanField(
      default = True 
    )
# calling time
    sun_time = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'sun call time', 
    )
# default call time: used by prepwork, to return call time to original value
# in case of "your schedule is not currently available" 
    sun_time_default = models.TimeField(
      choices = TIME_CHOICES,
      default = PM5,
      verbose_name = 'sun call time default', 
    )



# call time timezone
    EST = 'US/Eastern'
    CST = 'US/Central'
    MST = 'US/Mountain'
    PST = 'US/Pacific' 

    TIMEZONE_CHOICES = (
        (EST,'Eastern Time (ET)'), 
        (CST,'Central Time (CT)'), 
        (MST,'Mountain Time (MT)'), 
        (PST,'Pacific Time (PT)'), 
    ) 

    call_timezone = models.CharField(
        max_length = 30,
        choices = TIMEZONE_CHOICES,
        default = MST,
        verbose_name = 'Time zone', 
        help_text = 'Select your time zone'
    )


# call time interval (for testing purposes really...)
    t_start = models.TimeField(
      default = datetime.time(hour=5),
      verbose_name = 't_start', 
    )

    t_end = models.TimeField(
      default = datetime.time(hour=12),
      verbose_name = 't_end', 
    )

# is_time_ok  
# check if it is ok to call at the time of the check (timenow is in users selected time gap)
    @property
    def is_time_ok(self):
        timezone.activate(self.call_timezone)
        datetime_now = timezone.localtime(timezone.now()) # this is now in users local time
        timezone.deactivate()
        time_now = datetime_now.time() # just get the datetime.time object, it is naive now, can compare to autocall.time
        # self.time - start in local time (datetime.time object),  
        # datetime.time(hour = self.time.hour + 1) - end in local time (an hour later)   
        # doing datetime.time(hour=x) can cause ValueError if x is not in 0..23, but
        # in our case is ok because self.time.hour is max 21 and we're adding 1
        if time_now > self.time and time_now < datetime.time(hour = self.time.hour + 1 ):
            return True
        else:
            return False
        
# is_time_ok_custom  
# check if it is ok to call at the time of the check (timenow is in users selected time gap)
# for custom calls, check time depending of which day of the week it is 
    @property
    def is_time_ok_custom(self):
        timezone.activate(self.call_timezone)
        datetime_now = timezone.localtime(timezone.now()) # this is now in users local time
        timezone.deactivate()
        time_now = datetime_now.time() # just get the datetime.time object, it is naive now, can compare to autocall.time
        # self.time - start in local time (datetime.time object),  
        # datetime.time(hour = self.time.hour + 1) - end in local time (an hour later)   
        # doing datetime.time(hour=x) can cause ValueError if x is not in 0..23, but
        # in our case is ok because self.time.hour is max 21 and we're adding 1
        weekday_now = datetime_now.weekday()
        # we also need the day of the week for agg calls
        # returns integer: 0 for Mon, 1 for Tue, ..., 6 for Sun
        if weekday_now == 0:
        # if it is Mon
            if self.mon_call:
            # boolean, make call on this day?  
                call_time = self.mon_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 1:
        # if it is Tue 
            if self.tue_call:
            # boolean, make call on this day?  
                call_time = self.tue_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 2:
        # if it is Wed 
            if self.wed_call:
            # boolean, make call on this day?  
                call_time = self.wed_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 3:
        # if it is Thu 
            if self.thu_call:
            # boolean, make call on this day?  
                call_time = self.thu_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 4:
        # if it is Fri 
            if self.fri_call:
            # boolean, make call on this day?  
                call_time = self.fri_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 5:
        # if it is Sat 
            if self.sat_call:
            # boolean, make call on this day?  
                call_time = self.sat_time
            else:
            # if no call on this day, don't even check the time
                return False
        elif weekday_now == 6:
        # if it is Sun 
            if self.sun_call:
            # boolean, make call on this day?  
                call_time = self.sun_time
            else:
            # if no call on this day, don't even check the time
                return False
            
        # if above didn't return False, the we have call_time for today  
        #if time_now > self.time and time_now < datetime.time(hour = self.time.hour + 1 ):
        if time_now > call_time and time_now < datetime.time(hour = call_time.hour + 1 ):
            return True
        else:
            return False

# is_time_ok_test (for testing purposes) 
# check if it is ok to call at the time of the check (timenow is in testline's 'working hours' interval)
# t_start and t_end should be verified, but default values (5 and 12) should work in most cases  
    @property
    def is_time_ok_test(self):
        timezone.activate(self.call_timezone)
        datetime_now = timezone.localtime(timezone.now()) # this is now in users local time
        timezone.deactivate()
        time_now = datetime_now.time() # just get the datetime.time object, it is naive now, can compare to autocall.time
        # self.t_start - start of testline 'working hours' in local time (datetime.time object),  
        # self.t_end - end of testline 'working hours' in local time (datetime.time object),  
        if time_now > self.t_start and time_now < self.t_end:
            return True
        else:
            return False
        
# is_time_too_late  
# check now (time of check) is between midnight and  user's selected call interval end (minus 5 minutes, process/calls comes in every 2 min at least)
# used when user activates the call, to display one of the  messages: 
# "call active, will be made TODAY as scheduled" or "call active, next call will be made TOMOROW as scheduled"
    @property
    def is_time_too_late(self):
        timezone.activate(self.call_timezone)
        datetime_now = timezone.localtime(timezone.now()) # this is now in users local time
        timezone.deactivate()
        time_now = datetime_now.time() # just get the datetime.time object (users local time), it is naive now, can compare to autocall.time
 
        if time_now > datetime.time(hour = 0) and time_now < datetime.time(self.time.hour,55 ):
            return False # not too late 
        else:
            return True # too late today 
        

# call result 
    VERIFYING = 'verifying'
    ERR_PARAMS = 'err:params'
    NORESULT_INACTIVE = 'noresult:inactive'
    ERR_ATTEMPTS = 'err:attempts'
    NORESULT = 'noresult'
    CALLING = 'calling'
    READY2PROCESS = 'ready2process'
    PROCESSING = 'processing'
    DONOT = 'donot'
    REQ = 'req'
    EXP = 'exp'
    VER = 'ver'
    TRANSCRIPTION = 'transcription'

    RESULT_CHOICES = (
      (VERIFYING, 'Verifying call parameters...'),
      (ERR_PARAMS , 'Error: call parameters'),
      (NORESULT_INACTIVE , 'No result, call is inactive'),
      (ERR_ATTEMPTS , 'Error: no result after ' + str(MAX_ATTEMPTS) + ' call attempts'),
      (NORESULT , 'No result yet'),
      (CALLING,'Calling...'),
      (READY2PROCESS,'Ready to process'),
      (PROCESSING,'Processing...'),
      (DONOT , 'Do not test today'),
      (REQ , 'You are required to test today'),
      (EXP , 'Your testline id # has expired'),
      (VER , 'Error: Your testline id # could not be verified'),
      (TRANSCRIPTION , 'transcription'),
    )



    result = models.CharField(
      max_length = 255, 
      choices = RESULT_CHOICES,
      default = VERIFYING
    )

    # property to check if call is being processed at the moment: calling, ready2process, processing
    # used in dashboard for displaying result, in autocall/deactivate_call_view 
    @property
    def is_processing(self):
        if self.result == self.CALLING or self.result == self.READY2PROCESS or self.result == self.PROCESSING:
            return True 
        else:
            return False 
 
    result_timestamp = models.DateTimeField(
        default = timezone.now,
    )

    # property @is_verified 
    @property
    def is_verified(self):
        if self.result == self.VERIFYING or self.result == self.ERR_PARAMS:
            return False 
        else:
            return True
    # can use it as: if instance.is_verified 

# call attempts made     
    attempts = models.PositiveSmallIntegerField(
      default = 0
    )

# is call active 
    active = models.BooleanField(
      default = False 
    )


# already emailed/texted the result?
    emailed = models.BooleanField(
      default = False
    )

    texted = models.BooleanField(
      default = False
    )



    def __str__(self):
        return u'%s %s %s %s' %(self.user, self.testline_phone, self.testline_id, self.last_name)

    def copy_to_used(self):
        
        used_autocall = UsedAutocall.objects.create(testline_phone=self.testline_phone,
                                                    testline_id=self.testline_id,
                                                    last_name=self.last_name,
        ) 
        used_autocall.save()

#----------------------------------------------------------------------------------------------
# extending Autocall model to include Custom call and Scheduling Aggregate Industries
#---------------------------------------------------------------------------------------------

# ivr system type

    SENTRY = 'sentry'
    AGGREGATE = 'aggregate'
    CUSTOM = 'custom'

    IVR_TYPE_CHOICES = (
      (SENTRY, 'Sentry Notification System'),
      (AGGREGATE , 'Aggregate Industries Scheduling System'),
      (CUSTOM, 'Custom Call'),
    )

    ivr_type = models.CharField(
      max_length = 255, 
      choices = IVR_TYPE_CHOICES,
      default = SENTRY 
    )

# additional fields to support aggregate and custom calls

# number of steps (don't nkow if I'll use it, up to 3 at the moment)

    number_of_steps = models.PositiveSmallIntegerField(
      default = 0
    )

# steps, consisting of promt, enter_keys, wait_time
# i.e. directions for the call
    
    promt1 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter the promt you hear during the call',
    )


    promt2 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter the promt you hear during the call',
    )

    promt3 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter the promt you hear during the call',
    )

    enter_keys1 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter dial pad keys (1234567890#*)',
    )

    enter_keys2 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter dial pad keys (1234567890#*)',
    )

    enter_keys3 = models.CharField(
        max_length = 255,
        blank = True, # for CharField this means '' 
        help_text = 'Enter dial pad keys (1234567890#*)',
    )


    wait_time1 = models.PositiveSmallIntegerField(
      default = 0
    )

    wait_time2 = models.PositiveSmallIntegerField(
      default = 0
    )

    wait_time3 = models.PositiveSmallIntegerField(
      default = 0
    )

    TRANSCRIBE_SECONDS_MAX = 30


    transcribe_seconds = models.PositiveSmallIntegerField(
      default = 30
    )

    transcription = models.CharField(
      max_length = 500,
      blank = True,
    )


    # property @call_steps 
    @property
    def call_steps(self):
        if self.enter_keys3 != '' and self.enter_keys2 != '' and self.enter_keys1 != '':
            return ( (self.promt1, self.enter_keys1, self.wait_time1 ), (self.promt2, self.enter_keys2, self.wait_time2 ), (self.promt3, self.enter_keys3, self.wait_time3 ) ) 
        if self.enter_keys3 == '' and self.enter_keys2 != '' and self.enter_keys1 != '':
            return ( (self.promt1, self.enter_keys1, self.wait_time1 ), (self.promt2, self.enter_keys2, self.wait_time2 ), ) 
        if self.enter_keys3 == '' and self.enter_keys2 == '' and self.enter_keys1 != '':
            return ( (self.promt1, self.enter_keys1, self.wait_time1), ) 
        else:
            return None 


# transcription fixing--------------------------------------------------------------------------------
# fix twilio transcription
# common mistakes
 
    @property
    def fix_twilio_transcription(self):
        
        #transcription = self.transcription

        replace_words = {'bachelor':'batching',
                         'Bacha':'batching',
                         'bacha a':'batching',
                         'bacha':'batching',
                         ' yr ': ' your ',
                         ' here ': ' your ',
                         # this error always comes up (can't generalize it yet, since 'to' comes up elsewhere)
                         '17 to 0':'1720',
                         ' o ':'0',
                         ' O ':'0',
                         ' Oh ':'0',
                         ' oh ':'0',
                        # ' to ': '2',
                        # 'to' might show up instead of two(2) in 17 to 0, but it might not
                        # and if not, then next 'to' is ' to sign up press'
        }

 
        #for key, value in replace_words.items(): # for python 3
        for key, value in replace_words.iteritems():
            self.transcription =  self.transcription.replace(key,value,1)
        
            
        self.save()
        # don't think I need to return anything, but for testing purposes
        #return self.transcription

# check if transcription is valid 
# for example, this is not:
# "Sorry that password is not valid Please enter your password followed by the pound sign 
# Please enter your password followd by the pound sign Please enter your password followed by the pound sign. "
# use this in process/transcription_view

    @property
    def is_twilio_transcription_valid(self):
        
        wrong_words = ['password',]

        for word in wrong_words:
            # see if transcription has any of wrong_words
            # str.find() will return begining index of the word if looks for
            # if no word found, returns -1 
            if self.transcription.find(word) > -1:
                return False

        return True

#--------------------------------------------------------------------------------------------------
# Autocalls that were used before, for keeping records of phone #, id #, and lastname (first 3 letters)
#--------------------------------------------------------------------------------------------------

@python_2_unicode_compatible  # only if you need to support Python 2 (for __str__ method)
class UsedAutocall(models.Model):


    testline_phone = models.CharField(
        max_length=12,
        help_text = 'The phone number of the testline you have to call daily.'
    )


    testline_id = models.CharField(
        max_length=20,
        verbose_name = 'Testline ID',
        help_text = 'The ID number you enter during the call. '
    )

    last_name = models.CharField(
        max_length=3,
        verbose_name = 'Lastname',
        help_text = 'The first 3 letters of your lastname you have to confirm during the call.'
    )

    # this will be stamped upon creation (by admin after verification process)
    created_timestamp = models.DateTimeField(
        default = timezone.now,
    )

    deleted_timestamp = models.DateTimeField(
        blank = True,
        null = True,
    )

    def __str__(self):
        return u'%s %s %s' %( self.testline_phone, self.testline_id, self.last_name)



