import logging
# Get an instance of a logger
logger = logging.getLogger(__name__)

from django.shortcuts import render, redirect
from django.core.urlresolvers import reverse_lazy
from django.contrib.auth.decorators import login_required, user_passes_test
from django.utils.decorators import method_decorator
from django.views.decorators.csrf import csrf_exempt
from django.utils import timezone
from django.contrib import messages
from django.core.mail import send_mail
from django.core.mail import EmailMessage 
from django.core import mail # needed for opening SMTP connection
from django.views.decorators.http import require_POST, require_GET
from django.template.loader import render_to_string

from accounts.models import User
from autocall.models import Autocall 

import datetime
from datetime import timedelta 
import time
from django.http import HttpResponse
from django.http import JsonResponse

#from .forms import SubscriptionForm 

from django.conf import settings
import stripe
import json

import hashlib
import hmac
import base64

import requests

import os

import time

#twilio
from twilio.rest import Client
from twilio.base.exceptions import TwilioRestException
from twilio.twiml.voice_response import VoiceResponse
client = Client(settings.TWILIO_ACCOUNT_SID, settings.TWILIO_AUTH_TOKEN)


secret = settings.AUTOCALLBOT_SHARED_SECRET

# for processing recordings
from sys import path
#FPRINTING_DIR = '/Users/satkauskas/notify6/fprinting/'
path.append(settings.FPRINTING_DIR)
from fprinting import process_recording

#RECORDINGS_DIR = '/Users/satkauskas/tmp/'
#FPRINTS_DIR = '/Users/satkauskas/notify6/fprinting/obj/'


# decorator to check my signatures
def validate_request_hmac(view):
    def inner(request, *args, **kwargs):
        #hm = hmac.new(secret, request.body, hashlib.sha256)
        # add full url (with query string, if any) 
        hm = hmac.new(secret, request.build_absolute_uri()+request.body, hashlib.sha256)
        digest = base64.b64encode(hm.digest())
         
        if request.META.get("HTTP_AUTOCALLBOT_SIGNATURE", None) != digest:
            return HttpResponse(status=401)
 
        return view(request, *args, **kwargs)
 
    return inner

# decorator to check twilio POST signatures
# used for receiverec_view
def validate_twilio_post(view):
    def inner(request, *args, **kwargs):
        try:
            #hm = hmac.new(secret, request.body, hashlib.sha256)
            # add full url (with query string, if any) 
            sig_str = request.build_absolute_uri()
            
            # first need to sort all POST data by paramer
            post_params = []
            for key in request.POST:
                #print(key)
                post_params.append(key)
            #print 'post_params', post_params
            post_params.sort()
            for key in post_params:
                sig_str+=key
                sig_str+=request.POST.get(key)
            #print 'sig_str', sig_str
    
             
            #print request.META.get("HTTP_X_TWILIO_SIGNATURE", None) 
            # construct signature
            hm = hmac.new(settings.TWILIO_AUTH_TOKEN, sig_str, hashlib.sha1)
            digest = base64.b64encode(hm.digest())
            # compare it to twilio x-twilio-signature header 
            if request.META.get("HTTP_X_TWILIO_SIGNATURE", None) != digest:
                return HttpResponse(status=401)
        except Exception:
            logger.exception('message')

        return view(request, *args, **kwargs)
 
    return inner

# decorator to check allowed ips (easier and should be faster...)
# not used, to be considered ...
allowedIps = ['129.0.0.1', '127.0.0.1']
def allow_by_ip(view_func):
    def authorize(request, *args, **kwargs):
        user_ip = request.META['REMOTE_ADDR']
        for ip in allowedIps:
            if ip==user_ip:
                return view_func(request, *args, **kwargs)
        return HttpResponse(status=401)
    return authorize



#-------------------------------------------------------------------------------------------------
# test view
#-------------------------------------------------------------------------------------------------
@csrf_exempt
@validate_request_hmac
@require_POST
def index_view(request):
    # Retrieve the request's body and parse it as JSON
# request.META is a dictionary of headers and
# all headers are converted to upper case and get prefix HTTP_
#    print request.META.get('HTTP_AUTOCALLBOT_SIGNATURE') 
# decorator does this
#    hm = hmac.new(secret, request.body, hashlib.sha256)
#    digest = base64.b64encode(hm.digest())
# 
#    if request.META.get('HTTP_AUTOCALLBOT_SIGNATURE', None) != digest:
#        return HttpResponse(status=401)

#    FULL_URL_WITH_QUERY_STRINg:
    #print request.build_absolute_uri() 
    #print type(request.build_absolute_uri()) # string, dah
   
    #print 'HOSTNAME', request.get_host()
 
    event_json = json.loads(request.body) 
    #print event_json['key']
    if event_json['key'] == 'work':
        x = 0.
        for j in range(10**2):
            for i in range(10**4):
                x = x + 2.**2 + (-1)**2 * i 
        #print x    
 
#    # Do something with data here
    return HttpResponse(status=200)






# testline_id formating, insert 'w' every other number
def format_testline_id(testline_id, insert_string):
    if len(insert_string) > 0:
        new = ''
        for letter in testline_id:
            new = new + letter + insert_string

        return new[:-len(insert_string)]
    return testline_id




#------------------------------------------------------------------------------------
# view for making calls
#------------------------------------------------------------------------------------
@csrf_exempt
@validate_request_hmac
@require_POST
def calls_view(request):
    t0 = time.time()
    # get parameters
    #params = {'max_attempts':4, 'call_chunk_size':20, 'twilio_numbers':["+17207070348", "+17204282260"]}
    params = json.loads(request.body) 
    ivr_type = params['ivr_type']
    twilio_numbers = params['twilio_numbers']
    #max_attempts = params['max_attempts']
    max_attempts = Autocall.MAX_ATTEMPTS 
    call_chunk_size = params['call_chunk_size']
    # select calls that are active, have result='noresult', and attempts <  max_attempts
    # and take only the call_chunk_size number of autocalls  (next POST will take another chunk...) 
    # don't need to filter for attempts because recordings_view will set result = Autocall.ERR_ATTEMPTS if NORESULT and attempts = MAX_ATTEMPTS 
    # take only 'sentry' calls
    queryset = Autocall.objects.filter(active=True).filter(result=Autocall.NORESULT).filter(ivr_type=ivr_type)[:call_chunk_size]
    # update queryset with result=CALLING, so that next POST doesn't query these
    #queryset.update(result=Autocall.CALLING, result_timestamp=timezone.now())
    # doesn't work, have to for loop it
    for autocall in queryset:
        if autocall.is_time_ok: #or autocall.is_time_ok_test:
        # only the ones that will actually call
            autocall.result = Autocall.CALLING
            autocall.result_timestamp = timezone.now()
            autocall.save()
    autocall_counter = 0
    for autocall in queryset:
        #print 'call:', autocall, '--------------------------------------------------------------------' # dd
        # time check is done by properties in Autocall class
        #if time_now > autocall.time and time_now > autocall.t_start and time_now < autocall.t_end:
        # remove is_time_ok_test for production   dd
        if autocall.is_time_ok: # or autocall.is_time_ok_test:
        # time ok
            # first select a twilio number (alternate through the list of numbers)
            #print 'calling...' # dd
            twilio_number = twilio_numbers[autocall_counter % len(twilio_numbers)]           
            # fire a call 
            try:
                call = client.calls.create(
                    to = autocall.testline_phone,
                    #to = '+17207070348', # test, call to incoming/ and record
                    from_ = twilio_number,
                    #from_ = twilio_numbe, #test, dd
                    #from_ = '+13034434799', # test, dd
                    url = settings.HOST + "process/returnxml/?testlineid=" + format_testline_id(autocall.testline_id, 'w'),
                    # this type of error will be logged by root logger as 400 Not Found, also twilio will email
                    #url = settings.HOST + "rocess/returnxml/?testlineid=" + format_testline_id(autocall.testline_id, 'w'),
                    method="GET",
                    timeout = 40,
                    record=True,
                    #status_callback="http://autocallbot.com/twiliopost/statuscallback/",
                    recording_status_callback= settings.HOST + "process/receiverec/"
                )
            except TwilioRestException, e:
                #print 'e.status', e.status # (int) http status, 
                #print 'e.msg', e.msg # (str) human readable message 
                #print 'e.uri', e.uri # (str) uri that caused exception
                #print 'e.code', e.code #(int/None) twilio code fo the error, not available for all errors
                
                # although e.status and e.code are integers, %s works just fine...
                logger.error('twilio: user = %s, status = %s, uri = %s, code = %s, msg = %s',\
                             autocall.user, e.status, e.uri, e.code, e.msg  )
                # if call didn't go out, set result=NORESULT 
                autocall.result = Autocall.NORESULT
                autocall.result_timestamp = timezone.now()
                #autocall.attempts = autocall.attempts + 1 # don't need this - call did't go out, no charge, let it try again
                autocall.save()
                # looger will send an email...
                #send_mail('twilio error', 'in process/calls/ twilio error code:'+ str(e.code), 'service@autocallbot.com',['auto.callbot@gmail.com'],) 
                #return HttpResponse(status=500)
                continue
            except Exception, e:
                logger.exception("message") # exception traceback gets logged
                # if call didn't go out, set result=NORESULT 
                autocall.result = Autocall.NORESULT
                autocall.result_timestamp = timezone.now()
                #autocall.attempts = autocall.attempts + 1 # don't need this - call did't go out, no charge, let it try again
                autocall.save()
                continue

            # no exceptions ....
            #print call.sid # dd
            autocall_counter = autocall_counter + 1
            # update autocall instance 
            autocall.twilio_call_sid = call.sid
            autocall.attempts = autocall.attempts + 1
            autocall.save()

    t1 = time.time()
    #print 'total time ellapsed: ', t1-t0
    #print 'total calls made ', autocall_counter
    #return HttpResponse(status=200)
    response = JsonResponse({'total calls made' : autocall_counter})

#    return HttpResponse(status=200)
    return response


#--------------------------------------------------------------------------------------------------
# view to process recordings ----------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------

@csrf_exempt
@validate_request_hmac
@require_POST
def recordings_view(request):
    # get parameters
    t0 = time.time()

    params = json.loads(request.body) 

    rec_chunk_size = params['rec_chunk_size']
    #print 'recordings_view_----------------------------------------------------------'
    #print 'rec_chunk_size:', rec_chunk_size
    # query autocalls: active=True, result='ready2process' ('ready2process' is set by receiverec_view after POST from Twilio) 
    # and take only the rec_chunk_size number of autocalls  (next POST will take another chunk...) 
    autocall_counter = 0
    queryset = Autocall.objects.filter(active=True).filter(result=Autocall.READY2PROCESS)[:rec_chunk_size]
    # update this queryset with result='processing', so that next POST doesn't mess with this chunk
    #queryset.update(result=Autocall.PROCESSING, result_timestamp=timezone.now())
    # this does not work because 'Cannot update a query once a slice has been taken.'
    # workaround with subquery also doesn't, because "This version of MySQL does not support LIMIT & IN/ALL/ANY/SOME subquery"
    # so just for loop 
    t2 = time.time()
    for autocall in queryset:
        autocall.result = Autocall.PROCESSING
        autocall.result_timestamp = timezone.now()
        autocall.save()
    t3 = time.time()
    #print 'time for updating queryset: ', t3-t2 
    t_download = 0.
    for autocall in queryset:
        try:
            #print 'autocall ', autocall
            recording_sid = autocall.twilio_rec_sid 
            #print recording_sid
            # get recording.json to obtain duration 
            recording_url = 'https://api.twilio.com/2010-04-01/Accounts/' + settings.TWILIO_ACCOUNT_SID + '/Recordings/'+recording_sid
            rj = requests.get(recording_url+'.json')
            #call_sid = rj.json()['call_sid']
            duration = int(rj.json()['duration']) 
            #print 'duration:', duration
            # process recording
            if duration < 10:
            # if duration < 10, dont process, log error, result = noresult
                autocall_counter = autocall_counter + 1
                #print 'duration < 10' 
                result = Autocall.NORESULT
                # log an error ....
                logger.error('duration = %d, user = %s, rec_sid= %s ', duration,  autocall.user, autocall.twilio_rec_sid )
                #send_mail('error:duration','user='+autocall.user.__str__(), 'service@autocallbot.com', ['auto.callbot@gmail.com'])
    
            elif duration > 10:
                # download recording
                #print 'downloading rec...'
                t4=time.time()
                r = requests.get(recording_url)
                recording_file_path = settings.RECORDINGS_DIR + recording_sid 
                with open(recording_file_path, 'wb') as f:
                    f.write(r.content)
                    f.close()
                t5=time.time()
                t_download = t_download + t5-t4
                # process recording, result can only be one of the following strings: 'donot', 'req', 'exp', 'ver' 'noresult'
                fprints_file_path = settings.FPRINTING_DIR + 'obj/' + 'fprints_lib3.pkl'
                result = Autocall.NORESULT 
                #print 'processing ...'
                try:
                    result = process_recording(fprints_file_path, recording_file_path) 
                except:
                    logger.error('process_recording script exception')
                    pass
                # log 'noresult' 'ver'
                if result == Autocall.NORESULT or result == Autocall.VER:
                    # there might be a lot of these, hence 'info', mail_admins should only email on 'warning' and above. 
                    # we'll log 'error' only on NORESULT after MAX_ATTEMPTS 
                    logger.info('user = %s, result = %s, rec_sid= %s ',autocall.user, result, recording_sid )
                    
                #print('result= '+result)
                # after logging, change 'ver' to 'noresult', so that call can go out again
                # most of the time 'ver' result is false (DTNS didn't catch the id correctly)
                if result == Autocall.VER:
                    result = Autocall.NORESULT
                #delete wav file from disk?                                                                                                           
                os.remove(settings.RECORDINGS_DIR + recording_sid)
                autocall_counter = autocall_counter + 1
                
            # update autocall 
            autocall.result = result
            if autocall.attempts == Autocall.MAX_ATTEMPTS and autocall.result == Autocall.NORESULT:
            # if MAX_ATTEMPTS attempts  and 'noresult'
                autocall.result = Autocall.ERR_ATTEMPTS
                logger.error('user = %s, attempts = %s, result = %s', autocall.user, autocall.attempts, autocall.result)
                #send_mail('error:attempts','user='+autocall.user.__str__(), 'service@autocallbot.com', ['auto.callbot@gmail.com'])
            # set the time stamp for result
            autocall.result_timestamp = timezone.now()
            autocall.save()

        except Exception, e:
            #something went wrong, log, set result=NORESULT(so that next calls_view POST can get to it) and continue to the next autocall of the chunk
            logger.exception('message')
            autocall.result = Autocall.NORESULT
            autocall.result_timestamp = timezone.now()
            autocall.save()
            continue 
 
    t1 = time.time()
    #print 'total time ellapsed: ', t1-t0
    #print 't_download: ', t_download
    #print 'total recordings processed: ', autocall_counter
    #return HttpResponse(status=200)
    response = JsonResponse({'total recordings processed' : autocall_counter})
    return response


#---------------------------------------------------------------------------------------------------
# view to return xml directions for a call
#---------------------------------------------------------------------------------------------------

# twiml directions, returns xml
# needs formatted testlineid as GET parameter
# i.e. get process/xmlreturn/?testlineid=1w2ww3w4
# could be done differently, as with Agg calls (send autocall pk as GET param) 
@require_GET
def returnxml_view(request):
    # some kind of auth ...
    # implement twilio signatures later, need https for that...
    if request.GET.get('AccountSid') == settings.TWILIO_ACCOUNT_SID:
        return render(request,'process/call_directions.xml')
    else:
        return HttpResponse(status=401)




#---------------------------------------------------------------------------------------------------
# view for recording_callback
#---------------------------------------------------------------------------------------------------
 
# receive recording sid as POST and add it to the appropriate call
@csrf_exempt
@validate_twilio_post
@require_POST
def receiverec_view(request):
 
    recording_sid = request.POST['RecordingSid']
    call_sid = request.POST['CallSid']
    #print 'receiverec_view---------------------------------------------' 
    #print 'RecordingSid', recording_sid 
    #print 'CallSid', call_sid 
    # call_sid id unique
    autocall = Autocall.objects.get(twilio_call_sid = call_sid)
    autocall.twilio_rec_sid = recording_sid
    autocall.result = Autocall.READY2PROCESS 
    autocall.result_timestamp = timezone.now()
    autocall.save()
    return HttpResponse(status=200)




#---------------------------------------------------------------------------------------------------
# view for deactivation of autocalls whose user's have expired access_until 
# selects users that have free_trial or one_time5  or monthly4 subscriptions, have active call (not currently processing)
# deactivates the call, sets result = NORESULT_INACTIVE
# stripehook (upon customer.subscription.deleted) only sets access_until, then this script will take care of the rest
#---------------------------------------------------------------------------------------------------

@csrf_exempt
@validate_request_hmac
@require_POST
def deactivate_view(request):
    # get parameters
    params = json.loads(request.body) 

    deactivate_chunk_size = params['deactivate_chunk_size']
    #print 'deactivate_view_----------------------------------------------------------'
    #print 'deactivate_chunk_size:', deactivate_chunk_size
    # query users that have expired access_until
    # and take only the deactivate_chunk_size number of them  (next POST will take another chunk...)
    # although I think I could deactivate  'all' ... for now chunk size is 1000 - I'm an optimist!
    autocall_counter = 0
    email_list = []
    queryset = User.objects.filter(access_until__lt=timezone.now())[:deactivate_chunk_size]
    for user in queryset:
        if hasattr(user, 'autocall'):
        # if user has autocall
            autocall = user.autocall
            if autocall.active == True and not autocall.is_processing:
            # if autocall is active and autocall is not 'processing' 
                autocall.active = False
                autocall.result = Autocall.NORESULT_INACTIVE
                autocall.result_timestamp = timezone.now() # don't really need this, dashboard displays todays date when call is inactive
                autocall.save()

                # notify user
                host = settings.HOST + 'dashboard/'  
        
                subject = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/emails/subscription_expired_subject.txt', {})
                # subject can't contain new lines
                subject = ''.join(subject.splitlines())
                body = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/emails/subscription_expired_body.txt', {'host':host})

                #print 'body------------------------------------------------------'
                #print  body
                # a tuple or list of recipient email addresses
                to = [autocall.user.email,]
                #to = ['ignas.satkauskas@gmail.com',]
                # from email address, if omitted, DEFAULT_FROM_EMAIL is used
                #from_email = 'noreply@autocallbot.com'
                #from_email = 'Autocallbot <noreply@autocallbot.com>' # does not work with send_mail, but does with EmailMessage
                from_email = settings.NOTIFICATIONS_FROM_EMAIL # does not work with send_mail, but does with EmailMessage
                # email object 
                email = EmailMessage(
                   subject=subject,
                   body=body,
                   from_email=from_email,
                   to=to,
                )
                autocall_counter = autocall_counter + 1
                #email.send() # each email will open and close SMTP connection 
                # add email to list to be sent later with the rest of the chunk in one SMTP connection 
                email_list.append(email) 
                
                # also, send text message if user has mobile # set up
                if user.notify_by == User.TEXT or user.notify_by == User.BOTH: 
                
                    twilio_number = settings.TWILIO_FROM_NUMBER
                    # AUTOCALLBOT call on Fri Mar  2 10:51:40 2018
                    #text_body = 'AUTOCALLBOT: Your subscription has exprired, automatic call is now inactive. You will no longer receive any notifications.'   
            
                    text_body = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/texts/subscription_expired.txt', {})
                    #print 'len(text_body)', len(text_body)
                    #print 'body--------------------------------------------' 
                    #print text_body
                    try:
                        # send message
                        message = client.messages.create(
                            to = autocall.user.mobile_number,
                            from_ = twilio_number,
                            body = text_body 
                        ) 
                    except TwilioRestException, e:
                        #print 'e.status', e.status # (int) http status, 
                        #print 'e.msg', e.msg # (str) human readable message 
                        #print 'e.uri', e.uri # (str) uri that caused exception
                        #print 'e.code', e.code #(int/None) twilio code fo the error, not available for all errors
                        
                        # although e.status and e.code are integers, %s works just fine...
                        logger.error('twilio: user = %s, status = %s, uri = %s, code = %s, msg = %s',\
                                     autocall.user, e.status, e.uri, e.code, e.msg  )
                        # looger will send an email...
                        #send_mail('twilio error', 'in notifications/ twilio error code:'+ str(e.code), 'service@autocallbot.com',['auto.callbot@gmail.com'],) 
                        #return HttpResponse(status=500)
                        continue
                    except Exception, e:
                        logger.exception("message") # exception traceback gets logged
                        continue
       

    # send out a bunch of emails 
    try:
        connection = mail.get_connection()   # Use default email connection as specified in EMAIL_BACKEND
        connection.send_messages(email_list)
        # send_messages opens connection, sends messages, then closes the connection.
        # however, send_messages() will not manually open or close the connection if it is already open
    except Exception, e:
        logger.exception('message')
     
    #print 'total autocalls deactivated: ', autocall_counter
    #return HttpResponse(status=200)
    response = JsonResponse({'total autocalls deactivated' : autocall_counter})
    return response



#---------------------------------------------------------------------------------------------------
# prepwork to be done at the begging of new day:
# set result = NORESULT, emailed = False, texted = False, attempts=0 
#---------------------------------------------------------------------------------------------------
@csrf_exempt
@validate_request_hmac
@require_POST
def prepwork_view(request):
    #t0=time.time()
    # get parameters
    params = json.loads(request.body) 

    prepwork_timezone = params['prepwork_timezone']
    # prepwork_timezone can be 'MST', 'EST', 'CST' or 'PST'
    if prepwork_timezone == 'MST':
        queryset = Autocall.objects.filter(active=True).filter(call_timezone=Autocall.MST)
    elif prepwork_timezone == 'EST':
        queryset = Autocall.objects.filter(active=True).filter(call_timezone=Autocall.EST)
    elif prepwork_timezone == 'CST':
        queryset = Autocall.objects.filter(active=True).filter(call_timezone=Autocall.CST)
    elif prepwork_timezone == 'PST':
        queryset = Autocall.objects.filter(active=True).filter(call_timezone=Autocall.PST)

    #print queryset.count()
    # update autocalls
    # update() method is direct SQL statement
    # returns number of updated rows (can be smaller that queryset if values don't change ... but not really...)
    number_of_updates = queryset.update( result=Autocall.NORESULT, emailed=False, texted=False, attempts=0 )

    # adjust call times for agg calls only, i.e. return times to default values 
    for autocall in queryset:
        if autocall.ivr_type == Autocall.AGGREGATE:
            autocall.mon_time = autocall.mon_time_default
            autocall.tue_time = autocall.tue_time_default
            autocall.wed_time = autocall.wed_time_default
            autocall.thu_time = autocall.thu_time_default
            autocall.fri_time = autocall.fri_time_default
            autocall.sat_time = autocall.sat_time_default
            autocall.sun_time = autocall.sun_time_default
            
            autocall.save()
            
    #t1=time.time()
    #print 'prepwork time: ', t1-t0
    #print 'number_of_updates', number_of_updates
    #print 'type(number_of_updates)', type(number_of_updates) # its 'long'

    response = JsonResponse({'number of updates' : number_of_updates})
    return response


#------------------------------------------------------------------------------------
# view for making custom calls including aggregate 
# process/calls/custom/
#------------------------------------------------------------------------------------
@csrf_exempt
@validate_request_hmac
@require_POST
def calls_custom_view(request):
    t0 = time.time()
    # get parameters
    #params = {'max_attempts':4, 'call_chunk_size':20, 'twilio_numbers':["+17207070348", "+17204282260"]}
    params = json.loads(request.body) 
    ivr_type = params['ivr_type']
    twilio_numbers = params['twilio_numbers']
    #max_attempts = params['max_attempts']
    max_attempts = Autocall.MAX_ATTEMPTS 
    call_chunk_size = params['call_chunk_size']
    # select calls that are active, have result='noresult', and attempts <  max_attempts
    # and take only the call_chunk_size number of autocalls  (next POST will take another chunk...) 
    # don't need to filter for attempts because recordings_view will set result = Autocall.ERR_ATTEMPTS if NORESULT and attempts = MAX_ATTEMPTS 
    # take only 'sentry' calls
    queryset = Autocall.objects.filter(active=True).filter(result=Autocall.NORESULT).filter(ivr_type=ivr_type)[:call_chunk_size]
    # update queryset with result=CALLING, so that next POST doesn't query these
    #queryset.update(result=Autocall.CALLING, result_timestamp=timezone.now())
    # doesn't work, have to for-loop it
    for autocall in queryset:
        if autocall.is_time_ok_custom: #or autocall.is_time_ok_test:
        # only the ones that will actually call
            autocall.result = Autocall.CALLING
            autocall.result_timestamp = timezone.now()
            autocall.save()
    autocall_counter = 0
    for autocall in queryset:
        #print 'call:', autocall, '--------------------------------------------------------------------' # dd
        # time check is done by properties in Autocall class
        #if time_now > autocall.time and time_now > autocall.t_start and time_now < autocall.t_end:
        # remove is_time_ok_test for production   dd
        if autocall.is_time_ok_custom: # or autocall.is_time_ok_test:
        # time ok
            # first select a twilio number (alternate through the list of numbers)
            twilio_number = twilio_numbers[autocall_counter % len(twilio_numbers)]           
            # make query string for GET key=value
            GET_string = ''
            GET_string+= 'IvrType='+ivr_type
            # don't take pound sign at the end, don't want to code it as %23 ...
            # pass autocalls pk so that process/directions/ knows what to return to Twilio
            # should do this to aggregate calls as well, and sentry calls... I'm so stupid!
            GET_string+= '&AutocallPk=' + str(autocall.pk) 

 
            # fire a call 
            try:
                call = client.calls.create(
                    to = autocall.testline_phone,
                    #to = '+17207070348', # test, call to incoming/ and record
                    from_ = twilio_number,
                    #from_ = twilio_numbe, #test, dd
                    #from_ = '+13034434799', # test, dd
                    #url = settings.HOST + "process/returnxml/?testlineid=" + format_testline_id(autocall.testline_id, 'w'),
                    url = settings.HOST + "process/directions/?" + GET_string ,
                    # this type of error will be logged by root logger as 400 Not Found, also twilio will email
                    #url = settings.HOST + "rocess/returnxml/?testlineid=" + format_testline_id(autocall.testline_id, 'w'),
                    method="GET",
                    timeout = 40,
                    # record the whole call, for testing purposes only
                    #record = True,
                    # webhook with POST parameter CallStatus='completed' will be sent there
                    # also will have CallDuration parameter
                    # if call duration is about 17s (vs ~47), it means "Schedule is not currently available" 
                    status_callback = settings.HOST +  "process/statuscallback/", 
                    #status_callback="http://autocallbot.com/twiliopost/statuscallback/",
                    #recording_status_callback= settings.HOST + "process/receiverec/"
                )
            except TwilioRestException, e:
                #print 'e.status', e.status # (int) http status, 
                #print 'e.msg', e.msg # (str) human readable message 
                #print 'e.uri', e.uri # (str) uri that caused exception
                #print 'e.code', e.code #(int/None) twilio code fo the error, not available for all errors
                
                # although e.status and e.code are integers, %s works just fine...
                logger.error('twilio: user = %s, status = %s, uri = %s, code = %s, msg = %s',\
                             autocall.user, e.status, e.uri, e.code, e.msg  )
                # if call didn't go out, set result=NORESULT 
                autocall.result = Autocall.NORESULT
                autocall.result_timestamp = timezone.now()
                #autocall.attempts = autocall.attempts + 1 # don't need this - call did't go out, no charge, let it try again
                autocall.save()
                # looger will send an email...
                #send_mail('twilio error', 'in process/calls/ twilio error code:'+ str(e.code), 'service@autocallbot.com',['auto.callbot@gmail.com'],) 
                #return HttpResponse(status=500)
                continue
            except Exception, e:
                logger.exception("message") # exception traceback gets logged
                # if call didn't go out, set result=NORESULT 
                autocall.result = Autocall.NORESULT
                autocall.result_timestamp = timezone.now()
                #autocall.attempts = autocall.attempts + 1 # don't need this - call did't go out, no charge, let it try again
                autocall.save()
                continue

            # no exceptions ....
            #print call.sid # dd
            autocall_counter = autocall_counter + 1
            # update autocall instance 
            autocall.twilio_call_sid = call.sid
            autocall.attempts = autocall.attempts + 1
            autocall.save()

    #t1 = time.time()
    #print 'total time ellapsed: ', t1-t0
    #print 'total calls made ', autocall_counter
    #return HttpResponse(status=200)
    response = JsonResponse({'total calls made' : autocall_counter})

#    return HttpResponse(status=200)
    return response

#---------------------------------------------------------------------------------------------------
# view for status_callback, from custom_call_view fired call
# if call duration is less than 20s, reschedule a call up to 2 times
#---------------------------------------------------------------------------------------------------
# process/statuscallback/ 
@csrf_exempt
@validate_twilio_post
@require_POST
def status_callback_view(request):
 
    #recording_sid = request.POST['RecordingSid']
    call_sid = request.POST['CallSid']
    call_duration = request.POST.get('CallDuration')
    #print('call duration in process/status_callback_view/')
    #print(call_duration)
    if int(call_duration) < 20: 
    # if "your schedule is not currently available"
    # in this case, response.record step is not reached(transcription_view not hit), and no transcription is available
    # add an hour to call time (of the appropriate weekday) up to 2 times (3 attempts total)
    # if still < 20s, "no result, call manually"
        autocall = Autocall.objects.get(twilio_call_sid = call_sid)
        #transcription_text = "Your schedule is not currently available. Please try again later. Thank you, goodbye." 
        if autocall.attempts < 3:
            # which weekday it is in local time?
            timezone.activate(autocall.call_timezone)
            datetime_now = timezone.localtime(timezone.now()) # this is now in users local time
            timezone.deactivate()
            #time_now = datetime_now.time() # just get the datetime.time object, it is naive now, can compare to autocall.time
            weekday_now = datetime_now.weekday()
            # returns integer: 0 for Mon, 1 for Tue, ..., 6 for Sun
            # change call_time to 1 hour later, depending on which weekday it is:
            if weekday_now == 0:
                autocall.mon_time = datetime.time(hour = autocall.mon_time.hour + 1 )
            elif weekday_now == 1:
                autocall.tue_time = datetime.time(hour = autocall.tue_time.hour + 1 )
            elif weekday_now == 2:
                autocall.wed_time = datetime.time(hour = autocall.wed_time.hour + 1 )
            elif weekday_now == 3:
                autocall.thu_time = datetime.time(hour = autocall.thu_time.hour + 1 )
            elif weekday_now == 4:
                autocall.fri_time = datetime.time(hour = autocall.fri_time.hour + 1 )
            elif weekday_now == 5:
                autocall.sat_time = datetime.time(hour = autocall.sat_time.hour + 1 )
            elif weekday_now == 6:
                autocall.sun_time = datetime.time(hour = autocall.sun_time.hour + 1 )
        
            autocall.result = Autocall.NORESULT
            autocall.save()
        
        else:
        # attempts == 3
            autocall.result = Autocall.TRANSCRIPTION
            transcription_text = "No result after 3 call attempts, please call manually."
            autocall.transcription = transcription_text
            autocall.save()


    return HttpResponse(status=200)


#---------------------------------------------------------------------------------------------------
# view to return xml directions for a  aggregate call or custom call
# twiml directions, returns xml
# get parameters are passed by calls_aggregate_view
# i.e. get process/directions/?IvrType=aggregate&EnterKeys1=1w2ww3w4
#  get process/directions/?IvrType=custom&AutocallPk=26
#---------------------------------------------------------------------------------------------------

@require_GET
def directions_view(request):
    # some kind of auth ...
    # implement twilio signatures later, need https for that...
        print request.GET.get('IvrType')
#    if request.GET.get('AccountSid') == settings.TWILIO_ACCOUNT_SID:
        

        if request.GET.get('IvrType') == 'custom' or request.GET.get('IvrType') == 'aggregate':
            # first get the AutocallPk from GET params
            autocall_pk = long(request.GET.get('AutocallPk'))
            # retrieve autocall by pk
            autocall = Autocall.objects.get(pk=autocall_pk)
            # now get call directions
            number_of_steps = autocall.number_of_steps # don't need this, just use @property call_steps
            transcribe_seconds = autocall.transcribe_seconds 

            # action url for record verb in twiML,
            # after recording is completed, call will be directed to this url i.e. even hangup verb is unreachable after Record verb
            # only if Record receives empty recording, it continues disregards action url and continues with next verb
            action_url = settings.HOST + "process/directions/hangup/" 
            # transcription of the recording will be sent there once available
            # asynchronous POST with TranscriptionText as one of the params 
            transcribe_url = settings.HOST + "process/transcription/"



            response = VoiceResponse()

            for step in autocall.call_steps:
            #step = (promt, enter_keys, wait_time)
                response.pause(length = step[2]) 
                response.play(digits = step[1])

            # for testing purposes, short pause when recording the whole call, instead of response.record below
            #response.pause(length=10)
                
            response.record(action = action_url,
                            method = 'GET', 
                            # timeout after 30s of silence, stop recording
                            timeout = 30,
                            # max_length of the recording
                            max_length = transcribe_seconds,
                            # don't need beep sound before recording (default is True)
                            play_beep = False,
                            # do not trim silence
                            #trim = 'do-not-trim',
                            # if transcribeCallback is given, transcribe=True is implied
                            transcribe_callback = transcribe_url, 
            )

            response.hangup()
            # this verb will be reached only if recording is empty

            return HttpResponse(str(response))


        else:
            return HttpResponse(status=401)






#---------------------------------------------------------------------------------------------------
# view to return xml directions for a call to hangup
#---------------------------------------------------------------------------------------------------

# twiml directions, returns xml @
# /process/directions/hangup/
@require_GET
def hangup_view(request):

    response = VoiceResponse()
    response.hangup()

    return HttpResponse(str(response))



#---------------------------------------------------------------------------------------------------
# view for transcribeCallback used in Record verb of directions_view
# /process/transcription/
# receive recording sid as POST and add it to the appropriate call
# receive transcriptionText as well
#---------------------------------------------------------------------------------------------------


@csrf_exempt
@validate_twilio_post
@require_POST
def transcription_view(request):

    transcription_status = request.POST.get('TranscriptionStatus')
    if transcription_status == 'completed':
    # it can also be 'failed', in rare cases then recording is 1s long... 
    # if it is 'failed', then transcription_text will be null and fix_twilio_transcription will choke
        recording_sid = request.POST.get('RecordingSid')
        call_sid = request.POST.get('CallSid')
        transcription_text = request.POST.get('TranscriptionText')
        
        #print 'receiverec_view---------------------------------------------' 
        #print 'RecordingSid', recording_sid 
        #print 'CallSid', call_sid 
        # call_sid id unique
        autocall = Autocall.objects.get(twilio_call_sid = call_sid)
        autocall.twilio_rec_sid = recording_sid
        autocall.transcription = transcription_text 
        autocall.result_timestamp = timezone.now()
        autocall.save()

   
        # check is transcription is valid
        # i.e. will check if word 'password' is in the transcription
        if autocall.is_twilio_transcription_valid:
    
            # fix transcription (it saves autocall after fixing) 
            autocall.fix_twilio_transcription

            autocall.result = Autocall.TRANSCRIPTION 
            autocall.save()

        else:
        # transcription is not valid
            if autocall.attempts < 3:
            # if its not 3rd attempt, set NORESULT so that it calls again
                autocall.result = Autocall.NORESULT
                autocall.result_timestamp = timezone.now()
                autocall.save()
                logger.error('twilio transcription invalid: user = %s, call sid = %s, attempts = %s',\
                             str(autocall.user), call_sid, str(autocall.attempts) )
            else:
            # 3rd attempt
                autocall.result = Autocall.TRANSCRIPTION
                autocall.result_timestamp = timezone.now()
                transcription_text = "No result after 3 call attempts, please call manually."
                autocall.transcription = transcription_text
                autocall.save()

    
#    elif transcription_status == 'failed':
#    # this need more attention
#    # one special case: May 21, recording was 1s long, hence no transcription,
#    # call length was 18s, so statuscallback trigerred NORESULT and rescheduled it to 1h later
#    # but which one came first: statuscallback_url or transcription_url POST? (need attention)
#        autocall.result = Autocall.NORESULT 
#        autocall.result_timestamp = timezone.now()
#        autocall.save()
        
        
     
    return HttpResponse(status=200)



#---------------------------------------------------------------------------------------------------
# for testing purposes. record incoming call.
# webhook for call to (720) 707-0348
#---------------------------------------------------------------------------------------------------

#@require_GET
#def incoming_view(request):
#    """Returns TwiML which prompts the caller to record a message"""
#    # Start our TwiML response
#    response = VoiceResponse()
#    # Use <Say> to give the caller some instructions
#    response.say('Hello.')
#    # Use <Record> to record the caller's message
#    response.record(method="GET", timeout=30, maxLength=60, finishOnKey='#')
#    #print response
#    return HttpResponse(response)

# don't need the above, just return incoming.xml with directions 

    #return render(request,'process/incoming.xml')



