import logging                                                                                                                                                                                 
# Get an instance of a logger
logger = logging.getLogger(__name__)

from django.contrib import admin
from .models import Autocall
from .models import UsedAutocall
from django.core.mail import EmailMessage
from django.core import mail # needed for opening SMTP connection
from django.conf import settings
from django.contrib import messages
from django.template.loader import render_to_string

from django.urls import reverse
from django.utils.html import format_html

from django import forms

#twilio
from twilio.rest import Client
from twilio.base.exceptions import TwilioRestException
from twilio.twiml.voice_response import VoiceResponse
client = Client(settings.TWILIO_ACCOUNT_SID, settings.TWILIO_AUTH_TOKEN)


# custom actions -------------------------------------------------
# noresult
def set_noresult(modeladmin, request, queryset):
    queryset.update(result=Autocall.NORESULT)
set_noresult.short_description = "Set 'noresult'"

# calling 
def set_calling(modeladmin, request, queryset):
    queryset.update(result=Autocall.CALLING)
set_calling.short_description = "Set 'calling'"

# ready2process 
def set_ready2process(modeladmin, request, queryset):
    queryset.update(result=Autocall.READY2PROCESS)
set_ready2process.short_description = "Set 'ready2process'"

# processing
def set_processing(modeladmin, request, queryset):
    queryset.update(result=Autocall.PROCESSING)
set_processing.short_description = "Set 'processing'"

# zero out attempts
def zero_out_attempts(modeladmin, request, queryset):
    queryset.update(attempts=0)
zero_out_attempts.short_description = "Zero out attempts"

# activate
def activate(modeladmin, request, queryset):
    queryset.update(active=True)
activate.short_description = "Activate"

# deactivate
def deactivate(modeladmin, request, queryset):
    queryset.update(active=False)
deactivate.short_description = "Deactivate"

# clear call sid
def clear_call_sid(modeladmin, request, queryset):
    queryset.update(twilio_call_sid='')
clear_call_sid.short_description = "Clear CallSids"

# clear rec sid
def clear_rec_sid(modeladmin, request, queryset):
    queryset.update(twilio_rec_sid='')
clear_rec_sid.short_description = "Clear RecSids"

# clear rec sid
def set_not_emailed(modeladmin, request, queryset):
    queryset.update(emailed=False)
set_not_emailed.short_description = "Set not emailed"



# Register your models here.

class AutocallAdmin(admin.ModelAdmin):

    # override this method so that transcription CharField of Autocallbot model
    # is displayed with TextArea widget, on admin site only

    def formfield_for_dbfield(self, db_field, **kwargs):
        formfield = super(AutocallAdmin, self).formfield_for_dbfield(db_field, **kwargs)
        if db_field.name == 'transcription':
            formfield.widget = forms.Textarea(attrs=formfield.widget.attrs)
        return formfield

    list_display = ('active', 
                    'is_access_expired', 
                    #'user',
                    'link_to_user',
                    'get_ivr_type',
                    'testline_phone',
                    'testline_id', 
                    'last_name', 
                    'result',
                    'attempts',
                    't_start',
                    't_end',
                    'notify_by',
                    'emailed',
                    'texted', 
                    'get_twilio_call_sid',
                    'get_twilio_rec_sid', 
    )
    # fields that will link to change view, if none provided, defaults to the first element of the list_display tuple
    list_display_links = ('active', 'testline_phone', )
    #list_editable = ('result',)
    list_filter =('active','result','ivr_type')


    # if action is written as method, use quotes
    actions = [activate,
                deactivate, 
                'email_err_params',
                'email_err_params_agg',
                'text_err_params',
                'text_err_params_agg',
                'email_err_server',
                'text_err_server',
                'email_verified_params',
                'text_verified_params',
                'email_user_trial_used',
                'email_call_trial_used',
                'email_mobile_trial_used',
                'email_acb_is_back',
                'copy_to_used',
                set_noresult,
                set_calling, 
                set_ready2process, 
                set_processing, 
                zero_out_attempts,
                clear_call_sid, 
                clear_rec_sid, 
                set_not_emailed,
    ]

    #search_fields = ('result', )

    # make user@email.com in list display link to user change form
    def link_to_user(self, obj):
        link = reverse("admin:accounts_user_change", args=[obj.user.id])
        #return format_html('<a href="{}">Edit {}</a>', link, obj.user) if obj.user else None
        return format_html('<a href="{}">Edit {}</a>', link, obj.user) if obj.user else None
    link_to_user.short_description = 'User'



    # customize list_display callables
    # shorten ivr_type 
    def get_ivr_type(self, obj):
        if obj.ivr_type == 'sentry':
            name = 'sentry' 
        elif obj.ivr_type == 'aggregate':
            name = 'agg'
        elif obj.ivr_type == 'custom':
            name = 'custom'
        return name 
    get_ivr_type.short_description = 'IVR'

    # display only part of twilio sids 
    def get_twilio_call_sid(self, obj):
        return  obj.twilio_call_sid[:5]+'...'

    get_twilio_call_sid.short_description = 'CallSid'


    def get_twilio_rec_sid(self, obj):
        return  obj.twilio_rec_sid[:5]+'...'

    get_twilio_rec_sid.short_description = 'RecSid'

    def is_access_expired(self,obj):
        if obj.user.is_access_expired:
            return 'Expired'
        else:
            return 'OK'
    is_access_expired.short_description = 'Access'

    def notify_by(self,obj):
        return obj.user.notify_by
    notify_by.short_description = 'Notify'


# ACTIONS AS METHODS
# we'll write this action as method (self instead of modeladmin, move it to within class, and use quotes in actions list)
    def email_notification(self, request, queryset, subjectfile, bodyfile):
        #send emails
        email_list = []
        for autocall in queryset:
            #create context for use in rendering email templates(not all emails need it, but some do)
            context = {'user': autocall.user.__str__(),
                       'mobile_number': autocall.user.mobile_number,
                       # it could be empty string (default)
                       'testline_phone': autocall.testline_phone,
                       'testline_id': autocall.testline_id,
                       'lastname': autocall.last_name,
                       # for aggregate calls
                       'enter_keys1': autocall.enter_keys1,
                       'enter_keys2': autocall.enter_keys2,
            }
            user = autocall.user
            if user.notify_by == user.EMAIL or user.notify_by == user.BOTH:

                #subject = 'Error: incorrect call parameters'
                subject = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/emails/' + subjectfile,{} )
                # subject can't contain new lines
                subject = ''.join(subject.splitlines())
                body = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/emails/' +  bodyfile, context)

                # a tuple or list of recipient email addresses
                to = [user.email,]
                from_email = settings.NOTIFICATIONS_FROM_EMAIL # does not work with send_mail, but does with EmailMessage
                # email object 
                email = EmailMessage(
                   subject=subject,
                   body=body,
                   from_email=from_email,
                   to=to,
                )
         
                email_list.append(email)
        
        try: 
            connection = mail.get_connection()   # Use default email connection as specified in EMAIL_BACKEND
            delivered_emails = connection.send_messages(email_list)
            self.message_user(request, "%s emails sent successfully." % delivered_emails )
        except Exception, e:
            self.message_user(request, 'Error sending emails.', level = messages.ERROR)
            #messages.error(request, 'Error sending emails.')
            logger.exception("message") 
            

    def email_err_params(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'err_params_subject.txt',
            bodyfile = 'err_params_body.txt',
        )
    email_err_params.short_description = "Email: err params"


    def email_err_server(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'err_server_subject.txt',
            bodyfile = 'err_server_body.txt',
        )
    email_err_server.short_description = "Email: err server"

    def email_verified_params(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'verified_params_subject.txt',
            bodyfile = 'verified_params_body.txt',
        )
    email_verified_params.short_description = "Email: verified params"

# Autocallbot it back! (temporary email, to notify old users)
    def email_acb_is_back(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'acb_is_back_subject.txt',
            bodyfile = 'acb_is_back_body.txt',
        )
    email_acb_is_back.short_description = "Email: acb is back!"
    pass

# if same user email creates a new call after free trial expires(deletes acc then creates it again),
# set accout to expires and send this email
    def email_user_trial_used(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'user_trial_used_subject.txt',
            bodyfile = 'user_trial_used_body.txt',
        )
    email_user_trial_used.short_description = "Email: user trial used "
    pass

    def email_call_trial_used(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'call_trial_used_subject.txt',
            bodyfile = 'call_trial_used_body.txt',
        )
    email_call_trial_used.short_description = "Email: call trial used"
    pass


    def email_mobile_trial_used(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'mobile_trial_used_subject.txt',
            bodyfile = 'mobile_trial_used_body.txt',
        )
    email_mobile_trial_used.short_description = "Email: mobile trial used"
    pass

    def email_err_params_agg(self, request, queryset):
        return self.email_notification(
            request = request,
            queryset = queryset,
            subjectfile = 'err_params_agg_subject.txt',
            bodyfile = 'err_params_agg_body.txt',
        )
    email_err_params_agg.short_description = "Email: err params agg"
    pass


# texting notifications
# action as method (self instead of modeladmin, move it to within class, and use quotes in actions list)
    def text_notification(self, request, queryset, bodyfile):
        #send texts
        autocall_counter = 0 
        for autocall in queryset:
            user = autocall.user
            if user.notify_by == user.TEXT or user.notify_by == user.BOTH:

                #create context for use in rendering text templates(not all tests need it, but some do)
                context = {'user': autocall.user.__str__(),
                           'mobile_number': autocall.user.mobile_number,
                           # it could be empty string (default)
                           'testline_phone': autocall.testline_phone,
                           'testline_id': autocall.testline_id,
                           'lastname': autocall.last_name,
                           # for aggregate calls
                           'enter_keys1': autocall.enter_keys1,
                           'enter_keys2': autocall.enter_keys2,
                }
                text_body = render_to_string(settings.BASE_DIR + '/autocall/templates/autocall/texts/' +  bodyfile, context)

                try:
                    # send message
                    message = client.messages.create(
                        to = autocall.user.mobile_number,
                        from_ = settings.TWILIO_FROM_NUMBER,
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
                    self.message_user(request, 'Error texting: ' + autocall.user.__str__(), level = messages.ERROR)
                    # looger will send an email...
                    continue
                except Exception, e:
                    logger.exception("message") # exception traceback gets logged
                    self.message_user(request, 'Error texting: ' + autocall.user.__str__(), level = messages.ERROR)
                    continue
                
                # text fired ok
                autocall_counter = autocall_counter + 1
              
        self.message_user(request, "%s texts sent successfully." % autocall_counter )



    def text_err_params(self, request, queryset):
        return self.text_notification(
            request = request,
            queryset = queryset,
            bodyfile = 'err_params.txt',
        )
    text_err_params.short_description = "Text: err params"

    def text_err_server(self, request, queryset):
        return self.text_notification(
            request = request,
            queryset = queryset,
            bodyfile = 'err_server.txt',
        )
    text_err_server.short_description = "Text: err server"

    def text_verified_params(self, request, queryset):
        return self.text_notification(
            request = request,
            queryset = queryset,
            bodyfile = 'verified_params.txt',
        )
    text_verified_params.short_description = "Text: verified params"

    def text_err_params_agg(self, request, queryset):
        return self.text_notification(
            request = request,
            queryset = queryset,
            bodyfile = 'err_params_agg.txt',
        )
    text_err_params_agg.short_description = "Text: err params agg"

    def copy_to_used(self, request, queryset):
        exist_count = 0
        copied_count = 0 
        for autocall in queryset:
            #check if such autocall already exists as usedautocall
            if UsedAutocall.objects.filter(testline_phone=autocall.testline_phone, testline_id=autocall.testline_id, last_name=autocall.last_name):
                exist_count = exist_count + 1
            else:
                autocall.copy_to_used()
                copied_count = copied_count + 1

        if copied_count > 0:
            if copied_count == 1:
                self.message_user(request, "Successfully copied 1 autocall." )
            else:
                self.message_user(request, "Successfully copied %s autocalls." % copied_count )
        if exist_count > 0:
            if exist_count == 1:
                self.message_user(request, "1 autocall already exist, did not copy." , level=messages.ERROR )
            else:
                self.message_user(request, "%s autocalls already exist, did not copy." % exist_count, level=messages.ERROR )
            



class UsedAutocallAdmin(admin.ModelAdmin):

    list_display = ('testline_phone', 
                    'testline_id', 
                    'last_name',
                    'created_timestamp', 
                    'deleted_timestamp', 
    )
    # fields that will link to change view, if none provided, defaults to the first element of the list_display tuple
    #list_display_links = ('testline_phone', )
    #list_editable = ('result',)
    #list_filter =()


    search_fields = ('testline_id', )


admin.site.register(Autocall, AutocallAdmin)
admin.site.register(UsedAutocall, UsedAutocallAdmin)


