
from twilio.rest import Client
import importlib
import sys

def load_twilio_config():
    sys.path.append('../config')
    config_file = 'cfg_twilio'
    ctw = importlib.import_module(config_file)
    return (ctw.twilio_number, ctw.twilio_account_sid, ctw.twilio_auth_token, ctw.to_number)


class MessageClient(object):
    def __init__(self):
        (twilio_number, twilio_account_sid,
         twilio_auth_token, to_number) = load_twilio_config()

        self.twilio_number = twilio_number
        self.twilio_client = Client(twilio_account_sid,
                                              twilio_auth_token)
        self.to_number = to_number

    def send_message(self, body):
        self.twilio_client.messages.create(body=body, to=self.to_number,
                                           from_=self.twilio_number,
                                           )