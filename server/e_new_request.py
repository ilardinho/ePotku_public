# coding=utf-8
"""
defs __init__ and create_request are from new_request.py Potku with changes.

Created on 16.5.2020
Updated on 30.5.2020

ePotku is a web application version of Potku,
a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2020 Minja Hänninen, Ilari Jalli, Ville Kuokkanen,
Pasi Niininen and Tuomas Pitkänen
​
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
​
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
​
You should have received a copy of the GNU General Public License
along with this program (file named 'LICENCE').
"""
__author__ = "Minja Hänninen \n Ilari Jalli \n Ville Kuokkanen \n" \
             "Pasi Niininen \n ""Tuomas Pitkänen"
__version__ = "1.0"

import configparser
import logging
import time
import os
import re
import shutil

from pathlib import Path
from potku.modules.request import Request
import json
import uuid

from potku.modules.global_settings import GlobalSettings


# TODO: Use Request directly? This was originally for Potku UI dialogs
class NewRequest:
    """Class for creating a new request.
    """

    def __init__(self):
        # TODO: GlobalSettings should be a per-user thing. This is a
        #        temporary fix to prevent overwriting Potku files.
        self.settings = GlobalSettings(Path.cwd() / "testfiles")

        # TODO this path needs to be fixed. Path needs to be formed with user
        #      info which should be provided within the POST-method
        # TODO: change your own path/server path here
        # TODO: MAYBE THIS PATH SHOULD BE SET AS AN ENVIRONMENT VARIABLE
        self.folder = os.path.abspath("testfiles")
        self.directory = None
        self.tab_widgets = {}


    def create_request(self, rqname, request_owner, published):
        """Create new request.
        """
        try:
            # Adding .potku gives all requests the same ending.
            directory = Path(self.folder, f"{rqname}.potku")
            if not directory.exists():
                os.makedirs(directory)
                self.directory = directory
                self.request = Request(self.directory, rqname, self.settings,
                                       self.tab_widgets)

                # adds user's request to json-file
                self.request_json_helper(rqname, request_owner, published)
                print('Created a new request: ' + rqname + '.request')
                logging.getLogger("request").info("Created the request.")
            else:
                print(
                    'Seems like there is another request named as ' + rqname + '. Pick a new name for your request.')
                return
        except Exception as err:
            print(err)
            print("We've done something wrong. "
                  "Most likely invalid request name.")


    @staticmethod
    def request_json_helper(request_name, request_owner, published):
        """ Add's request to all_requests.json -file
        """
        # TODO: Change this to use database
        # Create UUID for the request
        # This creates string out of uuid as python dict doesn't like uuid's
        request_id = str(uuid.uuid1())

        # Wheter the request should be published or not. 1 for public,
        # 0 for non-public
        request_published = published

        request_info = {'requestUUID': request_id, 'requestName': request_name,
                        'requestOwner':
                            request_owner,
                        'published':
                            request_published}

        # TODO: change here the correct path of the JSON-file
        # TODO: path of the JSON-file should be received as an argument
        # TODO: MAYBE THIS PATH SHOULD BE SET AS AN ENVIRONMENT VARIABLE
        requests_file = os.path.abspath("testfiles/all_requests.json")

        if os.path.isfile(requests_file):
            # If requests_file exists
            with open(requests_file, 'a+') as outfile:
                # outfile.seek(-1, 2)
                outfile.seek(outfile.tell() - 1, 0)
                outfile.truncate()
                outfile.write(',')
                json.dump(request_info, outfile)
                outfile.write(']')
        else:
            # Creates requests file
            with open(requests_file, 'w') as outfile:
                array = []
                array.append(request_info)
                json.dump(array, outfile)


    def request_opener(self, requester):
        """ Reads UUID's of all the requests and compares those to the
        requester's UUID. Returns user_requests dict of user's requests.
        """
        user_requests = []

        # TODO: change here the correct path of the JSON-file
        # TODO: path of the JSON-file should be received as an argument
        # TODO: MAYBE THIS PATH SHOULD BE SET AS AN ENVIRONMENT VARIABLE
        requests_file = os.path.abspath("testfiles/all_requests.json")
        print(requests_file)

        # Check if all_requests -file exists
        if os.path.isfile(requests_file):
            # Opens requests_file for reading
            with open(requests_file, "r") as file:
                request_dict2 = json.load(file)

                # Comparison, if requester and a requestOwner match adds request
                # to dict
                for i in request_dict2:
                    # Comparison if owner id matches
                    if requester == i.get("requestOwner", ""):
                        user_requests.append(i)

            print("User's requests")
            print(user_requests)

        return user_requests


    def delete_request(self, request):
        """ Deletes request's folder and all the contents. Updates
        all_requests.json file.
        """
        # TODO: path of the JSON-file should be received as an argument
        # TODO: MAYBE THIS PATH SHOULD BE SET AS AN ENVIRONMENT VARIABLE
        requests_file = ""
        print(requests_file)
        requests_file = os.path.abspath("testfiles/all_requests.json")
        print(requests_file)
        all_requests_updated = []

        # TODO: ADD path to correct request. Should we get posix path, name of
        #  the request or request UUID? MAYBE THIS PATH SHOULD BE SET AS AN
        #  ENVIRONMENT VARIABLE
        dir_path = ' '
        dir_path = os.path.abspath("testfiles/")
        print(dir_path)
        dir_path += '/'
        dir_path += request
        dir_path += '.potku'
        print(dir_path)

        # Check if all_requests -file exists
        if os.path.isdir(dir_path):
            try:
                # Delete the request folder
                shutil.rmtree(dir_path)

                # Opens all_requests.json and updates requests
                with open(requests_file, "r+") as file:
                    request_list = json.load(file)

                    # Comparison, if requester and a requestOwner match adds
                    # request to dict
                    for i in request_list:
                        # Comparison if owner id matches
                        if request != i.get("requestName", ""):
                            all_requests_updated.append(i)

                with open(requests_file, "w") as file:
                    file.write(json.dumps(all_requests_updated))
                    return 'Request deleted'

            except OSError as e:
                print("Error: %s : %s" % (dir_path, e.strerror))
                return e

        return 'Nothing deleted'
