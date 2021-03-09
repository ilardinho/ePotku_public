# coding=utf-8
"""
Created on 4.5.2020
Updated on 30.6.2020

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


import os


def request_content(request_name):
    # TODO: take path as argument, not just request name
    """Gets all the samples and measurements in the request.

    Args:
        request_name: Requests name
    """
    # TODO: path to folder containing all the requests should be set as an
    #  environment variable
    # TODO: think this through again :D maybe all the measurement &
    #  simulation info will be stored in a database
    # TODO: use pathlib.Path instead of os.path
    all_requests_folder = os.path.abspath("testfiles")
    request_folder = os.path.join(all_requests_folder, request_name)

    request_dict = {}
    if not os.path.isdir(request_folder):
        # TODO: raise an exception
        return request_dict

    # TODO: Write logs to a log file
    try:
        # read the entries
        with os.scandir(request_folder) as entries:
            for entry in entries:
                # TODO: Load simulations too
                if entry.is_dir(follow_symlinks=False) and \
                        entry.name.startswith("Sample"):
                    print(f"Folder: {entry.name}")
                    request_dict[entry.name] = ""

            try:
                # read the entries in sample and measurement folders
                for key in request_dict:
                    print(f"Reading {key}")
                    folder_helper = os.path.join(request_folder, key)

                    sub_entry_dict = {}
                    with os.scandir(folder_helper) as sub_entries:
                        for entry in sub_entries:
                            if entry.is_dir(follow_symlinks=False) and \
                                    entry.name.startswith("Measurement"):
                                print(f"Subfolder: {entry.name}")
                                sub_entry_dict[entry.name] = entry.name

                    request_dict[key] = sub_entry_dict

            except OSError as e:
                print(f"Error: {request_folder} : {e.strerror}")
                return e

    except OSError as e:
        print(f"Error: {request_folder} : {e.strerror}")
        return e

    return request_dict
