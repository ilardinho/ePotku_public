# coding=utf-8
"""
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

import upload
import request_helper

from load_request import LoadRequest
from e_new_request import NewRequest

"""Initializing objects after getting requests from
   app.py
"""


def get_request_measurement():
    request = LoadRequest()
    data = request.get_measurement_data()
    return data


def load_testrequest(compression_x, compression_y, data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_data(compression_x, compression_y)


def load_selections(data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_all_selections()


def load_depth(depth_units, cuts, error, data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_default_depth(depth_units, cuts, error)


def load_element_losses(cuts, ref_cut, splits, data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_element_losses(cuts, ref_cut, splits)


def save_splits(selected_cuts, reference_cut, partition_count, data=None):
    request = LoadRequest()
    request.open_default()
    return request.save_splits(selected_cuts, reference_cut, partition_count)


def save_cuts(cuts, data=None):
    request = LoadRequest()
    request.open_default()
    return request.save_cuts(cuts)


def load_energy_spectra(cuts, bin_width, data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_energy_spectra(cuts, bin_width)


def get_cuts(data=None):
    request = LoadRequest()
    request.open_default()
    return request.get_cuts()


def new_request(requestName, requestOwner, published):
    new_requ = NewRequest()
    print('Creating a new request: ' + requestName)
    new_requ.create_request(requestName, requestOwner, published)
    success = 'Great success, request created!'
    return success


def get_requests(requestOwnerUUID):
    opener = NewRequest()
    requests = opener.request_opener(requestOwnerUUID)
    return requests


def delete_request(request_name):
    delete_requ= NewRequest()
    deleted_requ = delete_requ.delete_request(request_name)
    return deleted_requ


def request_content(request_name):
    request_files = request_helper.request_content(request_name)
    return request_files


def measurement_upload(file, request_name, sample_name, measurement_name):
    uploaded = upload.upload_measurement(file, request_name, sample_name,
                                         measurement_name)
    return uploaded

