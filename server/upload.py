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
import upload


from pathlib import Path
from load_request import LoadRequest
from e_load_measurement import LoadMeasurement

from werkzeug.utils import secure_filename
from potku.modules.selection import Selector


def upload_measurement(file, request_name, sample_name, measurement_name):
    """ Upload handler for samples and measurements to the request.

    Args:
        file: file to upload
        request_name: Requests name
        sample_name: Sample name where file will be uploaded
        measurement_name: Name for the uploaded measurement
    """
    # open current request
    current_request = LoadRequest()

    # form the file path
    filepath = Path("testfiles", f"{request_name}.potku",
                    request_name).resolve()
    # TODO: next conversion is needed as open_request is expecting a string
    filepath_string = os.path.abspath(filepath)
    current_request.open_request(filepath_string)

    helper = LoadMeasurement(current_request.request.samples.samples,
                             current_request.request.directory)

    helper.add_sample(name=sample_name)
    helper.load_measurement(helper, measurement_name, sample_name)

    sample_name = sample_name.replace(" ", "_")
    sample = current_request.request.samples.add_sample(name=sample_name)
    helper.tab_id = "measurement"

    # Uploaded file is temporarily saved to disk, as measurement.py's
    # add_measurement_file method does not work with file objects.
    try:
        # TODO: upload_tmp folder should be somewhere on server, atm testfiles
        #       folder is used for this
        upload_tmp = os.path.abspath("testfiles/")
        # securing the filename
        filename = secure_filename(file.filename)
        file.save(os.path.join(upload_tmp, filename))
        tmp_file_path = os.path.join(upload_tmp, filename)
    except OSError as e:
        print("Error: %s : %s" % (upload_tmp, e.strerror))
        return e

    upload.add_new_measur_simu(current_request,
                               "measurement", tmp_file_path, sample,
                               load_data=True,
                               object_name=helper.name)

    # Deleting file from upload_tmp (ATM testfiles folder)
    if os.path.isdir(upload_tmp) and os.path.isfile(os.path.join(upload_tmp,
                                                                 filename)):
        try:
            os.remove(os.path.join(upload_tmp, filename))
            print('tmp file removed')
        except OSError as e:
            print("Error: %s : %s" % (upload_tmp, e.strerror))
            return e

    return measurement_name


# TODO: This one is def add_new_tab from potku.py file from original Potku,
#  but renamed and customised for ePotku use.
def add_new_measur_simu(self, tab_type, filepath: Path, sample, file_current=0,
                        file_count=1, load_data=False, object_name="",
                        import_evnt_or_binary=False):
    """Add new Sample and measurement to Request

    Args:
        tab_type: Either "measurement" or "simulation".
        filepath: A Path representing measurement or simulation file
        path, or data path when creating a new measurement.
        sample: The sample under which the measurement or simulation is put.
        file_current: An integer representing which number is currently
        being read. (for GUI)
        file_count: An integer representing how many files will be loaded.
        load_data: A boolean representing whether to load data or not. This
            is to save time when loading a request and we do not want to
            load every measurement.
        object_name: When creating a new Measurement, this is the name
            for it.
        import_evnt_or_binary: Whether evnt or lst data is being imported
            or not.
        progress: a ProgressReporter object
    """
    try:
        cur_progress = (100 / file_count) * file_current
    except ZeroDivisionError:
        cur_progress = 0

    # TODO: This next line is disabled, because measurement.py line 145 uses
    #  endswith, and 'PosixPath' object has no attribute 'endswith'.
    # TODO: some versions of Potku use 'suffix' on measurement.py. If it
    #  crashes, try uncommenting/commenting next line.
    filepath = Path(filepath)
    rest = (100 - cur_progress) * 0.01

    if tab_type == "measurement":
        measurement = \
            self.request.samples.measurements.add_measurement_file(
                sample, filepath, self.tab_id, object_name,
                import_evnt_or_binary=import_evnt_or_binary,
                selector_cls=Selector)

        return measurement

    # TODO: This here is a placeholder for uploading new simulations.
    #  Uploading simulations can be done just like measurements
    if tab_type == "simulation":
        pass


if __name__ == '__main__':
    upload_measurement()
