# coding=utf-8
"""
Created on x.x.2020
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

import numpy as np

from pathlib import Path
from flask import jsonify

from potku.modules.global_settings import GlobalSettings
from potku.modules.request import Request
import potku.modules.math_functions as mf
from potku.modules.element_losses import ElementLosses
from potku.modules.element import Element
import load_element_losses
import load_energy_spectrum
import load_depth_profile

from potku.modules.selection import Selector


# Largely code from potku.py, modified to work here
class LoadRequest:
    """ Load_request class for loading a request."""

    def __init__(self):
        self.request = None
        # TODO: tab_widgets and tab_id are for Potku UI. Remove them.
        self.tab_widgets = {}
        self.tab_id = 0
        # TODO: GlobalSettings should be a per-user thing. This is a
        #       temporary fix to prevent overwriting Potku files.
        self.settings = GlobalSettings(Path.cwd() / "testfiles")
        self.l_measurements = {}
        self.current_measurement = None
        self.data = None

        self.__x_data = None
        self.__y_data = None

    def load_request_samples(self):
        """"Load sample files in the request. From potku.py"""
        sample_paths_in_request = self.request.get_samples_files()
        if sample_paths_in_request:
            for sample_path in sample_paths_in_request:
                sample = self.request.samples.add_sample(sample_path)
            self.request.increase_running_int_by_1()

    def get_request_samples_and_measurements(self):
        """Returns a dictionary of samples and measurements in the request.
        """
        s_and_m = {}
        a = self.request.samples.get_samples_and_measurements()
        for sample in a:
            sample_serial = sample.serial_number
            sample_name = sample.name
            measurements = []
            for pathto in a[sample]:
                p = os.path.split(pathto)[1].split('.')[0]
                measurements.append(p)
            s_and_m[sample_serial] = {'measurements' : measurements,
                                      'name' : sample_name}

        return jsonify(s_and_m)

    def get_user_requests(self, user_id):
        """Get a list of requests owned by a user"""
        direc = os.path.abspath("testfiles/"+user_id)
        requests = os.listdir(direc)
        requestnames = [i.split(".")[0] for i in requests]

        return jsonify(requestnames)

    # TODO: this should be more like load_request_measurements
    def open_measurement(self, sample_id, measurement):
        """Opens a measurement from a specified sample."""
        ssample = None
        for sample in self.request.samples.samples:
            if sample.serial_number == int(sample_id):
                ssample = sample

        measurements = []
        samples_with_measurements = \
            self.request.samples.get_samples_and_measurements()

        # sample_id:n mukaan voi laittaa suoraan indeksin tähän
        for sample, measurements in samples_with_measurements.items():
            for measurement_file in measurements:
                if sample.serial_number == int(sample_id) and \
                        os.path.split(measurement_file)[1].split('.')[0] == measurement:
                    self.current_measurement = \
                    self.request.samples.measurements.add_measurement_file(
                        sample, measurement_file, 0, "", 0, selector_cls=Selector)

        # this bit is unnecessary
        self.__x_data = [x[0] for x in self.current_measurement.data]
        self.__y_data = [x[1] for x in self.current_measurement.data]

        # TODO: are these tests?
        self.current_measurement.update_folders_and_selector()
        self.current_measurement._Measurement__check_for_old_selection()

    def load_request_measurements(self, data=None, measurements=None):
        """Load measurement files in the request.

        Args:
            measurements: A list representing loadable measurements when
            importing measurements to the request.
        """
        if measurements is None:
            measurements = []
        if measurements:
            samples_with_measurements = measurements
            load_data = True
        else:
            # a dict with the sample as a key, and measurements'
            # info file paths in the value as a list
            samples_with_measurements = \
                self.request.samples.get_samples_and_measurements()
            load_data = False

        # Adding measurements to a dict
        for i, (sample, measurements) in \
                enumerate(samples_with_measurements.items(), 0):
            for measurement_file in measurements:
                self.l_measurements[i] = ("measurement",
                                          measurement_file, sample)
                #dirtyinteger += 1

                self.current_measurement = \
                    self.request.samples.measurements.add_measurement_file(
                        sample, measurement_file, 0, "", 0,
                        selector_cls=Selector)
                if self.current_measurement is not None:
                    self.current_measurement.load_data()

        try:
            self.data = self.current_measurement.data
            self.current_measurement.update_folders_and_selector()
            self.current_measurement._Measurement__check_for_old_selection()
        except AttributeError:
            print('This request did not have any measurements yet')


    # TODO: is this used for anything?
    def get_measurement_data(self, measurements=None):
        print('getting measurement data')
        # file = os.path.abspath("testfiles/gradu_testi.potku/gradu_testi.request")
        file = os.path.abspath("testfiles/Potku_v2.potku/Potku_v2.request")
        try:
            request = Request.from_file(file, self.settings, self.tab_widgets)
        except Exception as e:
            print('something went wrong while getting measurement data')
            print(e)
            return
        self.request = request

        folder = os.path.split(file)[0]
        self.settings.set_request_directory_last_open(folder)

        self.load_request_samples()

        if measurements is None:
            measurements = []
        if measurements:
            samples_with_measurements = measurements
            load_data = True
        else:
            samples_with_measurements = \
                self.request.samples.get_samples_and_measurements()
            load_data = False

        dirtyinteger = 0
        for sample, measurements in samples_with_measurements.items():
            for measurement_file in measurements:
                self.l_measurements[dirtyinteger] = ("measurement", measurement_file, sample)
                dirtyinteger += 1

                self.current_measurement = \
                    self.request.samples.measurements.add_measurement_file(
                        sample, measurement_file, 0, "", 0, selector_cls=Selector)

                if self.current_measurement is not None:
                    measurement_data = self.current_measurement.load_data()
                    return measurement_data

    def __open_request(self, file, data=None):
        """Opens a request"""
        try:
            request = Request.from_file(file, self.settings, self.tab_widgets)
        except Exception as e:
            print('error while opening request')
            print(e)
            return
        self.request = request

        folder = os.path.split(file)[0]
        self.settings.set_request_directory_last_open(folder)

        self.load_request_samples()
        self.load_request_measurements(data)

    def get_data(self, compression_x, compression_y):
        """Creates a 2D histogram from the measurement using
           compression values. Does binning for the raw .asc
           data.

        Args:
            compression_x: value for compression of x-axis
            compression_y: value for compression of y-axis

        Return:
            The histogram as a bytes object
        """

        # TODO: finish constraining compression values if they exceed max/min
        constrained = False
        if compression_x < 4:
            compression_x = 4
            constrained = True
        if compression_y < 4:
            compression_y = 4
            constrained = True
        if compression_y > 20:
            compression_y = 20
            constrained = True
        if compression_x > 20:
            compression_x = 20
            constrained = True

        # adapted to work for ndarray from math_functions.calculate_bin_counts()
        x0 = int(np.amin(self.data[0]))
        xn = int(np.amax(self.data[0]))
        y0 = int(np.amin(self.data[1]))
        yn = int(np.amax(self.data[1]))
        bin_count_x = int(abs(x0 - xn) / compression_x)
        bin_count_y = int(abs(y0 - yn) / compression_y)
        bin_counts = (bin_count_x, bin_count_y)

        counts, xedges, yedges = np.histogram2d(self.data[0], self.data[1], bins=bin_counts)
        histogram = counts.astype(np.uint16)

        # replace unused part of the histogram to transmit metadata
        histogram[0,0] = xedges[0]
        histogram[0,1] = yedges[0]
        histogram[0,2] = bin_count_x
        histogram[0,3] = bin_count_y

        return histogram.tobytes()

    def get_all_selections(self):
        """Initializes values from existing selections

        Return:
            All selections from the existing folder structure
        """

        all_sel = {}
        key = 0
        for sel in self.current_measurement.selector.selections:
            p = sel.get_points()
            scatter = None

            # if RBS, element scatter is an Element object, otherwise str
            if sel.type == "RBS":
                scatter = sel.element_scatter.symbol
            else:
                scatter = sel.element_scatter

            # TODO: if RBS, isotope should be that of the scattered element, not the beam ion?
            all_sel[key] = {'element_type': sel.type,
                            'element': sel.element.symbol,
                            'isotope': sel.element.isotope,
                            'weight_factor': sel.weight_factor,
                            'scatter': scatter,
                            'color': sel.default_color,
                            'points': p}
            key += 1

        return all_sel

    def get_cuts(self):
        """Reads existing cut files

        Return:
            The names of all existing cut files from the folder structure
        """
        all_cuts = []
        cuts_directory = self.current_measurement.directory_cuts

        for cut in os.listdir(cuts_directory):
            if os.path.isfile(os.path.join(cuts_directory, cut)):
                all_cuts.append(cut)

        return jsonify(all_cuts)

    def get_element_losses(self, selected_cuts, reference_cut, splits):
        """Calculates the element loss data based on cuts

        Args:
            Cut file names as a list

            Return:
            Elemental loss data as JSON
        """

        reference_cut_file = \
            Path(self.current_measurement.directory_cuts, reference_cut)

        checked_cuts = []
        cuts_directory = self.current_measurement.directory_cuts
        for ct in selected_cuts:
            checked_cuts.append(Path(cuts_directory, ct))

        values = load_element_losses.get_losses(
            self.current_measurement, reference_cut_file,
            checked_cuts, splits)

        return values

    # is this the best way to save splits?
    # Creates a new ElementLosses object and saves the splits from that
    def save_splits(self, selected_cuts, reference_cut_file, partition_count):
        """Function for saving splits from selected cut files

        Args:
            Cut file names as a list

        Return:
            Message indicating the save was successful
        """

        checked_cuts = []
        cuts_directory = self.current_measurement.directory_cuts

        for ct in selected_cuts:
            checked_cuts.append(Path(cuts_directory, ct))

        losses = ElementLosses(self.current_measurement.directory_cuts,
                               self.current_measurement.directory_composition_changes,
                               Path(self.current_measurement.directory_cuts, reference_cut_file),
                               checked_cuts,
                               partition_count)

        losses.count_element_cuts()
        losses.save_splits()

        return jsonify("splits saved successfully")

    # TODO: adapt to work for specific requests/measurements
    def save_cuts(self, data):
        """Save all user selections, then save cuts based on the selections"""

        selection_strings = []

        key = 0
        for sel in data:
            x = data[key]['x']
            y = data[key]['y']
            x = ','.join(str(i) for i in x)
            y = ','.join(str(j) for j in y)
            coord_string = ';'.join((x, y))

            if data[key]['elementType'] == 'ERD':
                save_string = "{0}    {1}    {2}    {3}    {4}    {5}    {6}\n".\
                    format(
                        data[key]['elementType'],
                        data[key]['element'],
                        data[key]['isotope'],
                        float(data[key]['weightFactor']),
                        "",
                        data[key]['color'],
                        coord_string
                    )
            # TODO: element and isotope for RBS are those of the ion beam,
            #       and can be found in the relevant .MEASUREMENT file
            else:
                save_string = "{0}    {1}    {2}    {3}    {4}    {5}    {6}\n".\
                    format(
                        data[key]['elementType'],
                        'Cl',
                        '35',
                        float(data[key]['weightFactor']),
                        data[key]['element'],
                        data[key]['color'],
                        coord_string
                    )

            selection_strings.append(save_string)
            key += 1

        # Load default selections
        selections_file = os.path.abspath(
            "testfiles/Potku_v2.potku/Sample_01-s1/Measurement_01-ToF-E_815/Data/ToF-E_815.selections"
        )
        with open(selections_file, "w") as fp:
            for sel in selection_strings:
                fp.write(sel)

        self.current_measurement._Measurement__check_for_old_selection()
        print('saving cuts...')
        self.current_measurement.save_cuts()

        return jsonify("adding selection successful")

    def get_energy_spectra(self, selected_cuts, bin_width):
        """Function for getting data for drawing an energy spectra

        Args:
            selected_cuts: cuts from which the energy spectra is calculated
            bin_width: histogram bin width

        Return:
            Energy spectra data as JSON
        """

        checked_cuts = []
        cuts_directory = self.current_measurement.directory_cuts
        for ct in selected_cuts:
            checked_cuts.append(Path(cuts_directory, ct))

        results = load_energy_spectrum.get_spectra(self.current_measurement,
                                                   checked_cuts,
                                                   bin_width=bin_width)

        return jsonify(results)

    def get_default_depth(self, x_units, cuts, error):
        """Calculates data needed for drawing a depth profile

        Args:
            measurement: current measurement
            x_units: X axis units to be used in the depth profiles
            cuts: the cuts from which depth profiles are generated
            error: error margin

        Return:
            Depth profile data as JSON
        """

        data = load_depth_profile.get_depth(
            self.current_measurement, x_units, cuts, error)

        return jsonify(data)

    def open_default(self):
        # This implementation is vulnerable to path traversal attacks.
        # Do not use it as a reference for handling user-supplied input.
        working_directory = os.getcwd()
        default_request = os.path.normpath(
            'testfiles/Potku_v2.potku/Potku_v2.request')
        file = os.path.realpath(
            os.path.join(working_directory, default_request))
        print('opening request...')
        self.__open_request(file)

    def open_request(self, user_id):
        """Opens a request based (just) on user id"""
        path = Path("testfiles", user_id+".potku", user_id+".request")
        print(path)
        file = os.path.abspath(path)
        self.__open_request(file)

    def open_request2(self, user_id, request_name):
        """Opens a request based on user id and request name

        Args:
            user_id: id of the user
            request_name: e.g. (request_name).potku
        """
        path = Path("testfiles", user_id, request_name + ".potku", request_name + ".request")
        print(path)
        file = os.path.abspath(path)
        self.__open_request(file)


# main for debugging
if __name__ == '__main__':
    l = LoadRequest()
    l.open_request("1337")
    l.get_all_selections()
    a = l.get_element_losses(["m1.1H.ERD.0.cut"])
    print(a)
