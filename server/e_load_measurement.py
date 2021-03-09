# coding=utf-8
"""
This module is copied from Potku. Some modifications made by ePotku 
project team are be present, f.e. removing all QT and unnecessary
stuff.

Created on 26.2.2018
Updated on 30.6.2020

Potku is a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2013-2018 Jarkko Aalto, Severi Jääskeläinen, Samuel Kaiponen,
Timo Konu, Samuli Kärkkäinen, Samuli Rahkonen, Miika Raunio, Heta Rekilä and
Sinikka Siironen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (file named 'LICENCE').
"""
__author__ = "Severi Jääskeläinen \n Samuel Kaiponen \n Heta Rekilä " \
             "\n Sinikka Siironen"
__version__ = "2.0"

from pathlib import Path

from e_new_sample import NewSample
# from potku.modules.measurement import Measurements


class LoadMeasurement:
    """Class for loading a measurement.
    """
    def __init__(self, samples, directory):
        """Inits a measurement loader.

        Args:
            samples: Samples of request.
            directory: Directory where to open the file browser.
        """
        self.name = ""
        self.sample = None
        self.directory = directory
        self.filename = ""
        self.samples = samples
        # TODO: this might be wrong place for add_measurement_file
        # self.add_measurement_file = Measurements
        self.tab_id = ""

    def add_sample(self, name=""):
        sample_sample = NewSample(self.samples)
        sample_sample.create_sample(name)
        print(sample_sample)

    def load_measurement(self, m_path, m_name, sample_name):
        self.path = m_path
        self.name = m_name.replace(" ", "_")
        self.sample = sample_name.replace(" ", "_")
        if not self.path:
            return
        if not self.name:
            return
        if not self.sample:
            return

        sample = self.__find_existing_sample()

        if sample:
            # Check if measurement on the same name already exists.
            for key in sample.measurements.measurements.keys():
                if sample.measurements.measurements[key].name == self.name:
                    # TODO: Raise an exception instead
                    print(f"Measurement with name '{self.name}' already exists. "
                          f"Choose another name.")
                    self.__close = False
                    break

    def __find_existing_sample(self):
        """
        Find existing sample that matches the sample name in self.

        Return:
            Sample object or None.
        """
        for sample in self.samples:
            if "Sample " + "%02d" % sample.serial_number + " " + sample.name \
                    == self.sample:
                return sample
        return None
