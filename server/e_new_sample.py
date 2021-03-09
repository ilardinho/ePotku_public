# coding=utf-8
"""
This module is copied from Potku. Modifications:
- Removed Qt-related code
- (Does not do anything yet)

Created on 26.2.2018
Updated on 12.6.2018

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
__author__ = "Severi Jääskeläinen \n Samuel Kaiponen \n Heta Rekilä \n " \
             "Sinikka Siironen"
__version__ = "2.0"

from pathlib import Path


# TODO: Use Sample directly instead. NewSample was originally used for
#       dialogs. This does not seem to do anything.
class NewSample:
    """Class for creating a new sample.
    """

    def __init__(self, samples):
        """Init new sample.

        Args:
            samples: List of samples.
        """

        self.name = ""
        self.description = ""
        self.samples = samples

    # TODO: This should not probably be used. There are no dialogs here.
    # TODO: here we need name attribute !
    def create_sample(self, new_sample_name):
        """Read sample name from view and if it is accepted, close dialog.
        """
        self.name = new_sample_name.replace(" ", "_")
        for sample in self.samples:
            if sample.name == self.name:
                # TODO: Raise an exception instead
                print(f"Sample with name '{self.name}' already exists. "
                      f"Choose another name.")
                break
