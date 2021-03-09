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

import os
from potku.modules.energy_spectrum import EnergySpectrum
# TODO: The next import is cumbersome to use, but
#       "from potku.modules import cut_file" does not seem to work
import potku.modules.cut_file
from potku.modules.element import Element
from sort_element import sortt


# Customized code from EnergySpectrumWidget from potku
def get_spectra(measurement, use_cuts, bin_width=0.025, save_file_int=0):
    """Energy spectra data in selected measurement and from selected cut files

    Args:
        measurement: A measurement class object
        use_cuts: cut files to be used
        bin_width: bin width
        save_file_int: for something?

    Return:
        Energy spectra data: x, y, color, label
    """

    save_file = "widget_energy_spectrum.save"
    rbs_list = {}

    measurement.generate_tof_in()
    energy_spectrum = EnergySpectrum(measurement, use_cuts, bin_width)
    es_data = energy_spectrum.calculate_spectrum()
    selection_colors = measurement.selector.get_colors()

    for cut in use_cuts:
        filename = os.path.basename(cut)
        split = filename.split(".")
        if potku.modules.cut_file.is_rbs(cut):
            key = "{0}.{1}.{2}.{3}".format(split[1], split[2],
                                           split[3], split[4])
            rbs_list[key] = potku.modules.cut_file.get_scatter_element(cut)

    # From this point on code is from matplotlibenergyspectrumwidget
    # self, es_data, rbs_list, spectrum_type, changed yms.
    # siirtyminen suoraan em. on_draw metodiin

    all_lines = {}

    element_counts = {}
    keys = [item[0] for item in sorted(es_data.items(),
                                       key=lambda x: sortt(
                                            x[0]))]
    indx = 0
    for key in keys:
        cut_file = key.split('.')
        cut = es_data[key]
        element_object = Element.from_string(cut_file[0])
        element = element_object.symbol
        isotope = element_object.isotope
        #if key in ignore_elements:
        #    continue

        # Check RBS selection
        rbs_string = ""
        if len(cut_file) == 3:
            if key + ".cut" in rbs_list.keys():
                element_object = rbs_list[key + ".cut"]
                element = element_object.symbol
                isotope = element_object.isotope
                rbs_string = "*"
        else:
            if key in rbs_list.keys():
                element_object = rbs_list[key]
                element = element_object.symbol
                isotope = element_object.isotope
                rbs_string = "*"

        x, y = get_axis_values(cut)

        if isotope is None:
            isotope = ""

        # Get color for selection
        dirtyinteger = 0
        if rbs_string == "*":
            color_string = "{0}{1}{2}{3}".format("RBS_", isotope,
                                                 element, dirtyinteger)
        else:
            color_string = "{0}{1}{2}".format(isotope, element,
                                              dirtyinteger)

        while color_string in element_counts:
            dirtyinteger += 1
            if rbs_string == "*":
                color_string = "{0}{1}{2}{3}".format("RBS_", isotope,
                                                     element,
                                                     dirtyinteger)
            else:
                color_string = "{0}{1}{2}".format(isotope, element,
                                                  dirtyinteger)

        element_counts[color_string] = 1
        if color_string not in selection_colors:
            color = "red"
        else:
            color = selection_colors[color_string]

        if len(cut_file) == 3:
            label =  str(isotope) + element + rbs_string
        else:
            label = str(isotope) + element \
                    + rbs_string + cut_file[2]

        all_lines[indx]= {'x': x,
                          'y': y,
                          'color': color,
                          'label': label}

        indx+=1

    return all_lines


def get_axis_values(data):
    """Returns the x and y axis values from given data."""
    return (
        tuple(float(pair[0]) for pair in data),
        tuple(float(pair[1]) for pair in data)
    )


def fix_minimum(lst, minimum):
    if lst and lst[0] < minimum:
        return lst[0], True
    return minimum, False
