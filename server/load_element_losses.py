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
import potku.modules.cut_file as cut_file
from potku.modules.element_losses import ElementLosses
from potku.modules.element import Element
from sort_element import sortt


# Functions for loading element losses from cut files for load_request
def get_losses(measurement, reference_cut_file, checked_cuts, partition_count):
    """Elemental losses data in selected measurement from selected cut files

    Args:
        measurement: A measurement class object
        reference_cut_file: the reference cut file
        checked_cuts: selected cut files for calculating losses
        partition_count: amount of partitions

    Return:
        Elemental losses data: data, color, label
    """

    losses = ElementLosses(measurement.directory_cuts,
                           measurement.directory_composition_changes,
                           reference_cut_file,
                           checked_cuts,
                           partition_count)

    split_counts = losses.count_element_cuts()
    scale_mode = 1

    # From elementlosseswidget
    # Check for RBS selections.
    rbs_list = {}
    for cut in checked_cuts:
        filename = os.path.basename(cut)
        split = filename.split(".")
        if cut_file.is_rbs(cut):
            # This should work for regular cut and split.
            key = "{0}.{1}.{2}.{3}".format(split[1], split[2],
                                           split[3], split[4])
            rbs_list[key] = cut_file.get_scatter_element(cut)

    reference_cut_key = None
    if reference_cut_file:
        end_ref_cut = os.path.split(reference_cut_file)[-1]
        without_mes_name = end_ref_cut.split(".", 1)
        without_cut_suffix = without_mes_name[-1].rsplit(".", 1)
        reference_cut_key = without_cut_suffix[0]

    keys = [item[0] for item in sorted(split_counts.items(),
                                       key=lambda x: sortt(x[0]))]

    all_losses = {}
    for i, key in enumerate(keys, 0):
        c_file = key.split('.')
        # TODO provide elements as parameter rather than initializing
        #      them here
        element_object = Element.from_string(c_file[0].strip())
        element = element_object.symbol
        isotope = element_object.isotope
        if reference_cut_key and key == reference_cut_key:
            continue
        # Check RBS selection
        rbs_string = ""
        if len(c_file) == 3:
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

        # Get color for selection
        if isotope is None:
            isotope = ""
        if rbs_string == "*":
            color_string = "{0}{1}{2}".format("RBS_" + isotope, element,
                                              c_file[2])
        else:
            color_string = "{0}{1}{2}".format(isotope, element, c_file[2])
        if color_string not in measurement.selector.get_colors().keys():
            color = "red"
        else:
            color = measurement.selector.get_colors()[color_string]

        # Set label text
        if len(c_file) == 3:
            label = str(isotope) + element + rbs_string
        else:
            label =str(isotope) + element + rbs_string \
                    + c_file[3]
        # Modify data if scaled to 100.
        data = split_counts[key]
        if scale_mode == 2:
            n = None
            for val in data:
                if val == 0:
                    continue
                else:
                    n = val
                    break
            modifier = 100 / n  # self.split[key][0]
            data = [i * modifier for i in data]

        all_losses[i] = {'data': data,
                         'color': color,
                         'label': label}

    return all_losses
