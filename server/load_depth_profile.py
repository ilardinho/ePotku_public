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

from pathlib import Path
from potku.modules import depth_files
from potku.modules.depth_files import DepthProfileHandler
from potku.modules.element import Element


def get_depth(measurement, depth_units, selected_cuts, error):
    """Creates depth files and reads data from the same depth files.
       Parts of code from Matplotlibdepthprofilewidget

    Args:
        measurement: current measurement
        depth_units: X axis units to be used in the depth profiles
        selected_cuts: the cuts from which depth profiles are generated
        error: error

    Return:
        Depth profile data
    """

    cuts = []
    elements = []

    # cut name format "m1.12C.ERD.0.cut" or "m1.35Cl.RBS_Mn.0.cut"
    for cut in selected_cuts:
        # absolute path to a specific cut file
        cuts.append(Path(measurement.directory_cuts, cut))
        if 'RBS_' in cut:
            # e.g. get just 'Mn'
            element = (cut.split(".")[2]).replace('RBS_', '')
            elements.append(element)
        else:
            elements.append(cut.split(".")[1]) # e.g. get '12C'

    formatted_elements = [Element.from_string(element) for element in elements]

    path_d = Path(measurement.directory_depth_profiles)
    depth_files.generate_depth_files(cuts, path_d, measurement)

    dpf = DepthProfileHandler()
    dpf.read_directory(path_d, formatted_elements, depth_units=depth_units)
    # absolute_profiles: a dict with element or 'total' as a key
    absolute_profiles = dpf.get_absolute_profiles()
    relative_profiles = dpf.get_relative_profiles()

    depth_profiles = {
        "absolute": {},
        "relative": {}
    }

    # calculate_ratios() needs an iterable object,
    # but it's unnecessary here so we pass an empty set
    percentage_conc, error_margins = \
        dpf.calculate_ratios(set(), systematic_error=error)
    absolute_concentrations = dpf.integrate_concentrations()

    profiles_to_use = absolute_profiles
    # From matplotlibdepthprofilewidget
    sorted_profile_names = \
        sorted(filter(lambda x: x != "total", profiles_to_use),
               key=lambda x: profiles_to_use[x].element)

    # TODO: remove copy-paste code
    # From matplotlibprofilewidget as well
    for i, profile_name in enumerate(sorted_profile_names, 0):
        if profile_name == "total":
            continue
        element = profiles_to_use[profile_name].element

        x_ax = profiles_to_use[profile_name].depths
        y_ax = profiles_to_use[profile_name].concentrations
        label = str(element)

        depth_profiles['absolute'][i] = {
            'x': x_ax,
            'y': y_ax,
            'label': label,
            'percentage_concentration': percentage_conc[profile_name],
            'error_margin': error_margins[profile_name],
            'absolute_concentration': absolute_concentrations[profile_name]
        }

    profiles_to_use = relative_profiles
    # From matplotlibdepthprofilewidget
    sorted_profile_names = \
        sorted(filter(lambda x: x != "total", profiles_to_use),
               key=lambda x: profiles_to_use[x].element)

    # From matplotlibprofilewidget as well
    for i, profile_name in enumerate(sorted_profile_names, 0):
        if profile_name == "total":
            continue
        element = profiles_to_use[profile_name].element

        x_ax = profiles_to_use[profile_name].depths
        y_ax = profiles_to_use[profile_name].concentrations
        label = str(element)

        depth_profiles['relative'][i] = {
            'x': x_ax,
            'y': y_ax,
            'label': label,
            'percentage_concentration': percentage_conc[profile_name],
            'error_margin': error_margins[profile_name],
            'absolute_concentration': absolute_concentrations[profile_name]
        }

    return depth_profiles
