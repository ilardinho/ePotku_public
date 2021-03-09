# coding=utf-8
"""
This module is copied from Potku. Modifications:
- Copied sortt from matploblib

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
__author__ = "Jarkko Aalto \n Timo Konu \n Samuli Kärkkäinen " \
             "\n Samuli Rahkonen \n Miika Raunio \n Severi Jääskeläinen " \
             "\n Samuel Kaiponen \n Heta Rekilä \n Sinikka Siironen"
__version__ = "1.0"


from potku.modules.element import Element


# Copied from potku/widgets/maplotlib/measurement/element_losses.py
def sortt(key):
    cut_file = key.split('.')
    # TODO use RBS selection to sort (for example m1.35Cl.RBS_Mn.0.cut
    #  should be sorted by Mn, not Cl)
    # TODO provide elements as parameter instead of reinitializing them
    #  over and over again
    return Element.from_string(cut_file[0].strip())
