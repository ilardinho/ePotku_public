# coding=utf-8
"""
Created on 15.3.2018
Updated on 30.5.2018

Potku is a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2018 Severi Jääskeläinen, Samuel Kaiponen, Heta Rekilä and
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

import warnings

from PyQt5 import QtWidgets

from widgets.matplotlib.base import MatplotlibWidget

warnings.warn("widgets.matplotlib.simulation.energy_spectrum module is "
              "deprecated.", DeprecationWarning, stacklevel=2)


class MatplotlibSimulationEnergySpectrumWidget(MatplotlibWidget):
    """Energy spectrum widget
    """
    def __init__(self, parent, data):
        """Inits Energy Spectrum widget.

        Args:
            parent: EnergySpectrumWidget class object.
            data: Energy spectrum data.
        """
        super().__init__(parent)
        super().fork_toolbar_buttons()
        # self.draw_legend = legend
        self.energy_spectrum_data = data
        # self.__rbs_list = rbs_list
        self.__icon_manager = parent.icon_manager
        # self.__selection_colors = parent.measurement.selector.get_colors()
        
        self.__initiated_box = False
        self.__ignore_elements = []
        self.__log_scale = False
        
        self.canvas.manager.set_title("Energy Spectrum")
        self.axes.fmt_xdata = lambda x: "{0:1.2f}".format(x)
        self.axes.fmt_ydata = lambda y: "{0:1.0f}".format(y)
        
        self.mpl_toolbar.addSeparator()
        self.__button_toggle_log = QtWidgets.QToolButton(self)
        self.__button_toggle_log.clicked.connect(self.__toggle_log_scale)
        self.__button_toggle_log.setCheckable(True)
        self.__button_toggle_log.setToolTip("Toggle logarithmic Y axis scaling")
        self.__icon_manager.set_icon(self.__button_toggle_log,
                                     "monitoring_section.svg")
        self.mpl_toolbar.addWidget(self.__button_toggle_log)
        
        self.__button_ignores = QtWidgets.QToolButton(self)
        self.__button_ignores.setToolTip("Select elements which are "
                                         "included in the graph.")
        self.__icon_manager.set_icon(self.__button_ignores, "gear.svg")
        self.mpl_toolbar.addWidget(self.__button_ignores)
        
        self.on_draw()

    def on_draw(self):
        """Draw method for matplotlib.
        """
        # Values for zoom
        x_min, x_max = self.axes.get_xlim()
        y_min, y_max = self.axes.get_ylim()

        self.axes.clear()  # Clear old stuff

        self.axes.set_ylabel("Yield (counts)")
        self.axes.set_xlabel("Energy (MeV)")

        x = tuple(float(pair[0]) for pair in self.energy_spectrum_data)
        y = tuple(float(pair[1]) for pair in self.energy_spectrum_data)

        self.axes.plot(x, y)

        if 0.09 < x_max < 1.01:  # This works...
            x_max = self.axes.get_xlim()[1]
        if 0.09 < y_max < 1.01:
            y_max = self.axes.get_ylim()[1]

        # Set limits accordingly
        self.axes.set_ylim([y_min, y_max])
        self.axes.set_xlim([x_min, x_max])

        # Remove axis ticks
        self.remove_axes_ticks()

        # Draw magic
        self.canvas.draw()

    def __toggle_log_scale(self):
        """Toggle log scaling for Y axis in depth profile graph.
        """
        self.__log_scale = self.__button_toggle_log.isChecked()
        self.on_draw()
