# coding=utf-8
"""
Created on 4.5.2018
Updated on 24.5.2019

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

import json
import time

import dialogs.dialog_functions as df
import widgets.binding as bnd
import modules.general_functions as gf

from pathlib import Path

from modules.run import Run
from modules.simulation import Simulation

from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import uic

from widgets.detector_settings import DetectorSettingsWidget
from widgets.measurement.settings import MeasurementSettingsWidget


class SimulationSettingsDialog(QtWidgets.QDialog):
    """
    Dialog class for handling the simulation parameter input.
    """
    use_request_settings = bnd.bind("defaultSettingsCheckBox")

    def __init__(self, tab, simulation: Simulation, icon_manager):
        """
        Initializes the dialog.

        Args:
            tab: A SimulationTabWidget.
            simulation: A Simulation object whose parameters are handled.
            icon_manager: An icon manager.
        """
        super().__init__()
        uic.loadUi(Path("ui_files", "ui_specific_settings.ui"), self)

        self.tab = tab
        self.simulation = simulation
        self.icon_manager = icon_manager

        self.setWindowTitle("Simulation Settings")
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        screen_geometry = QtWidgets.QDesktopWidget.availableGeometry(
            QtWidgets.QApplication.desktop())
        self.resize(self.geometry().width() * 1.2,
                    screen_geometry.size().height() * 0.8)
        self.defaultSettingsCheckBox.stateChanged.connect(
            self.__change_used_settings)
        self.OKButton.clicked.connect(self.__save_settings_and_close)
        self.applyButton.clicked.connect(self.__update_parameters)
        self.cancelButton.clicked.connect(self.close)

        # Add measurement settings view to the settings view
        self.measurement_settings_widget = MeasurementSettingsWidget(
            self.simulation)
        self.tabs.addTab(self.measurement_settings_widget, "Measurement")

        self.measurement_settings_widget.beam_selection_ok.connect(
            lambda b: self.OKButton.setEnabled(b)
        )

        # Add detector settings view to the settings view
        detector_object = self.simulation.detector
        if detector_object is None:
            detector_object = self.simulation.request.default_detector

        self.detector_settings_widget = DetectorSettingsWidget(
            detector_object, self.simulation.request, self.icon_manager)

        # 2 is calibration tab that is not needed
        calib_tab_widget = self.detector_settings_widget.tabs.widget(2)
        self.detector_settings_widget.tabs.removeTab(2)
        calib_tab_widget.deleteLater()

        self.tabs.addTab(self.detector_settings_widget, "Detector")

        self.use_request_settings = self.simulation.use_request_settings

        self.measurement_settings_widget.nameLineEdit.setText(
            self.simulation.measurement_setting_file_name)
        self.measurement_settings_widget.descriptionPlainTextEdit \
            .setPlainText(
                self.simulation.measurement_setting_file_description)
        self.measurement_settings_widget.dateLabel.setText(time.strftime(
            "%c %z %Z", time.localtime(self.simulation.modification_time)))

        self.tabs.currentChanged.connect(lambda: df.check_for_red(self))

        self.exec()

    def __change_used_settings(self):
        """Set specific settings enabled or disabled based on the "Use
        request settings" checkbox.
        """
        check_box = self.sender()
        if check_box.isChecked():
            self.tabs.setEnabled(False)
        else:
            self.tabs.setEnabled(True)

    def __update_parameters(self):
        """
         Update Simulation's Run, Detector and Target objects. If simulation
         specific parameters are in use, save them into a file.
        """
        if self.measurement_settings_widget.isotopeComboBox.currentIndex()\
                == -1:
            QtWidgets.QMessageBox.critical(
                self, "Warning",
                "No isotope selected.\n\n"
                "Please select an isotope for the beam element.",
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)
            return False

        if not self.simulation.measurement_setting_file_name:
            self.simulation.measurement_setting_file_name = \
                self.simulation.name

        if not self.tabs.currentWidget().fields_are_valid:
            QtWidgets.QMessageBox.critical(
                self, "Warning",
                "Some of the setting values have not been set.\n"
                "Please input values in fields indicated in red.",
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)
            return False

        if (self.use_request_settings != self.simulation.use_request_settings) \
                or (not self.use_request_settings and self.values_changed()):
            # User has switched from simulation settings to request settings,
            # or vice versa. Confirm if the user wants to delete old simulations
            # OR
            # If values that require rerunning simulations, prompt user
            # delete previous or currently running simulation results (if those
            # exists)
            if not df.delete_element_simulations(
                    self, self.simulation, tab=self.tab,
                    msg="simulation settings"):
                return False

        try:
            # Update simulation settings
            self.simulation.use_request_settings = \
                self.use_request_settings
            measurement_settings_file_path = Path(
                self.simulation.directory,
                f"{self.simulation.measurement_setting_file_name}.measurement")
            target_file_path = Path(self.simulation.directory,
                                    f"{self.simulation.target.name}.target")
            det_folder_path = Path(self.simulation.directory, "Detector")

            if self.simulation.run is None:
                # Create a default Run for simulation
                self.simulation.run = Run()
            if self.simulation.detector is None:
                df.update_detector_settings(self.simulation,
                                            det_folder_path,
                                            measurement_settings_file_path)

            # Set Detector object to settings widget
            self.detector_settings_widget.obj = self.simulation.detector

            # Update settings
            self.measurement_settings_widget.update_settings()
            self.detector_settings_widget.update_settings()
            self.simulation.detector.path = \
                Path(det_folder_path,
                     f"{self.simulation.detector.name}.detector")

            # Save measurement settings parameters.
            new_measurement_settings_file_path = Path(
                self.simulation.directory,
                self.simulation.measurement_setting_file_name +
                ".measurement")
            general_obj = {
                "name": self.simulation.measurement_setting_file_name,
                "description":
                    self.simulation.measurement_setting_file_description,
                "modification_time":
                    time.strftime("%c %z %Z", time.localtime(
                        time.time())),
                "modification_time_unix": time.time(),
                "use_request_settings": self.simulation.use_request_settings
            }

            if new_measurement_settings_file_path.exists():
                with open(new_measurement_settings_file_path) as mesu:
                    obj = json.load(mesu)
                obj["general"] = general_obj
            else:
                obj = {
                    "general": general_obj
                }

            # Delete possible extra .measurement files
            gf.remove_files(self.simulation.directory,
                            exts={".measurement"})

            # Write measurement settings to file
            with open(new_measurement_settings_file_path, "w") as file:
                json.dump(obj, file, indent=4)

            self.simulation.to_file(self.simulation.path)

            # Save Run object to file
            self.simulation.run.to_file(new_measurement_settings_file_path)
            # Save Detector object to file
            self.simulation.detector.to_file(self.simulation.detector.path,
                                             new_measurement_settings_file_path)

            # Save Target object to file
            self.simulation.target.to_file(target_file_path,
                                           new_measurement_settings_file_path)
            return True

        except TypeError:
            QtWidgets.QMessageBox.question(
                self, "Warning",
                "Some of the setting values have not been set.\n"
                "Please input setting values to save them.",
                QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.Ok)

        return False

    def __save_settings_and_close(self):
        """Saves settings and closes the dialog if __update_parameters returns
        True.
        """
        if self.__update_parameters():
            self.close()

    def values_changed(self):
        """
        Check if measurement or detector settings have
        changed.

        Return:

            True or False.
        """
        if self.measurement_settings_widget.are_values_changed():
            return True
        if self.detector_settings_widget.values_changed():
            return True
        return False
