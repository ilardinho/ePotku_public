# coding=utf-8
"""
Created on 2.2.2020

Potku is a graphical user interface for analyzation and
visualization of measurement data collected from a ToF-ERD
telescope. For physics calculations Potku uses external
analyzation components.
Copyright (C) 2020 TODO

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
__author__ = "Juhani Sundell"
__version__ = "2.0"


import unittest
import tempfile
import os

import tests.mock_objects as mo
import tests.utils as utils

from pathlib import Path
from unittest.mock import patch
from unittest.mock import Mock
from modules.simulation import Simulation
from modules.enums import OptimizationType


class TestSimulation(unittest.TestCase):
    def test_slots(self):
        """Tests that __slots__ work correctly for Simulation."""

        with tempfile.TemporaryDirectory() as temp_dir:
            req = mo.get_request()
            sim = Simulation(os.path.join(temp_dir, "test.simu"), req)

            # Logging needs to be disabled, otherwise loggers retain file
            # handlers that prevent removing the temp_dir
            utils.disable_logging()

            self.assertRaises(AttributeError, lambda: utils.slots_test(sim))

        # Just in case make sure that the temp_dir got deleted
        self.assertFalse(os.path.exists(temp_dir))

    @patch("modules.mcerd.MCERD")
    @patch("modules.mcerd.MCERD.run")
    @patch("modules.mcerd.MCERD.stop_process")
    def test_get_active_simulation(self, mock1, mock2, mock3):
        sim = mo.get_simulation()

        sim.add_element_simulation(mo.get_recoil_element())
        self.assertEqual(([], [], [], []), sim.get_active_simulations())

        elem_sim = sim.element_simulations[0]
        elem_sim.use_default_settings = False
        elem_sim.start(1, 1).subscribe(Mock())
        self.assertEqual(([elem_sim], [], [], []), sim.get_active_simulations())

        elem_sim._set_flags(False)
        self.assertEqual(
            ([], [elem_sim], [], []), sim.get_active_simulations())

        elem_sim.optimization_recoils = [mo.get_recoil_element()]
        elem_sim.start(1, 1, optimization_type=OptimizationType.RECOIL)
        self.assertEqual(
            ([], [elem_sim], [elem_sim], []), sim.get_active_simulations())

        elem_sim._set_flags(False)
        self.assertEqual(
            ([], [elem_sim], [], [elem_sim]), sim.get_active_simulations())

    def test_serialization(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            fp = Path(tmp_dir, "test.simu")
            sim = Simulation(tmp_dir, mo.get_request(), name="foo",
                             description="bar", save_on_creation=False,
                             run=mo.get_run(), detector=mo.get_detector(),
                             target=mo.get_target())

            sim.to_file(fp)

            sim2 = Simulation.from_file(mo.get_request(), fp,
                                        detector=mo.get_detector(),
                                        target=mo.get_target(),
                                        run=mo.get_run())

            utils.disable_logging()

            self.assertEqual(sim.name, sim2.name)
            self.assertEqual(sim.description, sim2.description)
            self.assertEqual(sim.measurement_setting_file_description,
                             sim2.measurement_setting_file_description)
            self.assertEqual(sim.measurement_setting_file_name,
                             sim2.measurement_setting_file_name)
            self.assertEqual(sim.use_request_settings,
                             sim2.use_request_settings)

            self.assertEqual(sim.detector.name, sim2.detector.name)
            self.assertEqual(sim.detector.angle_offset,
                             sim2.detector.angle_offset)

            self.assertEqual(sim.target.scattering_element,
                             sim2.target.scattering_element)
            self.assertEqual(sim.target.target_type, sim2.target.target_type)

            self.assertEqual(sim.run.fluence, sim2.run.fluence)

            self.assertNotEqual(sim.modification_time, sim2.modification_time)
