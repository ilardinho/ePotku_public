# coding=utf-8
"""
Created on 25.4.2018
Updated on 8.2.2020

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

__author__ = "Severi Jääskeläinen \n Samuel Kaiponen \n Heta Rekilä \n" \
             "Sinikka Siironen \n Juhani Sundell"
__version__ = "2.0"

import os
import platform
import shutil
import subprocess
import re
import multiprocessing
import rx

from . import general_functions as gf
from . import observing

from pathlib import Path
from rx import operators as ops
from rx.scheduler import ThreadPoolScheduler

from .layer import Layer
from .concurrency import CancellationToken


class MCERD:
    """
    An MCERD class that handles calling the mcerd binary and creating the
    files it needs.
    """
    __slots__ = "__settings", "__rec_filename", "__filename", \
                "recoil_file", "sim_dir", "result_file", "__target_file", \
                "__command_file", "__detector_file", "__foils_file", \
                "__presimulation_file", "__seed"

    def __init__(self, seed, settings, filename, optimize_fluence=False):
        """Create an MCERD object.

        Args:
            settings: All settings that MCERD needs in one dictionary.
        """
        self.__seed = seed
        self.__settings = settings

        rec_elem = self.__settings["recoil_element"]

        if optimize_fluence:
            self.__rec_filename = f"{rec_elem.prefix}-optfl"
        else:
            self.__rec_filename = rec_elem.get_full_name()

        self.__filename = filename

        self.sim_dir = Path(self.__settings["sim_dir"])

        suffix = self.__settings["simulation_type"].get_recoil_suffix()

        res_file = f"{self.__rec_filename}.{self.__seed}.erd"

        # The recoil file and erd file are later passed to get_espe.
        self.recoil_file = self.sim_dir / f"{self.__rec_filename}.{suffix}"
        self.result_file = self.sim_dir / res_file

        # These files will be deleted after the simulation
        self.__command_file = self.sim_dir / self.__rec_filename
        self.__target_file = self.sim_dir / f"{self.__filename}.erd_target"
        self.__detector_file = self.sim_dir / f"{self.__filename}.erd_detector"
        self.__foils_file = self.sim_dir / f"{self.__filename}.foils"
        self.__presimulation_file = self.sim_dir / f"{self.__filename}.pre"

    def get_command(self):
        """Returns the command that is used to start the MCERD process.
        """
        if platform.system() == "Windows":
            cmd = str(gf.get_bin_dir() / "mcerd.exe")
        else:
            cmd = "./mcerd"

        return cmd, str(self.__command_file)

    def run(self, print_to_console=False, cancellation_token=None,
            poll_interval=10, first_check=0.2, max_time=None, ct_check=0.2):
        """Starts the MCERD process.

        Args:
            print_to_console: whether MCERD output is also printed to console
            cancellation_token: token that is checked periodically to see if
                the simulation should be stopped.
            poll_interval: seconds between each check to see if the simulation
                process is still running.
            first_check: seconds until the first time mcerd is polled.
            max_time: maximum running time in seconds.
            ct_check: how often cancellation is checked in seconds.

        Return:
            observable stream where each item is a dictionary. All dictionaries
            contain the same keys.
        """
        # Create files necessary to run MCERD
        self.__create_mcerd_files()

        cmd = self.get_command()
        if cancellation_token is None:
            cancellation_token = CancellationToken()

        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            cwd=gf.get_bin_dir())

        errs = rx.from_iterable(iter(process.stderr.readline, b""))
        outs = rx.from_iterable(iter(process.stdout.readline, b""))

        is_running = rx.timer(first_check, poll_interval).pipe(
            # TODO change this to run at an increasing interval, i.e:
            #       - first check after 0.0 seconds,
            #       - second check after 0.2 seconds,
            #       - third after 1.0, ... etc.
            #   MCERD is likely to crash early (?) so it makes sense to
            #   run the check more frequently at the beginning.
            ops.map(lambda _: {
                "is_running": MCERD.is_running(process)
            })
        )

        ct_check = rx.timer(0, ct_check).pipe(
            # Filter out all the positive values. We do not want to push
            # updates every time cancellation is checked as there is nothing
            # to update.
            ops.map(
                lambda _: MCERD._stop_if_cancelled(
                    process, cancellation_token)),
            ops.filter(lambda x: not x),
            ops.map(lambda _: {
                "is_running": False,
                "msg": "Simulation was stopped"
            })
        )

        if max_time is not None:
            timeout = rx.timer(max_time).pipe(
                ops.do_action(
                    # Request cancellation so all simulation processes that
                    # share the same cancellation_token are also stopped.
                    # TODO not working as intended if simulation is short enough
                    #   to stop before max_time has elapsed. Maybe let caller
                    #   implement its own timeout check when multiple processes
                    #   are being run.
                    on_next=lambda _:
                    cancellation_token.request_cancellation()),
                ops.map(lambda _: {
                    "is_running": bool(MCERD.stop_process(process)),
                    "msg": "Simulation timed out"
                })
            )
        else:
            timeout = rx.empty()

        thread_count = multiprocessing.cpu_count()
        pool_scheduler = ThreadPoolScheduler(thread_count)

        merged = rx.merge(errs, outs).pipe(
            ops.subscribe_on(pool_scheduler),
            MCERD.get_pipeline(self.__seed, self.__rec_filename),
            ops.combine_latest(rx.merge(
                is_running, ct_check, timeout
            )),
            ops.map(lambda x: {
                **x[0], **x[1],
                "is_running": x[0]["is_running"] and x[1]["is_running"]
            }),
            ops.take_while(lambda x: x["is_running"], inclusive=True),
        )

        if print_to_console:
            merged = merged.pipe(
                observing.get_printer(
                    f"simulation process with seed {self.__seed}.")
            )

        # on_completed does not get called if the take_while condition is
        # inclusive so this is a quick fix to get the files deleted.
        # TODO surely there is a way to get the on_completed called?
        def del_if_not_running(x):
            if not x["is_running"]:
                self.delete_unneeded_files()

        return merged.pipe(
            ops.do_action(
                on_next=del_if_not_running,
                on_error=lambda _: self.delete_unneeded_files(),
                on_completed=self.delete_unneeded_files)
        )

    @staticmethod
    def _stop_if_cancelled(process: subprocess.Popen,
                           cancellation_token: CancellationToken):
        """Stops the process if cancellation has been requested. Returns
        True if cancellation has not been requested and simulation is still
        running, False otherwise.
        """
        if cancellation_token.is_cancellation_requested():
            MCERD.stop_process(process)
            return False
        return True

    @staticmethod
    def is_running(process: subprocess.Popen) -> bool:
        """Checks if the given process is running. Raises SubprocessError if
        the process returns an error code.
        """
        res = process.poll()
        if res == 0:
            return False
        if res is None:
            return True
        raise subprocess.SubprocessError(
            f"MCERD stopped with an error code {res}.")

    @staticmethod
    def get_pipeline(seed: int, name: str) -> rx.pipe:
        """Returns an rx pipeline that parses the raw output from MCERD
        into dictionaries.

        Each dictionary contains the same keys. If certain value cannot be
        parsed from the output (i.e. the raw line does not contain it),
        either the value from the previous dictionary is carried over or a
        default value is used.

        Args:
            seed: seed used in the MCERD process
            name: name of the process (usually the name of the recoil element)
        """
        # TODO add handling for fatal error messages
        # The first line that MCERD prints out is either 'MCERD is alive'
        # or 'Initializing parameters' depending on whether debug mode is on.
        # We let either one of these messages pass to observers and start
        # reducing from the next line which is:
        first_line_init = "Reading input files."
        last_line_init = "Starting simulation."
        # Note: there are quite a bit of lines between those two, some of which
        # maybe of interest for the user. Maybe implement parsing for them too.

        first_line_end = "Opening target file "
        last_line_end = "angave "

        return rx.pipe(
            ops.map(lambda x: x.decode("utf-8").strip()),
            observing.get_printer(),
            observing.reduce_while(
                reducer=str_reducer,
                start_from=lambda x: x == first_line_init,
                end_at=lambda x: x == last_line_init
            ),
            observing.reduce_while(
                reducer=str_reducer,
                start_from=lambda x: x.startswith(first_line_end),
                end_at=lambda x: x.startswith(last_line_end)
            ),
            ops.scan(lambda acc, x: {
                "presim": acc["presim"] and x != "Presimulation finished",
                **parse_raw_output(
                    x, end_at=lambda y: y.startswith(first_line_end))
            }, seed={"presim": True}),
            ops.scan(lambda acc, x: {
                "seed": seed,
                "name": name,
                "presim": x["presim"],
                "calculated": x.get("calculated", acc["calculated"]),
                "total": x.get("total", acc["total"]),
                "percentage": x.get("percentage", acc["percentage"]),
                "msg": x.get("msg", ""),
                "is_running": x.get("is_running", True)
            }, seed={"calculated": 0, "total": 0, "percentage": 0}),
            ops.take_while(lambda x: x["is_running"], inclusive=True)
        )

    @staticmethod
    def stop_process(process: subprocess.Popen):
        """Stop the MCERD process and delete the MCERD object.
        """
        if platform.system() == "Windows":
            cmd = f"TASKKILL /F /PID {process.pid} /T"
            subprocess.call(cmd)
        else:
            process.kill()

    def __create_mcerd_files(self):
        """Creates the temporary files needed for running MCERD.
        """
        # Create the main MCERD command file
        with open(self.__command_file, "w") as file:
            file.write(self.get_command_file_contents())

        # Create the MCERD detector file
        with open(self.__detector_file, "w") as file:
            file.write(self.get_detector_file_contents())

        # Create the MCERD target file
        with open(self.__target_file, "w") as file:
            file.write(self.get_target_file_contents())

        # Create the MCERD foils file
        with open(self.__foils_file, "w") as file:
            file.write(self.get_foils_file_contents())

        # Create the recoil file
        with open(self.recoil_file, "w") as file:
            file.write(self.get_recoil_file_contents())

    def get_recoil_file_contents(self):
        """Returns the contents of the recoil file.
        """
        recoil_element = self.__settings["recoil_element"]
        return "\n".join(recoil_element.get_mcerd_params())

    def get_command_file_contents(self):
        """Returns the contents of MCERD's command file as a string.
        """
        beam = self.__settings["beam"]
        target = self.__settings["target"]
        recoil_element = self.__settings["recoil_element"]
        min_scat_angle = self.__settings['minimum_scattering_angle']
        min_main_scat_angle = self.__settings['minimum_main_scattering_angle']
        min_ene_ions = self.__settings['minimum_energy_of_ions']
        rec_count = self.__settings['number_of_recoils']
        sim_mode = self.__settings['simulation_mode']
        scale_ion_count = self.__settings['number_of_scaling_ions']
        ions_in_presim = self.__settings['number_of_ions_in_presimu']

        return "\n".join([
            f"Type of simulation: {self.__settings['simulation_type']}",
            *beam.get_mcerd_params(),
            f"Target description file: {self.__target_file}",
            f"Detector description file: {self.__detector_file}",
            f"Recoiling atom: {recoil_element.element.get_prefix()}",
            f"Recoiling material distribution: {self.recoil_file}",
            f"Target angle: {target.target_theta} deg",
            "Beam spot size: " + ("%0.1f %0.1f mm" % beam.spot_size) + "",
            f"Minimum angle of scattering: {min_scat_angle} deg",
            f"Minimum main scattering angle: {min_main_scat_angle} deg",
            f"Minimum energy of ions: {min_ene_ions} MeV",
            f"Average number of recoils per primary ion: {rec_count}",
            f"Recoil angle width (wide or narrow): {sim_mode}",
            f"Presimulation * result file: {self.__presimulation_file}",
            f"Number of real ions per each scaling ion: {scale_ion_count}",
            f"Number of ions: {self.__settings['number_of_ions']}",
            f"Number of ions in the presimulation: {ions_in_presim}",
            f"Seed number of the random number generator: {self.__seed}",
        ])

    def get_detector_file_contents(self):
        """Returns the contents of the detector file as a string.
        """
        detector = self.__settings["detector"]
        foils = "\n----------\n".join("\n".join(foil.get_mcerd_params())
                                      for foil in detector.foils)

        return "\n".join([
            *detector.get_mcerd_params(),
            f"Description file for the detector foils: {self.__foils_file}",
            "==========",
            foils
        ])

    def get_target_file_contents(self):
        """Returns the contents of the target file as a string.
        """
        target = self.__settings["target"]
        cont = []
        for layer in target.layers:
            for element in layer.elements:
                cont.append(element.get_mcerd_params())

        # First layer is used for target surface calculation.
        cont += Layer.get_default_mcerd_params()

        # An indexed list of all elements is written first.
        # Then layers and their elements referencing the index.
        count = 0
        for layer in target.layers:
            cont += layer.get_mcerd_params()
            for element in layer.elements:
                cont.append(f"{count} "
                            f"{element.get_mcerd_params(return_amount=True)}")
                count += 1

        return "\n".join(cont)

    def get_foils_file_contents(self):
        """Returns the contents of the foils file.
        """
        detector = self.__settings["detector"]
        cont = []
        for foil in detector.foils:
            for layer in foil.layers:
                # Write only one layer since mcerd soesn't know how to
                # handle multiple layers in a foil
                for element in layer.elements:
                    cont.append(element.get_mcerd_params())
                break

        # An indexed list of all elements is written first.
        # Then layers and their elements referencing the index.
        count = 0
        for foil in detector.foils:
            for layer in foil.layers:
                cont += layer.get_mcerd_params()
                for element in layer.elements:
                    cont.append(
                        f"{count} "
                        f"{element.get_mcerd_params(return_amount=True)}")
                    count += 1
                break

        return "\n".join(cont)

    def copy_results(self, destination):
        """Copies MCERD result file (.erd) and recoil file into given
        destination.

        Args:
            destination: Destination folder.
        """
        shutil.copy(self.result_file, destination)
        shutil.copy(self.recoil_file, destination)

    def delete_unneeded_files(self):
        """
        Delete mcerd files that are not needed anymore.
        """
        def delete_files(*files):
            for f in files:
                try:
                    f.unlink()
                except OSError:
                    pass

        delete_files(
            self.__command_file, self.__detector_file, self.__target_file,
            self.__foils_file)

        def filter_func(f):
            return f.startswith(self.__rec_filename)

        gf.remove_files(
            self.sim_dir, exts={".out", ".dat", ".range", ".pre"},
            filter_func=filter_func)


_pattern = re.compile(r"Calculated (?P<calculated>\d+) of (?P<total>\d+) ions "
                      r"\((?P<percentage>\d+)%\)")


def parse_raw_output(raw_line, end_at=None):
    """Parses raw output produced by MCERD into something meaningful.
    """
    m = _pattern.match(raw_line)
    try:
        return {
            "calculated": int(m.group("calculated")),
            "total": int(m.group("total")),
            "percentage": int(m.group("percentage"))
        }
    except AttributeError:
        if raw_line == "Presimulation finished":
            return {
                "calculated": 0,
                "percentage": 0,
                "msg": raw_line
            }
        elif end_at is not None and end_at(raw_line):
            return {
                "msg": raw_line,
                "percentage": 100,
                "is_running": False
            }
        return {
            "msg": raw_line
        }


def str_reducer(acc, x):
    """Helper function for reducing strings from multiple lines.
    Appends a newline and x to the previously accumulated string.
    """
    return f"{acc}\n{x}"
