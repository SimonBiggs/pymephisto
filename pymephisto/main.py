# Copyright (C) 2015 Simon Biggs
# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public
# License along with this program. If not, see
# http://www.gnu.org/licenses/.

import os

import numpy as np

from .mccread import (
    pull_mephisto_item, pull_mephisto_number, pull_mephisto_data)
from .csvoutput import file_output


def load_mephisto(filepath, output_to_file=True, output_directory=None,
                  sort=True):
    """Input the filepath of a mephisto .mcc file and return the data of the
    scans in four lists, distance, relative_dose, scan_curvetype, and
    scan_depth. Each respective element in these lists corresponds to an
    individual scan.
    """
    # Open the file and store the contents in file_contents
    with open(filepath) as file_pointer:
            file_contents = np.array(file_pointer.readlines())

    # Use the functions defined within mccread.py to pull the desired data
    distance, relative_dose = pull_mephisto_data(file_contents)
    scan_curvetype = pull_mephisto_item('SCAN_CURVETYPE', file_contents)
    scan_depth = pull_mephisto_number('SCAN_DEPTH', file_contents)

    if sort:
        sort_ref = np.hstack([
            np.where(scan_curvetype == 'PDD')[0],
            np.where(scan_curvetype == 'INPLANE_PROFILE')[0],
            np.where(scan_curvetype == 'CROSSPLANE_PROFILE')[0]
        ])

        assert len(sort_ref) == len(scan_curvetype)

        distance = list(np.array(distance)[sort_ref])
        relative_dose = list(np.array(relative_dose)[sort_ref])
        scan_curvetype = list(np.array(scan_curvetype)[sort_ref])
        scan_depth = list(np.array(scan_depth)[sort_ref])

    # Output csv's if "output_to_file" is True
    if output_to_file:

        # If user didn't define an output_directory use a default one
        if output_directory is None:
            # Define output directory as a mephisto folder
            filepath_directory = os.path.dirname(filepath)
            filename = os.path.splitext(os.path.basename(filepath))[0]
            output_directory = os.path.join(
                filepath_directory, filename)

        # If the output directory does not exist create it
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        file_output(
            output_directory, distance, relative_dose,
            scan_curvetype, scan_depth)

    return distance, relative_dose, scan_curvetype, scan_depth
