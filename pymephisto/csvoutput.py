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
import pandas as pd


def file_output(output_directory, distance, relative_dose,
                scan_curvetype, scan_depth):
    """Store the loaded mephisto data into csv files for easy user confirmation
    and use.
    """
    # Determines the filepaths for the output
    filepaths = determine_output_filepaths(
        output_directory, scan_curvetype, scan_depth)

    columns = ['distance (mm)', 'relative dose']

    # Loop over each curvetype and save the data to csv
    for i, curvetype in enumerate(scan_curvetype):
        # Stacks the data into one array and transposes into column orientation
        data = np.vstack([distance[i], relative_dose[i]]).T

        # Use pandas to save data to csv
        df = pd.DataFrame(data, columns=columns)
        df.to_csv(filepaths[i])


def determine_output_filepaths(output_directory, scan_curvetype, scan_depth):
    """Determine a useful filepath for the saving of each mephisto scan.
    """
    filepaths = []

    # Loop over each scan curvetype creating a relevant filepath
    for i, curvetype in enumerate(scan_curvetype):
        if curvetype == 'PDD':
            # Create the filename to be pdd_[number].csv
            filepaths.append(os.path.join(
                output_directory, "pdd_[{0:d}].csv".format(i)))

        elif curvetype == 'INPLANE_PROFILE':
            # Create the filename to be inplaneprofile_depth_[number].csv
            filepaths.append(os.path.join(
                output_directory,
                "inplaneprofile_{0:d}mm_[{1:d}].csv".format(
                    int(scan_depth[i]), i)))

        elif curvetype == 'CROSSPLANE_PROFILE':
            # Create the filename to be crossplaneprofile_depth_[number].csv
            filepaths.append(os.path.join(
                output_directory,
                "crossplaneprofile_{0:d}mm_[{1:d}].csv".format(
                    int(scan_depth[i]), i)))

        else:
            # Raise an error if the curve type was not as expected
            raise Exception("Unexpected scan_curvetype")

    return filepaths
