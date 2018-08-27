# ----------------------------------------------------------------------------
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import tempfile
import hashlib
import subprocess
import pkg_resources

import biom
import skbio
import qiime2.util
import pandas as pd
import q2templates
import qiime2

from q2_types.feature_data import TaxonomyFormat, TSVTaxonomyFormat
from qiime2 import Metadata

TEMPLATES = pkg_resources.resource_filename('q2_corncob', 'assets')

def run_commands(cmds, verbose=True):
    if verbose:
        print("Running external command line application(s). This may print "
              "messages to stdout and/or stderr.")
        print("The command(s) being run are below. These commands cannot "
              "be manually re-run as they will depend on temporary files that "
              "no longer exist.")
    for cmd in cmds:
        if verbose:
            print("\nCommand:", end=' ')
            print(" ".join(cmd), end='\n\n')
        subprocess.run(cmd, check=True)

def differentialtest(table: biom.Table, metadata: qiime2.Metadata,
        variable: str, taxonomy: TSVTaxonomyFormat) -> pd.DataFrame:

    if table.is_empty():
        raise ValueError("The provided table object is empty")
    ## run the R script on the file
    with tempfile.TemporaryDirectory() as temp_dir_name:
        ## write the biom table to file
        input_table = os.path.join(temp_dir_name, 'table.tsv')
        input_metadata = os.path.join(temp_dir_name, 'metadata.tsv')

        with open(input_table, 'w') as fh:
            fh.write(table.to_tsv())
        metadata.save(input_metadata)

        output = os.path.join(temp_dir_name, 'data.tsv')

        cmd = ['differentialtest.R', input_table, input_metadata, str(taxonomy), str(variable), str(output)]
        run_commands([cmd])
        data = pd.read_csv(output, sep = '\t')
        data.index.name = 'Feature ID'
    return data
