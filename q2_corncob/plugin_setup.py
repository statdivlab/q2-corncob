# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import importlib

import qiime2.plugin
from qiime2.plugin import (Plugin, Str, Properties, Choices, Int, Bool, Range, Float, Set, Visualization, Metadata, MetadataColumn, Categorical, Numeric, Citation, SemanticType)

import q2_types
from q2_types.sample_data import SampleData
from q2_types.feature_table import FeatureTable, Frequency, Composition
from q2_types.feature_data import FeatureData, Taxonomy, TSVTaxonomyFormat

import q2_corncob

citations = qiime2.plugin.Citations.load('citations.bib', package='q2_corncob')



plugin = qiime2.plugin.Plugin(
    name='corncob',
    version='0.1.0',
    website='http://bryandmartin.github.io/corncob/',
    package='q2_corncob',
    description=('This QIIME 2 plugin wraps corncob and supports '
                 'single-taxon regression using the corncob R library.'),
    short_description='Plugin for single-taxon regression with corncob.',
    citations=[citations['martin2018']]
)

plugin.methods.register_function(
    function=q2_corncob.differentialtest,
    inputs={'table': FeatureTable[Frequency],
            'taxonomy': FeatureData[Taxonomy]
    },
    parameters={'metadata': Metadata,
                'variable': Str,
    },
    
    outputs=[('output',FeatureData[Taxonomy % Properties(["Taxon", "DA", "DV"])])],
    input_descriptions={'table': ('A feature table.'),
                        'taxonomy': ('Your taxonomic classification by unique feature')
    },
                                 
    output_descriptions={'output': 'A table of FDR <0.05 p-values.'
    },
    parameter_descriptions={'metadata': ('Your metadata'),
                            'variable': ('A categorical variable in your metadata')},
    name='Run differential test',
    description='This method runs differential test.'
)

