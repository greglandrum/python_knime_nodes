# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------
#  Copyright by Greg Landrum
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License, Version 3, as
#  published by the Free Software Foundation.
#
#  This program is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <http://www.gnu.org/licenses>.
#
#  Additional permission under GNU GPL version 3 section 7:
#
#  KNIME interoperates with ECLIPSE solely via ECLIPSE's plug-in APIs.
#  Hence, KNIME and ECLIPSE are both independent programs and are not
#  derived from each other. Should, however, the interpretation of the
#  GNU GPL Version 3 ("License") under any applicable laws result in
#  KNIME and ECLIPSE being a combined program, KNIME AG herewith grants
#  you the additional permission to use and propagate KNIME together with
#  ECLIPSE with only the license terms in place for ECLIPSE applying to
#  ECLIPSE and the GNU GPL Version 3 applying for KNIME, provided the
#  license terms of ECLIPSE themselves allow for the respective use and
#  propagation of ECLIPSE together with KNIME.
#
#  Additional permission relating to nodes for KNIME that extend the Node
#  Extension (and in particular that are based on subclasses of NodeModel,
#  NodeDialog, and NodeView) and that only interoperate with KNIME through
#  standard APIs ("Nodes"):
#  Nodes are deemed to be separate and independent programs and to not be
#  covered works.  Notwithstanding anything to the contrary in the
#  License, the License does not apply to Nodes, you are not required to
#  license Nodes under the License, and you are granted a license to
#  prepare and propagate Nodes, in each case even if such Nodes are
#  propagated with or for interoperation with KNIME.  The owner of a Node
#  may freely choose the license terms applicable to such Node, including
#  when such Node is propagated with or for interoperation with KNIME.
# ------------------------------------------------------------------------
"""
Part of the RDKit Python extension. Node 'Normalize Molecule'.

@author Greg Landrum, ETH Zurich, Switzerland
"""

import logging
import knime.extension as knext
from rdkit import Chem
import pandas as pd
import knime.types.chemistry as cet  # To work with and compare against chemical data types like SMILES,...
import pyarrow as pa
from new_rdkit_nodes import utils
#import knime_arrow_pandas  # TODO Refactor once ticket AP-19209 is implemented

from rdkit.Chem.MolStandardize import rdMolStandardize

# The standardization code is very verbose, so disable the info log
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.info")


@knext.node(
    name="Normalize Molecule",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icon.png",
    category=utils.category,
)
@knext.input_table(name="Input Data", description="Input data with molecules")
@knext.output_table(name="Output Data",
                    description="Input data with normalized molecules")
class GetNormalizedMoleculeNode(knext.PythonNode):
    """Returns normalized structure of a molecule.

    This node returns the normalized structure of an input molecule.
    """

    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description=
        "Select the molecule column to normalize. The column has to be SMILES, SDF, or RDKit molecule.",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1: knext.Schema):
        return input_schema_1.append(
            knext.Column(ktype=Chem.Mol, name="Parent Molecule"))

    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table):
        if self.molecule_column_param is None:
            raise AttributeError(
                "Molecule column was not selected in configuration dialog.")

        # Prepare mols: If input column consists of rdkit molecules, use them;
        # if input column consists of SMILES, convert to rdkit molecules
        molecule_column_type = input_1.schema[self.molecule_column_param].ktype
        df = input_1.to_pandas()
        mols = utils.convert_column_to_rdkit_mol(df,
                                                 molecule_column_type,
                                                 self.molecule_column_param,
                                                 sanitizeOnParse=False)

        progress = 0.0
        add_to_progress = 1 / input_1.num_rows
        pmols = []
        normalizer = rdMolStandardize.Normalizer()
        for mol in mols:
            mol = normalizer.normalize(mol)
            pmols.append(mol)
            progress += add_to_progress
            exec_context.set_progress(progress=progress)

        # Add the parent molecule to the pandas dataframe as a new column
        df["Parent Molecule"] = pmols
        # Convert the processed table back from a pandas DataFrame to a KNIME table
        return knext.Table.from_pandas(df)
