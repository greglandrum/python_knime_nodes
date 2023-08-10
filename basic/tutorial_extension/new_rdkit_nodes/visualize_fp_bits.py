# -*- coding: utf-8 -*-
# ------------------------------------------------------------------------
#  Copyright by KNIME AG, Zurich, Switzerland
#  Website: http://www.knime.com; Email: contact@knime.com
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
Part of the RDKit Python extension. Node 'Visualize Morgan fingerprint bits'.

@author Alice Krebs, KNIME GmbH, Konstanz, Germany
@author Steffen Fissler, KNIME GmbH, Konstanz, Germany
"""

import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from new_rdkit_nodes import utils

from PIL import Image
from io import BytesIO

LOGGER = logging.getLogger(__name__)
IPythonConsole.UninstallIPythonRenderer()


class visualizefpbits(knext.PythonNode):
    """
    This node has a description, and I will change it once I figured out the code...
    """

    number_bits = knext.IntParameter("Number of bits",
                                     "Define the number of bits",
                                     1024,
                                     min_value=0)
    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description=
        "Choose the column from the first input table containing the molecules",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False,
    )
    bits_column_param = knext.ColumnParameter(
        label="Bits column",
        description=
        "Choose the column from the second input table containing the bits as integer",
        port_index=1,
        column_filter=lambda column: column.ktype in
        (knext.int64(), knext.int32()),
        include_row_key=False,
        include_none_column=False,
    )

    def configure(self, configure_context, input_schema_1: knext.Schema,
                  input_schema_2):
        return input_schema_1

    def draw_molecule_with_bit(self, mol, idx, bi):
        raise NotImplementedError("needs to be defined in derived class")

    def init_generator(self):
        raise NotImplementedError("needs to be defined in derived class")

    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table, input_2: knext.Table):
        if self.molecule_column_param is None or self.bits_column_param is None:
            raise AttributeError(
                "Molecule or Bits column was not selected in configuration dialog."
            )

        fpgen, ao = self.init_generator()

        df = input_1.to_pandas()

        molecule_column_type = input_1.schema[self.molecule_column_param].ktype
        df = input_1.to_pandas()
        mols = utils.convert_column_to_rdkit_mol(df,
                                                 molecule_column_type,
                                                 self.molecule_column_param,
                                                 sanitizeOnParse=True)

        # create list of FP bits from input table 2
        fp_ids = input_2.to_pandas()[self.bits_column_param]

        cols = [None] * len(fp_ids)
        for i in range(len(cols)):
            cols[i] = []

        for mol in mols:
            fp = fpgen.GetFingerprint(mol, additionalOutput=ao)
            Chem.Kekulize(mol)  # kekulize molecules
            for i, idx in enumerate(
                    fp_ids
            ):  # if rendering fails, append an empty cell. Don't make if-else to catch the error
                if fp[idx]:
                    try:
                        # img = Draw.DrawMorganBit(mol, idx, bi, useSVG=True)
                        img = self.draw_molecule_with_bit(mol, idx, ao)
                        # sio = BytesIO(img)
                        # img = Image.open(sio)
                        cols[i].append(img)
                    except:
                        import traceback
                        traceback.print_exc()
                        cols[i].append(None)
                else:
                    cols[i].append(None)

        for i, idx in enumerate(fp_ids):
            df[f"bit{idx}"] = cols[i]

        return knext.Table.from_pandas(df)
