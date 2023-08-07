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
Part of the RDKit Python extension. Node 'Visualize RDKit fingerprint bits'.

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
from new_rdkit_nodes.visualize_fp_bits import visualizefpbits
from PIL import Image
from io import BytesIO

LOGGER = logging.getLogger(__name__)
IPythonConsole.UninstallIPythonRenderer()


@knext.node(
    name="Visualize bits of RDKit fingerprints",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icon.png",
    category=utils.category)
@knext.input_table(name="Input table 1", description="Input table 1 with molecules")
@knext.input_table(name="Input table 2", description="Input table 2 with RDKit fingerprint bits")
@knext.output_table(
    name="Highlighted bits",
    description="Output tables including images of the highlighted bits",
)

class visualizerdkitfpbits(visualizefpbits):
    """
    This node will take a first input table with molecules and a second input table with bits and hightlight them in the output table, which is still a shitty description. 
    """
    def get_fp(self, mol, maxPath, fpSize, bitInfo):
        Chem.RDKFingerprint(mol, maxPath, fpSize, bitInfo)
    
    def draw_molecule_with_bit(self, mol, idx, bi):
        img = Draw.DrawRDKitBit(mol, idx, bi)
        return img


# class visualizerdkitfpbits:
#     """
#     This node will take a first input table with molecules and a second input table with bits and hightlight them in the output table, which is still a shitty description. 
#     """

#     fp_size = knext.IntParameter(
#         "size of fingerprint",
#         "Define the fingerprint size aka number of bits",
#         2048,
#         min_value=0,
#     )
#     # min_path = knext.IntParameter("minimum path length", "Define the path length", 2, min_value=0)
#     max_path = knext.IntParameter(
#         "maximum path length", "Define the path length", 2, min_value=0
#     )

#     # def is_molecule(column):  # Filter columns visible in the column_param for chemistry type; column_filter=is_molecule
#     #     return (
#     #         column.ktype == knext.smiles()
#     #         or column.ktype == knext.smarts()
#     #         or column.ktype == knext.sdf()
#     #     )

#     molecule_column = knext.ColumnParameter(
#         label="Molecule column",
#         description="Choose the column from the first input table containing the molecules",
#         port_index=0,
#     )
#     bits_column = knext.ColumnParameter(
#         label="Bits column",
#         description="Choose the column from the second input table containing the bits as integer",
#         port_index=1,
#     )

#     def configure(self, configure_context, input_schema_1: knext.Schema,
#                   input_schema_2):
#         return input_schema_1


#     def execute(self, exec_context, input_1, input_2):
#         a = len(self.bits_column)

#         input_1_pandas = input_1.to_pandas()
#         input_2_pandas = input_2.to_pandas()

#         mols = []
#         for smi in input_1_pandas[self.molecule_column]:
#             mols.append(Chem.MolFromSmiles(smi))


#         # create list of FP bits from input table 2
#         fp_ids = input_2_pandas[self.bits_column]

#         cols = [None] * len(fp_ids)
#         for i in range(len(cols)):
#             cols[i] = []

#         # defining the output table
#         output_table_1 = input_1_pandas.copy()

#         for mol in mols:
#             rdkbi = {}  # defining a dictionary
#             fp = Chem.RDKFingerprint(
#                 mol, maxPath=self.max_path, fpSize=self.fp_size, bitInfo=rdkbi
#             )  # calculate fingerprint with user-defined path length and nr of bits aka fp size
#             Chem.Kekulize(mol)  # kekulize molecules
#             for i, idx in enumerate(
#                 fp_ids
#             ):  # if rendering fails, append an empty cell. Don't make if-else to catch the error
#                 if fp[idx]:
#                     try:
#                         cols[i].append(Draw.DrawRDKitBit(mol, idx, rdkbi))
#                     except:
#                         cols[i].append(None)
#                 else:
#                     cols[i].append(None)

#         for i, idx in enumerate(fp_ids):
#             output_table_1[f"bit{idx}"] = cols[i]

#         return knext.Table.from_pandas(output_table_1)
