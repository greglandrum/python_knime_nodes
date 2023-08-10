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
@author Greg Landrum, ETH Zurich, Zurich, Switzerland
@author Alice Krebs, KNIME GmbH, Konstanz, Germany
"""

import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from new_rdkit_nodes import utils
from new_rdkit_nodes.visualize_fp_bits import visualizefpbits
from PIL import Image
from io import BytesIO

LOGGER = logging.getLogger(__name__)
IPythonConsole.UninstallIPythonRenderer()


@knext.node(name="Visualize Morgan fingerprint bits",
            node_type=knext.NodeType.MANIPULATOR,
            icon_path="icon.png",
            category=utils.category)
@knext.input_table(name="Input table 1",
                   description="Input table 1 with molecules")
@knext.input_table(name="Input table 2",
                   description="Input table 2 with Morgan fingerprint bits")
@knext.output_table(
    name="Highlighted bits",
    description="Output tables including images of the highlighted bits",
)
class visualizemorganfpbits(visualizefpbits):
    radius = knext.IntParameter("FP radius",
                                "Define the number of the FP radius",
                                2,
                                min_value=0)

    def init_generator(self):
        self._generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=self.radius, fpSize=self.number_bits)
        ao = rdFingerprintGenerator.AdditionalOutput()
        ao.AllocateBitInfoMap()
        return self._generator, ao

    def draw_molecule_with_bit(self, mol, idx, additionalOutput):
        bi = additionalOutput.GetBitInfoMap()
        img = Draw.DrawMorganBit(mol, idx, bi, useSVG=True)
        return img
