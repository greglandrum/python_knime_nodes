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
Part of the RDKit Python extension. Node 'Mols to Grid Image'.
@author Greg Landrum, ETH Zurich, Switzerland
@author Alice Krebs, KNIME AG, Konstanz, Germany
"""

import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from new_rdkit_nodes import utils
from PIL import Image
from io import BytesIO
LOGGER = logging.getLogger(__name__)
IPythonConsole.UninstallIPythonRenderer()

@knext.node(
    name="Molecules to Image Grid",
    node_type=knext.NodeType.VISUALIZER,
    icon_path="icons/MolsToGridImage.png",
    category=utils.category
    )
@knext.input_table(
    name="Input Data", 
    description="Input data with molecules."
    )
@knext.output_image(
    name="Output Image",
    description="Image with molecules arranged in the defined grid."
    )
class MolsToGrid:
    """
    This node outputs a grid of molecule images as PNG or SVG.  
    
    This node outputs a grid of molecule images as PNG or SVG. 
    """
    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description="Select the molecule column to render the images from. The column has to be SMILES, SDF, or RDKit molecule.",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False)

    image_label_column = knext.ColumnParameter(label="Column with image label", 
                                                description="Select a column from the input table that should serve as image label.", 
                                                port_index=0)

    mols_per_row = knext.IntParameter(label="Molecules per row", 
                                        description="Set the number of molecule images per row.", 
                                        default_value=3, 
                                        min_value=1)

    image_width = knext.IntParameter(label="Image width", 
                                        description="Define the overall image width in pixels. The default is 200.", 
                                        default_value=200, 
                                        min_value=1)

    image_height = knext.IntParameter(label="Image height", 
                                        description="Define the overall image height in pixels. The default is 200.", 
                                        default_value=200, 
                                        min_value=1)
    output_type = knext.StringParameter(label="Output image type",
                                        description="Select if the returned output image type should be PNG or SVG",
                                        default_value="SVG",
                                        enum = ["SVG", "PNG"]
                                        )
    scale = knext.BoolParameter(label="Draw molecules on same scale",
                                description="Select if the molecules should be drawn on the same scale taking the largest molecule of the input data as reference. Default value is true",
                                default_value=True
                                )

    def configure(self, config_context, input_schema_1: knext.Schema):
        if self.output_type == "PNG": 
            return (knext.ImagePortObjectSpec(knext.ImageFormat.PNG))
        else: 
            return (knext.ImagePortObjectSpec(knext.ImageFormat.SVG))

 
    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table):
    
        if self.molecule_column_param is None:
            raise AttributeError(
            "Molecule column was not selected in configuration dialog."
            )
        
        molecule_column_type = input_1.schema[self.molecule_column_param].ktype 
        
        df = input_1.to_pandas()
    
        mols = utils.convert_column_to_rdkit_mol(df,
                                                molecule_column_type,
                                                self.molecule_column_param,
                                                sanitizeOnParse=True)
        dOpts = Draw.MolDrawOptions()
        dOpts.drawMolsSameScale=self.scale

        if self.output_type == "PNG":
            SVG_bool = False

        elif self.output_type == "SVG":
            SVG_bool = True
        
        img = Draw.MolsToGridImage(mols, 
                                   molsPerRow=self.mols_per_row, 
                                   subImgSize=(self.image_width,self.image_height), 
                                   legends=[str(x) for x in df[self.image_label_column]], 
                                   useSVG=SVG_bool,
                                   returnPNG=True,
                                   drawOptions=dOpts
                                   )

        return img

