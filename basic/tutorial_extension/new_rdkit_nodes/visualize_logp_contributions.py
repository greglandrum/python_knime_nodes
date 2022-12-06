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
Part of the RDKit Python extension. Node 'Visualize atomic LogP contributions'.

@author Martyna Pawletta, KNIME GmbH, Konstanz, Germany
@author Steffen Fissler, KNIME GmbH, Konstanz, Germany
"""

import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem import rdMolDescriptors
from new_rdkit_nodes.utils import category
#import chem_extension_types as cet  # To work with and compare against chemical data types like SMILES,...


LOGGER = logging.getLogger(__name__)


@knext.node(name="Visualize atomic LogP contributions", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category)
@knext.input_table(name="Input Data", description="Table containing a SMILES column")
@knext.output_table(name="Output Data", description="Input including LogP contributions")
class LogPContribNode:
    """Visualizes atomic LogP contributions.

    This node adds a column to the output table containing visualized 
    contributions to LogP with the RDKit molecules.
    """

    def column_is_smiles_or_mol(column: knext.Column):
        c_type = column.ktype
        smiles_type = knext.logical(cet.SmilesValue)
        smiles_adapter_type = knext.logical(cet.SmilesAdapterValue)
        rdkit_mol_type = knext.logical(Chem.rdchem.Mol)
        return (
            c_type == smiles_type
            or c_type == smiles_adapter_type
            or c_type == rdkit_mol_type
        )


    molecule_column = knext.ColumnParameter(label="Select column containing SMILES data", description="SMILES column needs to be an RDKit Canonical SMILES Column Type", column_filter=column_is_smiles_or_mol)


    def configure(self, configure_context, input_schema_1: knext.Schema):
        return input_schema_1.append(knext.Column(knext.string(), "ALogP_contributions"))
 


    def execute(self, exec_context, input_1: knext.Table):
        if self.core_col is None or self.molecule_column is None:
            raise AttributeError(
               "SMILES column was not selected in configuration dialog."
            )


        # Prepare mols: If input column consists of rdkit molecules, use them;
        # if input column consists of SMILES, convert to rdkit molecules
        molecule_column_type = input_1.schema[self.molecule_column_param].ktype
        mols = []
        df = input_1.to_pandas()
        if molecule_column_type == knext.logical(Chem.rdchem.Mol):
            LOGGER.warning("rdkit mols detected")
            mols = df[self.molecule_column_param]
        else:
            for smi in df[self.molecule_column_param]:
                mols.append(Chem.MolFromSmiles(smi, sanitize=False))

        svg = []

        progress = 0.0
        add_to_progress = 1 / input_1.num_rows

        for mol in mols:
            d = Draw.MolDraw2DSVG(600, 600)
            contribs = rdMolDescriptors._CalcCrippenContribs(mol)
            SimilarityMaps.GetSimilarityMapFromWeights(mol,[x for x,y in contribs], draw2d=d, colorMap="PiYG")
            d.FinishDrawing()
            svg.append(d.GetDrawingText())
            progress = progress + add_to_progress
            exec_context.set_progress(progress=progress)

            
        df['ALogP_contributions'] = svg
        return knext.Table.from_pandas(df)


       