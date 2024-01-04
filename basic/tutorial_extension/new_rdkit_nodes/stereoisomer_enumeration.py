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
Part of the RDKit Python extension. Node 'Enumerate Stereoisomers'.
@author Greg Landrum, ETH Zurich, Switzerland
@author Alice Krebs, KNIME AG, Konstanz, Germany
"""

import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import EnumerateStereoisomers
import knime.types.chemistry as cet
from new_rdkit_nodes import utils
import pandas as pd
import numpy as np
LOGGER = logging.getLogger(__name__)


@knext.node(
    name="Stereoisomer Enumeration",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icons/StereoisomerEnumeration.png",
    category=utils.category
    )
@knext.input_table(
    name="Input Table", 
    description="Input data with molecules"
    )
@knext.output_table(
    name="Output Table",
    description=""
    )
class StereoisomerEnumeration(knext.PythonNode):
    """
    This node does stereoisomer enumeration on the selected input molecules.

    This node does stereoisomer enumeration on the selected input molecules. 
    """

    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description="Select the molecule column to standardize. The column has to be SMILES, SDF, or RDKit molecule.",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False)  
    
    identifier_column_param = knext.ColumnParameter(
        label="Identifier column",
        description="Select the column with the unique identifier. The column has to be string or numeric.",
        port_index=0,
        include_row_key=True,
        include_none_column=False) 
    
    onlyunassigned = knext.BoolParameter(
        label="Unspecified centers only",
        description="If set (the default), stereocenters which have specified stereochemistry will not be perturbed unless they are part of a relative stereo group.",
        default_value=True
        )
    
    unique = knext.BoolParameter(
        label="Unique stereoisomers only",
        description="Remove duplicate stereoisomers.",
        default_value=True
        )
    
    onlystereogroups = knext.BoolParameter(
        label="Stereo groups only",
        description="Only enumerate stereocenters involved in stereo groups (enhanced stereochemistry).",
        default_value=False
        )

    tryembedding = knext.BoolParameter(
        label="tryEmbedding",
        description="The process attempts to generate a standard RDKit distance geometry conformation for the stereisomer (computationally expensive).",
        default_value=False
        )
    
    def configure(self, config_context, input_schema_1: knext.Schema):
        identifier_column_type = get_ktype_for_column(input_schema_1, self.identifier_column_param)
        index_type = knext.int32()
        molecule_column_type = knext.logical(Chem.rdchem.Mol)
        schema_1 = knext.Schema([identifier_column_type, index_type, molecule_column_type], [input_schema_1[self.identifier_column_param].name, "index", "stereoisomer"])
        return schema_1
 
    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table):
        if self.molecule_column_param is None:
            raise AttributeError(
            "Molecule column was not selected in configuration dialog."
            )
        
        if self.identifier_column_param is None:
            raise AttributeError(
            "Identifier column was not selected in configuration dialog."
            )
        
        molecule_column_type = input_1.schema[self.molecule_column_param].ktype
        identifier_column_name = input_1.schema[self.identifier_column_param].name

        df = input_1.to_pandas()
    
        mols = utils.convert_column_to_rdkit_mol(df,
                                                molecule_column_type,
                                                self.molecule_column_param,
                                                sanitizeOnParse=True)
        
        esopt = EnumerateStereoisomers.StereoEnumerationOptions(tryEmbedding=self.tryembedding, onlyUnassigned=self.onlyunassigned, maxIsomers=1024, rand=None, unique=self.unique, onlyStereoGroups=self.onlystereogroups)
        
        progress = 0.0
        add_to_progress = 1 / input_1.num_rows
        results = []
        for id, mol in zip(df[self.identifier_column_param], mols): 
            if mol is None:
                continue
            enum_res = EnumerateStereoisomers.EnumerateStereoisomers(mol, options=esopt)
            for i, iso in enumerate(enum_res):
                results.append((id, i, iso))
            progress += add_to_progress
            exec_context.set_progress(progress=progress)
        df_results = pd.DataFrame(results)
        df_results.rename(columns={df_results.columns[0]: identifier_column_name, df_results.columns[1]: 'index', df_results.columns[2]: 'stereoisomer'}, inplace=True)
        
        if df_results[identifier_column_name].dtypes == 'int64': 
            df_results = df_results.astype({identifier_column_name: 'int32'})
        df_results = df_results.astype({'index': 'int32'})

        return knext.Table.from_pandas(df_results)
    
# Returns the ktype of the first column which has the col_name as name
def get_ktype_for_column(schema: knext.Schema, col_name: str):
        for column in schema._columns:
            if column.name == col_name:
                return column.ktype
            