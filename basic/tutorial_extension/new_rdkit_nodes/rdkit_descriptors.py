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
Part of the RDKit Python extension. Node 'Descriptors'.
@author Greg Landrum, ETH Zurich, Switzerland
@author Alice Krebs, KNIME AG, Konstanz, Germany
"""
import os
import sys
import logging
import knime_extension as knext
import pandas as pd
from rdkit import Chem
from rdkit.Chem import  Descriptors
from new_rdkit_nodes import utils
# raise ValueError(os.environ)
# sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'Library', 'share','RDKit','Contrib'))
# import SA_Score
# from SA_Score import sascorer
# import NP_Score
# from NP_Score import npscorer
LOGGER = logging.getLogger(__name__)


@knext.node(
    name="RDKit Descriptor Calculator",
    node_type=knext.NodeType.MANIPULATOR,
    icon_path="icon.png",
    category=utils.category
    )
@knext.input_table(
    name="Input Table", 
    description="Input data with molecules"
    )
@knext.output_table(
    name="Output Table",
    description="Input table with an additional column for each descriptor"
    )
class ExtensiveDescriptorCalculator:
    """ This new descriptor calculator node calculates almost all descriptors available in the RDKit Python package (in comparison to the existing node that only implemented a selection of them. 
    It also includes the SA and NP score, as well as the TPSA including polar sulfur and phosphate. It does not calculare the Ipc descriptor. 
    """

    molecule_column_param = knext.ColumnParameter(
        label="Molecule column",
        description="Select the molecule column to calculate the descriptors on. The column has to be SMILES, SDF, or RDKit molecule.",
        port_index=0,
        column_filter=utils.column_is_convertible_to_mol,
        include_row_key=False,
        include_none_column=False)   

    def configure(self, config_context, input_schema_1: knext.Schema):
        mol_dict = getDescriptorDataTypes()
        for col_name in mol_dict: 
            if mol_dict[col_name] is int:
                ktype = knext.int32()
            else:
                ktype = knext.double()
            input_schema_1.append(knext.Column(ktype,col_name))
        return input_schema_1
 
    def execute(self, exec_context: knext.ExecutionContext,
                input_1: knext.Table):
        LOGGER.warning(os.environ)
        for i in os.environ:
            LOGGER.warning(i)
        sys.path.append(os.path.join(os.environ['CONDA_PREFIX'],'Library', 'share','RDKit','Contrib'))
        return None
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
        # calculate RDKit descriptors excluding Ipc
        descriptors = {key:[] for key in getDescriptorDataTypes().keys()}

        for mol in mols:             
            if mol is None: 
                for k in descriptors.keys():
                    descriptors[k].append(None)
            else: 
                descrs = getMolDescriptors(mol)
                for k,v in descrs.items():
                    descriptors[k].append(v)
        
        for nm,col in descriptors.items():
            df[nm] = col

        # calculate the NP score
        # fscore = npscorer.readNPModel()
        # NP = []
        # for mol in mols:
        #     if mol is None: 
        #        NP.append(None) 
        #     else: 
        #         NP.append(npscorer.scoreMol(mol,fscore))
        # df['NP_Score'] = NP

        # calculate the SA score
        # SA = []
        # for mol in mols:
        #     if mol is None: 
        #         SA.append(None)
        #     else: 
        #         SA.append(sascorer.calculateScore(mol))
        # df['SA_Score'] = SA

        # calculate TPSA including polar S and P
        TPSA_SP = []
        for mol in mols:
            if mol is None: 
                TPSA_SP.append(None)
            else: 
                TPSA_SP.append(Descriptors.TPSA(mol), includeSandP=True)
        df['TPSA_includeSandP'] = TPSA_SP 

        return knext.Table.from_pandas(df)

def getMolDescriptors(mol, missingVal=None):
    res = {}
    for nm,fn in Descriptors._descList:
        if nm == "Ipc":
            continue
        try:
            val = fn(mol)
        except:
            import traceback
            traceback.print_exc()
            val = missingVal
        res[nm] = val
    # calculate SA, NP AND TPSA here
    return res

def getDescriptorDataTypes(missingVal=None):
    res = {}
    mol = Chem.MolFromSmiles('CCO')
    for nm,fn in Descriptors._descList:
        if nm == "Ipc":
            continue
        try:
            tp = type(fn(mol))
        except:
            import traceback
            traceback.print_exc()
            tp = missingVal
        res[nm] = tp
    return res

#################
# TODO: 
# progress bar -> how to implement best? 
# import SA and NP Score crashes node (not shown in node repo) -> KeyError: 'CONDA_PREFIX', where is conda install, ask python kids