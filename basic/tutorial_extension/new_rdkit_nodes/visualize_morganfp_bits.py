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
from turtle import clear
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from new_rdkit_nodes.utils import category

LOGGER = logging.getLogger(__name__)
IPythonConsole.UninstallIPythonRenderer()

#import chem_extension_types as cet  # To work with and compare against chemical data types like SMILES,...


@knext.node(name="Visualize bits of Morgan fingerprints", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category)
@knext.input_table(name="Input table 1", description="Input table 1 with molecules")
@knext.input_table(name="Input table 2", description="Input table 2 with Morgan fingerprint bits")
@knext.output_table(
    name="Highlighted bits",
    description="Output tables including images of the highlighted bits",
)
class visualizemorganfpbits(knext.PythonNode):
    """
    This node has a description, and I will change it once I figured out the code...
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

    number_bits = knext.IntParameter(
        "Number of bits", "Define the number of bits", 1024, min_value=0
    )
    radius = knext.IntParameter(
        "FP radius", "Define the number of the FP radius", 2, min_value=0
    )

    molecule_column = knext.ColumnParameter(
        label="Molecule column",
        description="Choose the column from the first input table containing the molecules",
        port_index=0,
        column_filter=column_is_smiles_or_mol,
    )
    bits_column = knext.ColumnParameter(
        label="Bits column",
        description="Choose the column from the second input table containing the bits as integer",
        port_index=1,
        column_filter=lambda column: column.ktype == knext.int64
        or column.ktype == knext.int32f,
    )

    def configure(
        self, configure_context, input_schema_1: knext.Schema, input_schema_2
    ):
        return input_schema_1

    def execute(self, exec_context, input_1: knext.Table, input_2: knext.Table):
        if self.molecule_column is None or self.bits_column is None:
            raise AttributeError(
                "Molecule or Bits column was not selected in configuration dialog."
            )

        df = input_1.to_pandas()
        mols = self.get_mols(input_1, df)

        # create list of FP bits from input table 2
        fp_ids = input_2.to_pandas()[self.bits_column]

        cols = [None] * len(fp_ids)
        for i in range(len(cols)):
            cols[i] = []

        for mol in mols:
            bi = {}  # defining a dictionary
            fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(
                mol, radius=self.radius, nBits=self.number_bits, bitInfo=bi
            )  # calculate fingerprint with user-defined radius and nr of bits
            Chem.Kekulize(mol)  # kekulize molecules
            for i, idx in enumerate(
                fp_ids
            ):  # if rendering fails, append an empty cell. Don't make if-else to catch the error
                if fp[idx]:
                    try:
                        cols[i].append(Draw.DrawMorganBit(mol, idx, bi, useSVG=True))
                    except:
                        cols[i].append(None)
                else:
                    cols[i].append(None)

        for i, idx in enumerate(fp_ids):
            df[f"bit{idx}"] = cols[i]

        return knext.Table.from_pandas(df)

    def get_mols(self, input_1, df):
        # Prepare mols: If input column consists of rdkit molecules, use them;
        # if input column consists of SMILES, convert to rdkit molecules
        molecule_column_type = input_1.schema[self.molecule_column_param].ktype
        mols = []
        if molecule_column_type == knext.logical(Chem.rdchem.Mol):
            LOGGER.warning("rdkit mols detected")
            mols = df[self.molecule_column_param]
        else:
            for smi in df[self.molecule_column_param]:
                mols.append(Chem.MolFromSmiles(smi, sanitize=False))
        return mols
