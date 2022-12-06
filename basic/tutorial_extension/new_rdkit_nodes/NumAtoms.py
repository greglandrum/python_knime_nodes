import logging
import knime_extension as knext
from rdkit import Chem
import pandas as pd
import numpy as np
from new_rdkit_nodes.utils import category
LOGGER = logging.getLogger(__name__)


@knext.node(name="Get Number of Atoms", node_type=knext.NodeType.LEARNER, icon_path="icon.png", category=category)
@knext.input_table(name="Input Data", description="Input table containing SMILES")
@knext.output_table(name="Output Data", description="Input table appended with a column containing the number of atoms")
class AtomNode:
    """
    This node calculates the number of atoms from SMILES and appends a column.
    """

    molecule_column = knext.ColumnParameter(label="Select the column containing SMILES", description="Choose the column from the input table containing the SMILES")

    def configure(self, configure_context, input_schema_1):
        return input_schema_1.append(knext.Column(knext.int32(), "NumAtoms"))
 
 
    def execute(self, exec_context, input_1):    
        input_1_pandas = input_1.to_pandas() 

        LOGGER.warning(input_1_pandas[self.molecule_column].dtype)

        mols = []
        for smi in input_1_pandas[self.molecule_column]:
            mols.append(Chem.MolFromSmiles(smi))

        NumAtoms = []
        for mol in mols:
            NumAtoms.append(mol.GetNumAtoms())

        input_1_pandas['NumAtoms'] = pd.Series(NumAtoms, dtype=np.int32, index=input_1_pandas.index)
        return knext.Table.from_pandas(input_1_pandas)
