import logging
import knime_extension as knext
from rdkit import Chem
import rdkit.Chem.rdMolDescriptors
import pandas as pd
from new_rdkit_nodes.utils import category
LOGGER = logging.getLogger(__name__)

@knext.node(name="Get Molecular Formula", node_type=knext.NodeType.LEARNER, icon_path="new_rdkit_nodes/chem.png", category=category)
@knext.input_table(name="Input Data", description="Input table containing SMILES")
@knext.output_table(name="Output Data", description="Input table appended with a column containing the molecular formula")
class MolFormula(knext.PythonNode):
    """
    This node returns the molecular formula from SMILES in an appended column.
    """

    molecule_column = knext.ColumnParameter(label="Select the column containing SMILES", description="Choose the column from the input table containing the SMILES")

    def configure(self, configure_context, input_schema_1):
        return input_schema_1.append(knext.Column(knext.string(), "MolFormula"))
 
 
    def execute(self, exec_context, input_1):    
        input_1_pandas = input_1.to_pandas() 

        mols = input_1_pandas[self.molecule_column]
        #mols = [Chem.MolFromSmiles(smi) for smi in input_1_pandas[self.molecule_column]]

        MolFormula = []
        for x in mols: 
            if x is None: 
                MolFormula.append(None)
            else: 
                MolFormula.append(Chem.rdMolDescriptors.CalcMolFormula(x))

        input_1_pandas['MolFormula'] = MolFormula
        return knext.Table.from_pandas(input_1_pandas)
