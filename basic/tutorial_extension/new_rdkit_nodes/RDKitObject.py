import logging
import knime_extension as knext
import knime_schema as ks
from rdkit import Chem
import pandas as pd
from new_rdkit_nodes.utils import category
LOGGER = logging.getLogger(__name__)

@knext.node(name="RDKitMol from SMILES", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category)
@knext.input_table(name="Input Data", description="Input table containing SMILES")
@knext.output_table(name="Output Data", description="Input table appended with a column containing RDKit molecule")
class RDKitObject:
    """
    This node returns the RDKit molecule in an appended column.
    """

    molecule_column = knext.ColumnParameter(label="Select the column containing SMILES", description="Choose the column from the input table containing the SMILES")

    def configure(self, configure_context, input_schema_1):
        return input_schema_1.append(knext.Column(knext.logical(Chem.Mol), "RDKitMol"))
 
 
    def execute(self, exec_context, input_1):    
        input_1_pandas = input_1.to_pandas() 

        mols = [Chem.MolFromSmiles(smi) for smi in input_1_pandas[self.molecule_column]]

        #LOGGER.warning(mols)
        dtype = knext.logical(Chem.Mol).to_pandas()
        #LOGGER.warning(dtype)
        
        input_1_pandas['RDKitMol'] = pd.Series(mols, dtype=dtype, index=input_1_pandas.index)
        return knext.Table.from_pandas(input_1_pandas)
