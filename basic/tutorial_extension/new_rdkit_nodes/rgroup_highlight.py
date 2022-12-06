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
Part of the RDKit Python extension. Node 'RGroupHighlightNode'.

@author Daria Goldmann, KNIME GmbH, Konstanz, Germany
@author Steffen Fissler, KNIME GmbH, Konstanz, Germany
"""

from typing import List, Tuple
from pyparsing import col
import pyarrow as pa
import numpy as np
import logging
import knime_extension as knext
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdRGroupDecomposition as rdRGD
from rdkit.Chem import rdDepictor
from rdkit.Chem import rdqueries
from rdkit import Geometry
from rdkit import RDLogger
from new_rdkit_nodes.utils import category
#import chem_extension_types as cet  # To work with and compare against chemical data types like SMILES,...

# Could be removed?
LOGGER = logging.getLogger(__name__)

# Define the node, its inputs, and its outputs below
# Check https://docs.knime.com/latest/pure_python_node_extensions_guide/#_defining_the_nodes_configuration_dialog for further intformation
@knext.node(    name="Highlight R-Groups", node_type=knext.NodeType.MANIPULATOR, icon_path="icon.png", category=category)
# input table with R-groups from R-group decomposition
@knext.input_table(name="R-groups", description="Input Table with R-groups")
# input table with the SMARTS of the core of the molecules
@knext.input_table(name="Core", description="SMARTS of the core of the molecules")
# output table where molecules are oriented as the core and the R-groups are colored
@knext.output_table(
    name="Output Data",
    description="Output Table with a string column containing an svg with colored R-groups",
)
class RGroupHighlightNode(knext.PythonNode):
    """abc

    jklasjdfklsajlkf
    Huiuiuiui
    """

    # select columns with R-groups
    rgroups = knext.MultiColumnParameter(
        "R-Groups",
        "Select the columns with the R-Groups to be highlighted.",
        column_filter=lambda column: column.ktype == knext.logical(Chem.rdchem.Mol),
    )

    def column_is_smiles(column: knext.Column):
        c_type = column.ktype
        smiles_type = knext.logical(cet.SmilesValue)
        smiles_adapter_type = knext.logical(cet.SmilesAdapterValue)
        return c_type == smiles_type or c_type == smiles_adapter_type

    # select a column with the smiles of the molecules
    smiles_col = knext.ColumnParameter(
        "SMILES column",
        "Select the column containing SMILES of the molecules",
        column_filter=column_is_smiles,
        include_row_key=False,
        include_none_column=False,
        port_index=0,
    )

    def column_is_smarts(column: knext.Column):
        c_type = column.ktype
        smarts = knext.logical(cet.SmartsValue)
        smarts_adapter = knext.logical(cet.SmartsAdapterValue)
        return c_type == smarts or c_type == smarts_adapter

    # select a column with the MCS of the Core
    core_col = knext.ColumnParameter(
        "Core column",
        "Select the column containing SMARTS of the core",
        column_filter=column_is_smarts,
        include_row_key=False,
        port_index=1,
    )

    def configure(
        self,
        configure_context,
        input_schema_1: knext.Schema,
        input_schema_2: knext.Schema,
    ):  # the multiple column selection parameter needs to be provided the list of columns of an input table
        if (
            self.rgroups is None
        ):  # needs the check, otherwise the columns are not displayed when opening/executing the node for the first time
            self.rgroups = input_schema_1.column_names

        return input_schema_1.append(knext.Column(ktype=knext.string(), name="rgd_svg"))

    def execute(
        self,
        exec_context: knext.ExecutionContext,
        input_1: knext.Table,
        input_2: knext.Table,
    ):
        if len(self.rgroups) < 1:
            raise AttributeError("At least one R-Group needs to be selected.")
        if self.core_col is None or self.smiles_col is None:
            raise AttributeError(
                "SMILES or SMARTS column was not selected in configuration dialog."
            )

        # adapted from: https://greglandrum.github.io/rdkit-blog/tutorial/prototypes/drawing/2021/08/07/rgd-and-highlighting.html

        r_groups_df = input_1.to_pandas()
        smarts_core_df = input_2.to_pandas()

        # code for the r-group decomp highlight
        subImageSize = (250, 200)

        rdDepictor.SetPreferCoordGen(True)

        # column with the core
        # Example from the blog post
        core_smarts = smarts_core_df[self.core_col][0]

        # get labels of the colums with the R-groups
        lbls = tuple(col_name for col_name in self.rgroups)
        core = Chem.MolFromSmarts(core_smarts)
        rdDepictor.Compute2DCoords(core)

        ps = Chem.AdjustQueryParameters.NoAdjustments()
        ps.makeDummiesQueries = True
        qcore = Chem.AdjustQueryProperties(core, ps)

        ms = []
        progress = 0.0
        add_to_progress = 1 / (input_1.num_rows * 10)
        for smi in r_groups_df[self.smiles_col]:
            ms.append(Chem.MolFromSmiles(smi))
            progress = progress + add_to_progress
            exec_context.set_progress(progress=progress)

        mhs = [Chem.AddHs(x, addCoords=True) for x in ms]
        mms = [x for x in mhs if x.HasSubstructMatch(qcore)]

        for m in mms:
            for atom in m.GetAtoms():
                atom.SetIntProp("SourceAtomIdx", atom.GetIdx())
            progress = progress + add_to_progress
            exec_context.set_progress(progress=progress)

        RDLogger.DisableLog("rdApp.warning")
        groups, _ = rdRGD.RGroupDecompose([qcore], mms, asSmiles=False, asRows=True)
        svgs = []

        for i, m in enumerate(mhs):
            svg = self.highlight_rgroups(
                m,
                groups[i],
                qcore,
                lbls=lbls,
                width=subImageSize[0],
                height=subImageSize[1],
            )
            svgs.append(svg)
            progress = progress + add_to_progress * 8
            exec_context.set_progress(progress=progress)

        r_groups_df["rgd_svg"] = svgs
        return knext.Table.from_pandas(r_groups_df)

    def highlight_rgroups(
        self,
        mol,
        row,
        core,
        width=350,
        height=200,
        fillRings=True,
        legend="",
        sourceIdxProperty="SourceAtomIdx",
        lbls=("R1", "R2", "R3", "R4"),
    ):
        from collections import defaultdict

        # copy the molecule and core
        mol = Chem.Mol(mol)
        core = Chem.Mol(core)

        # -------------------------------------------
        # include the atom map numbers in the substructure search in order to
        # try to ensure a good alignment of the molecule to symmetric cores
        for at in core.GetAtoms():
            if at.GetAtomMapNum():
                at.ExpandQuery(
                    rdqueries.IsotopeEqualsQueryAtom(200 + at.GetAtomMapNum())
                )

        for lbl in row:
            if lbl == "Core":
                continue
            rg = row[lbl]
            for at in rg.GetAtoms():
                if (
                    not at.GetAtomicNum()
                    and at.GetAtomMapNum()
                    and at.HasProp("dummyLabel")
                    and at.GetProp("dummyLabel") == lbl
                ):
                    # attachment point. the atoms connected to this
                    # should be from the molecule
                    for nbr in at.GetNeighbors():
                        if nbr.HasProp(sourceIdxProperty):
                            mAt = mol.GetAtomWithIdx(nbr.GetIntProp(sourceIdxProperty))
                            if mAt.GetIsotope():
                                mAt.SetIntProp("_OrigIsotope", mAt.GetIsotope())
                            mAt.SetIsotope(200 + at.GetAtomMapNum())
        # remove unmapped hs so that they don't mess up the depiction
        rhps = Chem.RemoveHsParameters()
        rhps.removeMapped = False
        tmol = Chem.RemoveHs(mol, rhps)
        rdDepictor.GenerateDepictionMatching2DStructure(tmol, core)

        oldNewAtomMap = {}
        # reset the original isotope values and account for the fact that
        # removing the Hs changed atom indices
        for i, at in enumerate(tmol.GetAtoms()):
            if at.HasProp(sourceIdxProperty):
                oldNewAtomMap[at.GetIntProp(sourceIdxProperty)] = i
                if at.HasProp("_OrigIsotope"):
                    at.SetIsotope(at.GetIntProp("_OrigIsotope"))
                    at.ClearProp("_OrigIsotope")
                else:
                    at.SetIsotope(0)

        # ------------------
        #  set up our colormap
        #   the three choices here are all "colorblind" colormaps

        # "Tol" colormap from https://davidmathlogic.com/colorblind
        colors = [
            (51, 34, 136),
            (17, 119, 51),
            (68, 170, 153),
            (136, 204, 238),
            (221, 204, 119),
            (204, 102, 119),
            (170, 68, 153),
            (136, 34, 85),
        ]
        # "IBM" colormap from https://davidmathlogic.com/colorblind
        colors = [
            (100, 143, 255),
            (120, 94, 240),
            (220, 38, 127),
            (254, 97, 0),
            (255, 176, 0),
        ]
        # Okabe_Ito colormap from https://jfly.uni-koeln.de/color/
        colors = [
            (230, 159, 0),
            (86, 180, 233),
            (0, 158, 115),
            (240, 228, 66),
            (0, 114, 178),
            (213, 94, 0),
            (204, 121, 167),
        ]
        for i, x in enumerate(colors):
            colors[i] = tuple(y / 255 for y in x)

        # ----------------------
        # Identify and store which atoms, bonds, and rings we'll be highlighting
        highlightatoms = defaultdict(list)
        highlightbonds = defaultdict(list)
        atomrads = {}
        widthmults = {}

        rings = []
        for i, lbl in enumerate(lbls):
            color = colors[i % len(colors)]
            rquery = row[lbl]
            Chem.GetSSSR(rquery)
            rinfo = rquery.GetRingInfo()
            for at in rquery.GetAtoms():
                if at.HasProp(sourceIdxProperty):
                    origIdx = oldNewAtomMap[at.GetIntProp(sourceIdxProperty)]
                    highlightatoms[origIdx].append(color)
                    atomrads[origIdx] = 0.4
            if fillRings:
                for aring in rinfo.AtomRings():
                    tring = []
                    allFound = True
                    for aid in aring:
                        at = rquery.GetAtomWithIdx(aid)
                        if not at.HasProp(sourceIdxProperty):
                            allFound = False
                            break
                        tring.append(oldNewAtomMap[at.GetIntProp(sourceIdxProperty)])
                    if allFound:
                        rings.append((tring, color))
            for qbnd in rquery.GetBonds():
                batom = qbnd.GetBeginAtom()
                eatom = qbnd.GetEndAtom()
                if batom.HasProp(sourceIdxProperty) and eatom.HasProp(
                    sourceIdxProperty
                ):
                    origBnd = tmol.GetBondBetweenAtoms(
                        oldNewAtomMap[batom.GetIntProp(sourceIdxProperty)],
                        oldNewAtomMap[eatom.GetIntProp(sourceIdxProperty)],
                    )
                    bndIdx = origBnd.GetIdx()
                    highlightbonds[bndIdx].append(color)
                    widthmults[bndIdx] = 2

        d2d = Draw.MolDraw2DSVG(width, height)
        dos = d2d.drawOptions()
        dos.useBWAtomPalette()

        # ----------------------
        # if we are filling rings, go ahead and do that first so that we draw
        # the molecule on top of the filled rings
        if fillRings and rings:
            # a hack to set the molecule scale
            d2d.DrawMoleculeWithHighlights(
                tmol,
                legend,
                dict(highlightatoms),
                dict(highlightbonds),
                atomrads,
                widthmults,
            )
            d2d.ClearDrawing()
            conf = tmol.GetConformer()
            for (aring, color) in rings:
                ps = []
                for aidx in aring:
                    pos = Geometry.Point2D(conf.GetAtomPosition(aidx))
                    ps.append(pos)
                d2d.SetFillPolys(True)
                d2d.SetColour(color)
                d2d.DrawPolygon(ps)
            dos.clearBackground = False

        # ----------------------
        # now draw the molecule, with highlights:
        d2d.DrawMoleculeWithHighlights(
            tmol,
            legend,
            dict(highlightatoms),
            dict(highlightbonds),
            atomrads,
            widthmults,
        )
        d2d.FinishDrawing()
        svg = d2d.GetDrawingText()
        return svg
