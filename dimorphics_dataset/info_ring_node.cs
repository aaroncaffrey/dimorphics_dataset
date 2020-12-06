using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace dimorphics_dataset
{
    internal class info_ring_node
    {
        //http://protein.bio.unipd.it/ring/help
        //"NodeId       Chain   Position   Residue   Dssp    Degree   Bfactor_CA    x         y         z         pdbFileName              Rapdf     Tap"
        // A:31:_:SER   A       31         SER       (tab)   1        0.000         12.051    -15.082   -33.581   1AJYA_Repair.pdb#31.A    -26.022   0.315

        internal (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId; //The node ID.See the corresponding edge attributes for more details about the format.
        internal char Chain;//Three columns that reports the residue (or ligand) chain, position (PDB index) and residue 3 letters code.
        internal int Position;
        internal char Residue1;
        internal string Residue3;
        internal string Dssp;//The secondary structure calculated with the DSSP algorithm [5] re-implemented in RING-2.0. 
        internal double Degree;//The node degree, i.e. the number of directly connected nodes
        internal double Bfactor_CA;//The B factor of the C-alpha (when available).
        internal double x;//The 3D coordinates of the C-alpha as reported in the original PDB.
        internal double y;
        internal double z;
        internal string pdbFileName;//<pdb_file_name>#<index>.<chain> This is necessary to bind RINanylezer/StructureViz with Chimera visualization of the structure.
        internal double Rapdf;//(Optional). The RAPDF energy. It is calculated based on statistical potentials, see [2]. Tosatto, S.C.E., 2005. The victor/FRST function for model quality estimation. J. Comput. Biol. 12, 1316–1327. doi:10.1089/cmb.2005.12.1316
        internal double Tap;//(Optional). The TAP energy. It is calculated based on statistical potentials, see [3]. Tosatto, S.C.E., Battistutta, R., 2007. TAP score: torsion angle propensity normalization applied to local protein structure evaluation. BMC Bioinformatics 8, 155. doi:10.1186/1471-2105-8-155

        internal static List<info_ring_node> load(string filename)
        {
            var result = new List<info_ring_node>();

            var data = io_proxy.ReadAllLines(filename, nameof(info_ring_node), nameof(load)).Skip(1).Where(a => !string.IsNullOrWhiteSpace(a)).Select(a => a.Split('\t').Select(b => b.Trim()).Select(b => string.Equals(b, $@"-999.9", StringComparison.Ordinal) ? $@"0" : b).ToList()).ToList();

            foreach (var d in data)
            {
                var NodeId = (text: d[0], chain: d[0].Split(':')[0][0], res_id: int.Parse(d[0].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo), icode: d[0].Split(':')[2][0] == '_' ? ' ' : d[0].Split(':')[2][0],
                    amino_acid1: atom.Aa3To1(d[0].Split(':')[3]), amino_acid3: d[0].Split(':')[3]);

                result.Add(new info_ring_node()
                {
                    NodeId = NodeId,
                    Chain = d[1][0],
                    Position = d[2].Length > 0 ? int.Parse(d[2], NumberStyles.Integer, NumberFormatInfo.InvariantInfo) : 0,
                    Residue1 = atom.Aa3To1(d[3]),
                    Residue3 = d[3],
                    Dssp = d[4],
                    Degree = d[5].Length > 0 ? double.Parse(d[5], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Bfactor_CA = d[6].Length > 0 ? double.Parse(d[6], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    x = d[7].Length > 0 ? double.Parse(d[7], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    y = d[8].Length > 0 ? double.Parse(d[8], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    z = d[9].Length > 0 ? double.Parse(d[9], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    pdbFileName = d[10],
                    Rapdf = d[11].Length > 0 ? double.Parse(d[11], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Tap = d[12].Length > 0 ? double.Parse(d[12], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                });
            }

            return result;
        }
    }
}
