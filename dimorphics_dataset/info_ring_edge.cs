using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace dimorphics_dataset
{
    internal class info_ring_edge
    {
        //http://protein.bio.unipd.it/ring/help
        //"NodeId1     Interaction   NodeId2      Distance   Angle    Energy   Atom1   Atom2   Donor        Positive	Cation	Orientation"
        //A:31:_:SER   HBOND:SC_SC   A:43:_:LYS   2.911      12.333   17.000   OG      NZ      A:43:_:LYS   		

        internal (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId1; // These two columns report the source and the target node. Nodes can be either a residue or a ligand molecule. The node ID format follows the RIN Analyzer and RING (first version) standards.  <chain> : <index> : <insertion_code> : <residue_3_letter_code> 
        internal (string text, string interaction_type, string subtype_node1, string subtype_node2) Interaction; // This attribute reports the type of interaction and the interaction subtype. <interaction_type> : <subtype_node1> _ <subtype_node2> Where subtypes values are: main chain (MC), side chain (SC) and ligand (LIG).
        internal (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId2;
        internal double Distance;//The distance in Å between atom centers / mass centers / barycenters depending on the type of interaction and the type of residue.
        internal double? Angle;//The angle in degree. Depending on the type of interaction it is calculated differently as described in the previous section. For VDW interactions and IAC it is not possible to infer any angle and the -999.9 (NULL) value is assigned.

        internal double Energy;//The average bond free energy in Kj/mol according to literature. For hydrogen bonds it varies according to the involved donor/acceptor atoms. For disulphide bonds dissociation enthalpy is reported. For VDW the enrgy corresponds to the average attractiveness component of the van der Waals force. All energy values reported by RING-2.0 are listed in the following table.
        internal string Atom1;
        internal string Atom2;
        internal (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) Donor;
        internal string Positive;
        internal string Cation;
        internal string Orientation;


        internal static List<info_ring_edge> load(string filename)
        {
            var result = new List<info_ring_edge>();

            var data = io_proxy.ReadAllLines(filename, nameof(info_ring_edge), nameof(load)).Skip(1).Where(a => !string.IsNullOrWhiteSpace(a)).Select(a => a.Split('\t').Select(b => b.Trim()).Select(b => string.Equals(b, /*program.string_debug*/($@"-999.9"), StringComparison.Ordinal) ? /*program.string_debug*/($@"0") : b).ToList()).ToList();

            foreach (var d in data)
            {

                var NodeId1_a = (text: d[0], chain: d[0].Split(':')[0][0], res_id: int.Parse(d[0].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo), icode: d[0].Split(':')[2][0] == '_' ? ' ' : d[0].Split(':')[2][0], amino_acid1: atom.Aa3To1(d[0].Split(':')[3]), amino_acid3: d[0].Split(':')[3]);
                var Interaction_a = (text: d[1], interaction_type: d[1].Split(':')[0], subtype_node1: d[1].Split(':')[1].Split('_')[0], subtype_node2: d[1].Split(':')[1].Split('_')[1]);

                var NodeId2_a = (text: d[2], chain: d[2].Split(':')[0][0], res_id: int.Parse(d[2].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo), icode: d[2].Split(':')[2][0] == '_' ? ' ' : d[2].Split(':')[2][0], amino_acid1: atom.Aa3To1(d[2].Split(':')[3]), amino_acid3: d[2].Split(':')[3]);
                var Donor_a = (text: d[8], chain: d[8].Length > 0 ? d[8].Split(':')[0][0] : ' ', res_id: d[8].Length > 0 ? int.Parse(d[8].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo) : 0,

                    icode: d[8].Length > 0 ? (d[8].Split(':')[2][0] == '_' ? ' ' : d[8].Split(':')[2][0]) : ' ', amino_acid1: d[8].Length > 0 ? atom.Aa3To1(d[8].Split(':')[3]) : ' ', amino_acid3: d[8].Length > 0 ? d[8].Split(':')[3] : /*program.string_debug*/($@""));

                var NodeId2_b = (text: d[0], chain: d[0].Split(':')[0][0], res_id: int.Parse(d[0].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo), icode: d[0].Split(':')[2][0] == '_' ? ' ' : d[0].Split(':')[2][0], amino_acid1: atom.Aa3To1(d[0].Split(':')[3]), amino_acid3: d[0].Split(':')[3]);
                var Interaction_b = (text: d[1], interaction_type: d[1].Split(':')[0], subtype_node1: d[1].Split(':')[1].Split('_')[1], subtype_node2: d[1].Split(':')[1].Split('_')[0]);

                var NodeId1_b = (text: d[2], chain: d[2].Split(':')[0][0], res_id: int.Parse(d[2].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo), icode: d[2].Split(':')[2][0] == '_' ? ' ' : d[2].Split(':')[2][0], amino_acid1: atom.Aa3To1(d[2].Split(':')[3]), amino_acid3: d[2].Split(':')[3]);
                var Donor_b = (text: d[8], chain: d[8].Length > 0 ? d[8].Split(':')[0][0] : ' ', res_id: d[8].Length > 0 ? int.Parse(d[8].Split(':')[1], NumberStyles.Integer, NumberFormatInfo.InvariantInfo) : 0,
                    icode: d[8].Length > 0 ? (d[8].Split(':')[2][0] == '_' ? ' ' : d[8].Split(':')[2][0]) : ' ', amino_acid1: d[8].Length > 0 ? atom.Aa3To1(d[8].Split(':')[3]) : ' ', amino_acid3: d[8].Length > 0 ? d[8].Split(':')[3] : /*program.string_debug*/($@""));

                result.Add(new info_ring_edge()
                {
                    NodeId1 = NodeId1_a,
                    Interaction = Interaction_a,
                    NodeId2 = NodeId2_a,
                    Distance = d[3].Length > 0 ? double.Parse(d[3], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Angle = d[4].Length > 0 && !string.Equals(d[4], /*program.string_debug*/($@"-999.9"), StringComparison.Ordinal) ? (double?)double.Parse(d[4], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : null,
                    Energy = d[5].Length > 0 ? double.Parse(d[5], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Atom1 = d[6],
                    Atom2 = d[7],
                    Donor = Donor_a,
                    Positive = d[9],
                    Cation = d[10],
                    Orientation = d[11],
                });

                result.Add(new info_ring_edge()
                {
                    NodeId1 = NodeId1_b,
                    Interaction = Interaction_b,
                    NodeId2 = NodeId2_b,
                    Distance = d[3].Length > 0 ? double.Parse(d[3], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Angle = d[4].Length > 0 && !string.Equals(d[4], /*program.string_debug*/($@"-999.9"), StringComparison.Ordinal) ? (double?)double.Parse(d[4], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : null,
                    Energy = d[5].Length > 0 ? double.Parse(d[5], NumberStyles.Float, NumberFormatInfo.InvariantInfo) : 0,
                    Atom1 = d[6],
                    Atom2 = d[7],
                    Donor = Donor_b,
                    Positive = d[9],
                    Cation = d[10],
                    Orientation = d[11],
                });

                //io_proxy.WriteLine();
            }

            return result;
        }
    }
}
