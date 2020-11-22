using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class sequence
    {
        internal class sequence_id
        {
            internal string pdb_id;
            internal char chain_id;
            internal string mol;
            internal string len;
            internal string description;

            internal sequence_id(string sequence_id)
            {
                if (String.IsNullOrWhiteSpace(sequence_id))
                {
                    return;
                }

                if (sequence_id != null && sequence_id.FirstOrDefault() == '>') sequence_id = sequence_id.Substring(1);

                const string molMarker = "mol:";
                const string lenMarker = "length:";

                if (sequence_id.Contains($" {molMarker}", StringComparison.InvariantCulture) && sequence_id.Contains($" {lenMarker}", StringComparison.InvariantCulture))
                {
                    var idStrings = sequence_id.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                    var pdbId = idStrings[0][0] == '>' ? idStrings[0].Substring(1) : idStrings[0];
                    var chainId = "";
                    var mol = "";
                    var len = "";

                    foreach (var token in idStrings)
                    {
                        if (!String.IsNullOrWhiteSpace(mol) && !String.IsNullOrWhiteSpace(len))
                        {
                            break;
                        }

                        if (String.IsNullOrWhiteSpace(mol) && token.Length >= molMarker.Length && token.Substring(0, molMarker.Length) == molMarker)
                        {
                            mol = token.Replace(molMarker, "", StringComparison.InvariantCulture);
                            continue;
                        }

                        if (String.IsNullOrWhiteSpace(len) && token.Length >= lenMarker.Length && token.Substring(0, lenMarker.Length) == lenMarker)
                        {
                            len = token.Replace(lenMarker, "", StringComparison.InvariantCulture);
                            continue;
                        }

                    }

                    if (pdbId != null && pdbId.Contains("_", StringComparison.InvariantCulture)) //(mol == "protein")
                    {
                        chainId = pdbId.Substring(pdbId.IndexOf('_', StringComparison.InvariantCulture) + 1);
                        pdbId = pdbId.Substring(0, pdbId.IndexOf('_', StringComparison.InvariantCulture));
                    }

                    var description = sequence_id.Substring(pdbId.Length + 1 + chainId.Length + 1 + mol.Length + 1 + len.Length + 1);

                    if (description.Length > 0 && description.IndexOf(' ', StringComparison.InvariantCulture) > -1)
                    {
                        description = description.Substring(description.IndexOf(' ', StringComparison.InvariantCulture) + 1);
                    }

                    Init(pdbId.ToUpperInvariant(), chainId/*.ToUpperInvariant()*/[0], mol, len, description);
                }
                else
                {

                    var split = sequence_id.Split(new char[] { '_', '-', ':', ';', ',' }, StringSplitOptions.RemoveEmptyEntries);

                    if (split.Length < 2) return;// Init("", ' ', "", "", "");


                    var pdbId = split[0].ToUpperInvariant();
                    var chainId = split[1]/*.ToUpperInvariant()*/;

                    var description = split.Length > 2 ? string.Join(" ", split.Skip(2).ToList()) : "";

                    Init(pdbId, chainId[0], null, null, description);
                }

            }

            internal void Init(string pdb_id, char chain_id, string mol, string len, string description)
            {
                this.pdb_id = pdb_id;
                this.chain_id = chain_id;
                this.mol = mol;
                this.len = len;
                this.description = description;
            }
        }

        private string _id;

        internal string id
        {
            get { return _id; }
            set
            {
                _id = value;
                id_split = new sequence_id(this._id);
            }
        }
        internal sequence_id id_split;
        internal string full_sequence;

        internal sequence(string id, string full_sequence)
        {
            this.id = id;
            this.full_sequence = full_sequence;
        }

        internal static void Save(string outputFastaFile, List<sequence> sequenceList)
        {
            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            io_proxy.WriteAllText(outputFastaFile, string.Join("", sequenceList.Select(a => $"{a.id}\r\n{a.full_sequence}\r\n").ToList()), nameof(sequence), nameof(Save));
        }

        internal static void Save(string outputFastaFile, sequence sequence)
        {
            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            if (sequence == null)
            {
                throw new ArgumentNullException(nameof(sequence));
            }

            io_proxy.WriteAllText(outputFastaFile, $"{sequence.id}\r\n{sequence.full_sequence}\r\n", nameof(sequence), nameof(Save));
        }
    }
}
