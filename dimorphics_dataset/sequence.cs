using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace dimorphics_dataset
{
    public class sequence
    {
        public class sequence_id
        {
            public string PdbId;
            public char ChainId;
            public string Mol;
            public string Len;
            public string Description;

            public sequence_id(string sequence_id)
            {
                if (String.IsNullOrWhiteSpace(sequence_id))
                {
                    return;
                }

                if (sequence_id != null && sequence_id.FirstOrDefault() == '>') sequence_id = sequence_id.Substring(1);

                const string molMarker = "mol:";
                const string lenMarker = "length:";

                if (sequence_id.Contains(" " + molMarker) && sequence_id.Contains(" " + lenMarker))
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
                            mol = token.Replace(molMarker, "");
                            continue;
                        }

                        if (String.IsNullOrWhiteSpace(len) && token.Length >= lenMarker.Length && token.Substring(0, lenMarker.Length) == lenMarker)
                        {
                            len = token.Replace(lenMarker, "");
                            continue;
                        }

                    }

                    if (pdbId != null && pdbId.Contains("_")) //(mol == "protein")
                    {
                        chainId = pdbId.Substring(pdbId.IndexOf('_') + 1);
                        pdbId = pdbId.Substring(0, pdbId.IndexOf('_'));
                    }

                    var description = sequence_id.Substring(pdbId.Length + 1 + chainId.Length + 1 + mol.Length + 1 + len.Length + 1);

                    if (description.Length > 0 && description.IndexOf(' ') > -1)
                    {
                        description = description.Substring(description.IndexOf(' ') + 1);
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

            public void Init(string pdbId, char chainId, string mol, string len, string description)
            {
                ChainId = chainId;
                PdbId = pdbId;
                Mol = mol;
                Len = len;
                Description = description;
            }
        }

        private string _id;

        public string Id
        {
            get { return _id; }
            set
            {
                _id = value;
                IdSplit = new sequence_id(this._id);
            }
        }
        public sequence_id IdSplit;
        public string FullSequence;

        public sequence(string id, string fullSequence)
        {
            Id = id;
            FullSequence = fullSequence;
        }

        public static void Save(string outputFastaFile, List<sequence> sequenceList)
        {
            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            io.WriteAllText(outputFastaFile, string.Join("", sequenceList.Select(a => a.Id + "\r\n" + a.FullSequence + "\r\n").ToList()), nameof(sequence), nameof(Save));
        }

        public static void Save(string outputFastaFile, sequence sequence)
        {
            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            io.WriteAllText(outputFastaFile, sequence.Id + "\r\n" + sequence.FullSequence + "\r\n", nameof(sequence), nameof(Save));
        }
    }
}
