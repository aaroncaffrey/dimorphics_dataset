using System;
using System.Linq;

namespace dimorphics_dataset
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

            const string molMarker = @"mol:";
            const string lenMarker = @"length:";

            if (sequence_id.Contains(/*program.string_debug*/($@" {molMarker}"), StringComparison.Ordinal) && sequence_id.Contains(/*program.string_debug*/($@" {lenMarker}"), StringComparison.Ordinal))
            {
                var idStrings = sequence_id.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);

                var pdbId = idStrings[0][0] == '>' ? idStrings[0].Substring(1) : idStrings[0];
                var chainId = /*program.string_debug*/($@"");
                var mol = /*program.string_debug*/($@"");
                var len = /*program.string_debug*/($@"");

                foreach (var token in idStrings)
                {
                    if (!String.IsNullOrWhiteSpace(mol) && !String.IsNullOrWhiteSpace(len))
                    {
                        break;
                    }

                    if (String.IsNullOrWhiteSpace(mol) && token.Length >= molMarker.Length && string.Equals(token.Substring(0, molMarker.Length), molMarker, StringComparison.Ordinal))
                    {
                        mol = token.Replace(molMarker, /*program.string_debug*/($@""), StringComparison.Ordinal);
                        continue;
                    }

                    if (String.IsNullOrWhiteSpace(len) && token.Length >= lenMarker.Length && string.Equals(token.Substring(0, lenMarker.Length), lenMarker, StringComparison.Ordinal))
                    {
                        len = token.Replace(lenMarker, /*program.string_debug*/($@""), StringComparison.Ordinal);
                        continue;
                    }

                }

                if (pdbId != null && pdbId.Contains(/*program.string_debug*/($@"_"), StringComparison.Ordinal)) //(mol == /*program.string_debug*/($@"protein")
                {
                    chainId = pdbId.Substring(pdbId.IndexOf('_', StringComparison.Ordinal) + 1);
                    pdbId = pdbId.Substring(0, pdbId.IndexOf('_', StringComparison.Ordinal));
                }

                var description = sequence_id.Substring(pdbId.Length + 1 + chainId.Length + 1 + mol.Length + 1 + len.Length + 1);

                if (description.Length > 0 && description.IndexOf(' ', StringComparison.Ordinal) > -1)
                {
                    description = description.Substring(description.IndexOf(' ', StringComparison.Ordinal) + 1);
                }

                Init(pdbId.ToUpperInvariant(), chainId/*.ToUpperInvariant()*/[0], mol, len, description);
            }
            else
            {

                var split = sequence_id.Split(new char[] { '_', '-', ':', ';', ',' }, StringSplitOptions.RemoveEmptyEntries);

                if (split.Length < 2) return;// Init(/*program.string_debug*/($@""), ' ', /*program.string_debug*/($@""), /*program.string_debug*/($@""), /*program.string_debug*/($@""));


                var pdbId = split[0].ToUpperInvariant();
                var chainId = split[1]/*.ToUpperInvariant()*/;

                var description = split.Length > 2 ? string.Join(/*program.string_debug*/($@" "), split.Skip(2).ToList()) : /*program.string_debug*/($@"");

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
}
