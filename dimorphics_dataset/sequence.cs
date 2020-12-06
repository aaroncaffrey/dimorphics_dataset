using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class sequence
    {
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

            io_proxy.WriteAllText(outputFastaFile, string.Join($@"", sequenceList.Select(a => $"{a.id}\r\n{a.full_sequence}\r\n").ToList()), nameof(sequence), nameof(Save));
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
