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
            const string module_name = nameof(sequence);
            const string method_name = nameof(Save);

            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            io_proxy.WriteAllText(outputFastaFile, string.Join(/*program.string_debug*/($@""), sequenceList.Select(a => /*program.string_debug*/($"{a.id}\r\n{a.full_sequence}\r\n")).ToList()), module_name, method_name);
        }

        internal static void Save(string outputFastaFile, sequence sequence)
        {
            const string module_name = nameof(sequence);
            const string method_name = nameof(Save);

            if (String.IsNullOrWhiteSpace(outputFastaFile)) return;

            if (sequence == null)
            {
                throw new ArgumentNullException(nameof(sequence));
            }

            io_proxy.WriteAllText(outputFastaFile, /*program.string_debug*/($"{sequence.id}\r\n{sequence.full_sequence}\r\n"), module_name, method_name);
        }
    }
}
