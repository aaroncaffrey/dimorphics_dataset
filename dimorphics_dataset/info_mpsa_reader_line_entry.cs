using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class info_mpsa_reader_line_entry
    {
        internal info_mpsa_reader reader;
        internal int index = -1;
        internal char amino_acid = ' ';
        internal char predicted_ss_code = ' ';
        internal List<char> ss_column_headers = new List<char>() { 'H', 'E', 'C', 'T' };
        internal List<(char ss, char amino_acid, double value)> line_prob_values;

        internal info_mpsa_reader_line_entry()
        {

        }

        internal info_mpsa_reader_line_entry(info_mpsa_reader_line_entry mpsa_line_entry)
        {
            if (mpsa_line_entry == null) return;

            this.reader = mpsa_line_entry.reader;
            this.index = mpsa_line_entry.index;
            this.amino_acid = mpsa_line_entry.amino_acid;
            this.predicted_ss_code = mpsa_line_entry.predicted_ss_code;
            this.ss_column_headers = mpsa_line_entry.ss_column_headers.ToList();
            this.line_prob_values = mpsa_line_entry.line_prob_values.ToList();
        }
    }
}
