using System;

namespace dimorphics_dataset
{
    internal class info_blast_pssm_entry
    {
        internal int query_sequence_aa_pos;

        internal char query_sequence_aa;

        //internal int PositionAaPos;
        internal char position_aa;
        internal double score;
        internal int matrix_row_index;
        internal int matrix_column_index;

        internal info_blast_pssm_entry()
        {

        }

        internal info_blast_pssm_entry(info_blast_pssm_entry pssm_entry)
        {
            const string module_name = nameof(pssm_entry);
            const string method_name = nameof(pssm_entry);

            if (pssm_entry == null)
            {
                throw new ArgumentNullException(nameof(pssm_entry));
            }

            this.matrix_column_index = pssm_entry.matrix_column_index;
            this.matrix_row_index = pssm_entry.matrix_row_index;
            this.position_aa = pssm_entry.position_aa;
            this.score = pssm_entry.score;
            this.query_sequence_aa = pssm_entry.query_sequence_aa;
            this.query_sequence_aa_pos = pssm_entry.query_sequence_aa_pos;
        }
    }
}
