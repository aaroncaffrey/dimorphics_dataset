using System.Collections.Generic;

namespace dimorphics_dataset
{
    public class protein_subsequence_info
    {
        internal int class_id;
        internal string class_name;
        
        internal string dimer_type;
        internal string parallelism;
        internal string symmetry_mode;
        
        internal string pdb_id;
        internal char chain_id = ' ';
        internal List<(int res_id, char i_code, char amino_acid)> res_ids = new List<(int res_id, char i_code, char amino_acid)>();
        internal string aa_subsequence;
    }
}