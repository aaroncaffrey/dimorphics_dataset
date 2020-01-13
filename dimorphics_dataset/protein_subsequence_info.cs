using System.Collections.Generic;

namespace dimorphics_dataset
{
    public class protein_subsequence_info
    {
        public int class_id;
        public string class_name;
        
        public string dimer_type;
        public string parallelism;
        public string symmetry_mode;

        public string pdb_id;
        public char chain_id = ' ';
        public List<(int res_id, char i_code, char amino_acid)> res_ids = new List<(int res_id, char i_code, char amino_acid)>();
        public string aa_subsequence;
    }
}