using System;
using System.Collections.Generic;
using System.Text;

namespace dimorphics_dataset
{
    public class instance_meta_data
    {
        public string protein_pdb_id;
        public char protein_chain_id;
        public string protein_dimer_type;

        public int subsequence_class_id;
        public string subsequence_class_name;

        public string subsequence_parallelism;
        public string subsequence_symmetry_mode;

        public string subsequence_aa_seq;
        public List<(int residue_index, char i_code, char amino_acid)> subsequence_res_ids;

        public string subsequence_dssp_monomer;
        public string subsequence_dssp3_monomer;
        public string subsequence_dssp_multimer;
        public string subsequence_dssp3_multimer;
        public List<(string format, string prediction)> ss_predictions;

        public string protein_aa_seq;
    }
}
