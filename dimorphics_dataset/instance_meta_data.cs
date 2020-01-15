using System;
using System.Collections.Generic;
using System.Text;

namespace dimorphics_dataset
{
    public class instance_meta_data
    {
        internal string protein_pdb_id;
        internal char protein_chain_id;
        internal string protein_dimer_type;
        
        internal int subsequence_class_id;
        internal string subsequence_class_name;
        
        internal string subsequence_parallelism;
        internal string subsequence_symmetry_mode;
        
        internal string subsequence_aa_seq;
        internal List<(int residue_index, char i_code, char amino_acid)> subsequence_res_ids;
        
        internal string subsequence_dssp_monomer;
        internal string subsequence_dssp3_monomer;
        internal string subsequence_dssp_multimer;
        internal string subsequence_dssp3_multimer;
        internal List<(string format, string prediction)> ss_predictions;
        
        internal string protein_aa_seq;
    }
}
