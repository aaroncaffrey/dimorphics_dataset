namespace dimorphics_dataset
{
    internal class foldx_ala_scanning_result
    {
        // LYS 86 to ALA energy change is 0.0306484
        // THR 87 to ALA energy change is -0.727468
        // GLU 88 to ALA energy change is -0.923282

        internal string pdb_id;
        internal char chain_id;
        internal int residue_index;

        //internal string res_name;
        //internal string scan_res_name;

        internal char original_foldx_amino_acid_1;
        internal string original_foldx_amino_acid_3;
        internal char original_standard_amino_acid_1;
        internal string original_standard_amino_acid_3;

        internal char mutant_foldx_amino_acid_1;
        internal string mutant_foldx_amino_acid_3;
        internal char mutant_standard_amino_acid_1;
        internal string mutant_standard_amino_acid_3;

        internal double ddg;
    }
}
