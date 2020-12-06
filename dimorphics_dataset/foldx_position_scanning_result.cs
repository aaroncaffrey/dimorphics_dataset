namespace dimorphics_dataset
{
    internal class foldx_position_scanning_result
    {
        // FoldX --command=PositionScan --pdb-dir=c:\betastrands_dataset\pdb\ --pdb=1AJYA.pdb --positions MA30d

        // META30M 0
        // META30G 0.0440333
        // META30A 0.000395327
        // META30L -0.063261

        internal string pdb_id;
        internal char chain_id;
        internal int residue_index;

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
