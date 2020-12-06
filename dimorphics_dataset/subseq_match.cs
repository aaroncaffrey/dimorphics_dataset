namespace dimorphics_dataset
{
    internal class subseq_match
    {
        internal string pdb_id;
        internal int chain_number;
        internal int alphabetical_chain_number;
        //internal List<Atom> pdb_model_chain_atoms;
        internal char chain_id;
        internal int pdb_model_index;
        internal string seq;
        internal string subsequence;
        internal string dssp;
        internal string stride;
        internal int total_matches;
        internal int array_index;
        internal int residue_index;

        //internal List<(int array_index, int residue_index)> subseq_indexes;
    }
}
