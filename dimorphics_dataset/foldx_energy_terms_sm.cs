﻿using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal class foldx_energy_terms_sm : foldx_energy_terms
    {
        internal List<(char original_amino_acid1, char chain_id, int residue_index, char mutant_foldx_amino_acid1, string mutant_foldx_amino_acid3, char mutant_standard_amino_acid1, string mutant_standard_amino_acid3)> mutation_positions_data;
    }
}
