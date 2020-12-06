using System.Collections.Generic;

namespace dimorphics_dataset
{
    internal class info_aaindex_entry
    {
        internal string H_Accession_Number;
        internal string D_Data_Description;
        internal string R_PMID;
        internal string A_Authors;
        internal string T_Title_Of_Article;
        internal string J_Journal_Reference;
        internal List<(string entry_name, double similarity_score)> C_Accession_numbers_of_similar_entries = new List<(string, double)>();
        internal List<(char amino_acid, double index_value)> I_Amino_Acid_Index_Data = new List<(char, double)>();
        internal List<(char amino_acid, double index_value)> I_Amino_Acid_Index_Data_Normalised = new List<(char, double)>();
    }
}
