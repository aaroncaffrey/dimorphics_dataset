using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_1d
    {
        internal bool pse_aac_sequence_classification_data;
        internal bool sable_classification_data; // sable is providing some kind of 2d or 3d geometric data... but from 1d sequence.
        //internal bool mpsa_classification_data_subsequence = false;
        internal bool blast_pssm_subsequence_classification_data;
        internal bool aa_index_classification_data;
        internal bool sequence_geometry_classification_data;
        internal bool intrinsically_unordered_data;
        internal bool dna_binding_prediction_data;
        internal bool r_peptides;
        internal bool r_protr;

        internal feature_types_1d()
        {

        }
        internal feature_types_1d(bool enable)
        {
            pse_aac_sequence_classification_data = enable;
            sable_classification_data = enable;
            //mpsa_classification_data_subsequence = enable;
            blast_pssm_subsequence_classification_data = enable;
            aa_index_classification_data = enable;
            sequence_geometry_classification_data = enable;
            intrinsically_unordered_data = enable;
            dna_binding_prediction_data = enable;
            r_peptides = enable;
            r_protr = enable;
        }

        internal List<(string key, bool value)> AsArray()
        {
            var data = new List<(string key, bool value)>()
                {
                    ( nameof(pse_aac_sequence_classification_data         ),  pse_aac_sequence_classification_data              ) ,
                    ( nameof(sable_classification_data                    ),  sable_classification_data                         ) ,
                    //( nameof(mpsa_classification_data_subsequence         ),  mpsa_classification_data_subsequence              ) ,
                    ( nameof(blast_pssm_subsequence_classification_data   ),  blast_pssm_subsequence_classification_data        ) ,
                    ( nameof(aa_index_classification_data                 ),  aa_index_classification_data                      ) ,
                    ( nameof(sequence_geometry_classification_data        ),  sequence_geometry_classification_data             ) ,
                    ( nameof(intrinsically_unordered_data                 ),  intrinsically_unordered_data                      ) ,
                    ( nameof(dna_binding_prediction_data                  ),  dna_binding_prediction_data                       ) ,
                    ( nameof(r_peptides                                   ),  r_peptides                                        ) ,
                    ( nameof(r_protr                                      ),  r_protr                                           )
                };

            return data;
        }

        internal string key_value_list_hr()
        {
            var data = AsArray();

            var ret = $@"({string.Join(", ", data.Select(a => $"{a.key} = {a.value}").ToList())})";

            return ret;
        }

        public override string ToString()
        {
            return key_value_list_hr();
            //throw new NotImplementedException();
        }
    }
}
