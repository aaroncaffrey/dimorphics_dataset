using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace dimorphics_dataset
{
    public class feature_types_1d
    {
        internal bool pse_aac_sequence_classification_data = false;
        internal bool sable_classification_data = false;
        internal bool mpsa_classification_data_subsequence = false;
        internal bool blast_pssm_subsequence_classification_data = false;
        internal bool aa_index_classification_data = false;
        internal bool sequence_geometry_classification_data = false;
        internal bool intrinsically_unordered_data = false;
        internal bool dna_binding_prediction_data = false;
        internal bool r_peptides = false;
        internal bool r_protr = false;

        public feature_types_1d()
        {

        }
        public feature_types_1d(bool enable)
        {
            pse_aac_sequence_classification_data = enable;
            sable_classification_data = enable;
            mpsa_classification_data_subsequence = enable;
            blast_pssm_subsequence_classification_data = enable;
            aa_index_classification_data = enable;
            sequence_geometry_classification_data = enable;
            intrinsically_unordered_data = enable;
            dna_binding_prediction_data = enable;
            r_peptides = enable;
            r_protr = enable;
        }

        public List<(string key, bool value)> AsArray()
        {
            var data = new List<(string key, bool value)>()
                {
                    ( nameof(pse_aac_sequence_classification_data         ),  pse_aac_sequence_classification_data              ) ,
                    ( nameof(sable_classification_data                    ),  sable_classification_data                         ) ,
                    ( nameof(mpsa_classification_data_subsequence         ),  mpsa_classification_data_subsequence              ) ,
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

        public override string ToString()
        {
            var data = AsArray();

            var ret = $@"({string.Join(", ", data.Select(a => $"{a.key} = {a.value}").ToList())})";

            return ret;
        }
    }
}
