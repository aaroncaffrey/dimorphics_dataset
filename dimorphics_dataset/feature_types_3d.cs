using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_3d
    {
        internal bool pse_ssc_dssp_classification_data;
        //internal bool dssp_dist_classification_data = false;
        internal bool foldx_classification_data;
        internal bool ring_classification_data;
        internal bool sasa_classification_data;
        internal bool tortuosity_classification_data;
        internal bool intramolecular_classification_data;
        //internal bool atom_distance_classification_data = false;
        //internal bool aa_aa_distances = false;

        internal feature_types_3d()
        {

        }
        internal feature_types_3d(bool enable)
        {
            pse_ssc_dssp_classification_data = enable;
            //dssp_dist_classification_data = enable;
            foldx_classification_data = enable;
            ring_classification_data = enable;
            sasa_classification_data = enable;
            tortuosity_classification_data = enable;
            intramolecular_classification_data = enable;
            //atom_distance_classification_data = enable;
            //aa_aa_distances = enable;
        }

        internal List<(string key, bool value)> AsArray()
        {
            var data = new List<(string key, bool value)>()
                {
                    ( nameof(pse_ssc_dssp_classification_data         ),  pse_ssc_dssp_classification_data                ) ,
                    ( nameof(foldx_classification_data         ),  foldx_classification_data              ) ,
                    ( nameof(ring_classification_data         ),  ring_classification_data              ) ,
                    ( nameof(sasa_classification_data         ),  sasa_classification_data              ) ,
                    ( nameof(tortuosity_classification_data         ),  tortuosity_classification_data             ) ,
                    ( nameof(intramolecular_classification_data         ),  intramolecular_classification_data              ) ,
                    //( nameof(atom_distance_classification_data         ),  atom_distance_classification_data              ) ,
                    //( nameof(aa_aa_distances         ),  aa_aa_distances               ) ,
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
