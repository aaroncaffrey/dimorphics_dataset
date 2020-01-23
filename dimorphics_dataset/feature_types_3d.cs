using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace dimorphics_dataset
{
    public class feature_types_3d
    {
        internal bool pse_ssc_dssp_classification_data = false;
        //internal bool dssp_dist_classification_data = false;
        internal bool foldx_classification_data = false;
        internal bool ring_classification_data = false;
        internal bool sasa_classification_data = false;
        internal bool tortuosity_classification_data = false;
        internal bool intramolecular_classification_data = false;
        //internal bool atom_distance_classification_data = false;
        //internal bool aa_aa_distances = false;

        public feature_types_3d()
        {

        }
        public feature_types_3d(bool enable)
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

        public List<(string key, bool value)> AsArray()
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

        public override string ToString()
        {
            var data = AsArray();

            var ret = $@"({string.Join(", ", data.Select(a => $"{a.key} = {a.value}").ToList())})";

            return ret;
        }
    }
}
