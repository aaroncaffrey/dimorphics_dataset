using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_3d
    {
        internal enum_protein_data_source pds;

        internal bool pse_ssc_dssp;
        //internal bool dssp_dist_classification_data = false;
        internal bool foldx;
        internal bool ring;
        internal bool sasa;
        internal bool tortuosity;
        internal bool intramolecular;
        //internal bool atom_distance_classification_data = false;
        //internal bool aa_aa_distances = false;

        public feature_types_3d()
        {
            
        }

        internal feature_types_3d(enum_protein_data_source pds, bool? enable = null)
        {
            if (pds == enum_protein_data_source.interface_3d) { }
            else if (pds == enum_protein_data_source.neighbourhood_3d) { }
            else if (pds == enum_protein_data_source.chain_3d) { }
            else throw new Exception();

            this.pds = pds;

            if (enable != null) set_enable(enable.Value);
        }

        internal void set_enable(bool enable, string name = null)
        {
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(pse_ssc_dssp), StringComparison.OrdinalIgnoreCase)) pse_ssc_dssp = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(foldx), StringComparison.OrdinalIgnoreCase)) foldx = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(ring), StringComparison.OrdinalIgnoreCase)) ring = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(sasa), StringComparison.OrdinalIgnoreCase)) sasa = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(tortuosity), StringComparison.OrdinalIgnoreCase)) tortuosity = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(intramolecular), StringComparison.OrdinalIgnoreCase)) intramolecular = enable;

        }

        internal void set_enable((bool enable, string name)[] items)
        {
            foreach (var item in items)
            {
                set_enable(item.enable, item.name);
            }
        }

        internal static readonly string[] keys = new string[]
            {
               nameof(pse_ssc_dssp    ),
               nameof(foldx           ),
               nameof(ring            ),
               nameof(sasa            ),
               nameof(tortuosity      ),
               nameof(intramolecular  )
            }.OrderBy(a => a).ToArray();

        internal List<(string key, bool value)> key_value_list()
        {
            var data = new List<(string key, bool value)>()
            {
                ( nameof(pse_ssc_dssp         ),  pse_ssc_dssp                ) ,
                ( nameof(foldx         ),  foldx              ) ,
                ( nameof(ring         ),  ring              ) ,
                ( nameof(sasa         ),  sasa              ) ,
                ( nameof(tortuosity         ),  tortuosity             ) ,
                ( nameof(intramolecular         ),  intramolecular              ) ,
                //( nameof(atom_distance_classification_data         ),  atom_distance_classification_data              ) ,
                //( nameof(aa_aa_distances         ),  aa_aa_distances               ) ,
            }.OrderBy(a => a.key).ThenBy(a => a.value).ToList();

            return data;
        }

        internal string key_value_list_hr()
        {
            var data = key_value_list();

            var ret = $@"({string.Join($@", ", data.Select(a => $@"{a.key} = {a.value}").ToList())})";

            return ret;
        }

        public override string ToString()
        {
            return key_value_list_hr();
            //throw new NotImplementedException();
        }
    }
}
