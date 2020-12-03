using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_2d
    {
        internal enum_protein_data_source pds;

        //internal bool sable_classification_data = false;
        internal bool mpsa;


        public feature_types_2d()
        {

        }

        internal feature_types_2d(enum_protein_data_source pds, bool? enable = null)
        {
            if (pds == enum_protein_data_source.interface_2d) { }
            else if (pds == enum_protein_data_source.neighbourhood_2d) { }
            else if (pds == enum_protein_data_source.chain_2d) { }
            else throw new Exception();

            this.pds = pds;

            if (enable != null) set_enable(enable.Value);
        }

        internal void set_enable(bool enable, string name = null)
        {
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(mpsa), StringComparison.InvariantCultureIgnoreCase)) mpsa = enable;
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
            nameof(mpsa)
        }.OrderBy(a => a).ToArray();

        internal List<(string key, bool value)> key_value_list()
        {
            var data = new List<(string key, bool value)>()
            {
                //(nameof(sable_classification_data), sable_classification_data),
                (nameof(mpsa), mpsa),
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
