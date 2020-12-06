using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_1d
    {
        internal enum_protein_data_source pds;

        internal bool pse_aac;
        internal bool sable; // sable is providing some kind of 2d or 3d geometric data... but from 1d sequence.
        //internal bool mpsa_classification_data_subsequence = false;
        internal bool blast_pssm;
        internal bool aaindex;
        internal bool geometry;
        internal bool iupred2a; // intrinsically unordered protein
        internal bool stackdppred; // dna binding
        internal bool r_peptides;
        internal bool r_protr;

        public feature_types_1d()
        {

        }

        internal feature_types_1d(enum_protein_data_source pds, bool? enable = null)
        {
            if (pds == enum_protein_data_source.interface_1d) { }
            else if (pds == enum_protein_data_source.neighbourhood_1d) { }
            else if (pds == enum_protein_data_source.chain_1d) { }
            else throw new Exception();

            this.pds = pds;

            if (enable != null) set_enable(enable.Value);
        }

        internal void set_enable(bool enable, string name = null)
        {
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(pse_aac), StringComparison.OrdinalIgnoreCase)) pse_aac = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(sable), StringComparison.OrdinalIgnoreCase)) sable = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(blast_pssm), StringComparison.OrdinalIgnoreCase)) blast_pssm = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(aaindex), StringComparison.OrdinalIgnoreCase)) aaindex = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(geometry), StringComparison.OrdinalIgnoreCase)) geometry = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(iupred2a), StringComparison.OrdinalIgnoreCase)) iupred2a = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(stackdppred), StringComparison.OrdinalIgnoreCase)) stackdppred = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(r_peptides), StringComparison.OrdinalIgnoreCase)) r_peptides = enable;
            if (string.IsNullOrEmpty(name) || string.Equals(name, nameof(r_protr), StringComparison.OrdinalIgnoreCase)) r_protr = enable;
        }

        internal void set_enable((bool enable, string name)[] items)
        {
            foreach (var item in items)
            {
                set_enable(item.enable, item.name);
            }
        }

        internal static readonly string[] keys = new string[] {
            nameof(pse_aac),
            nameof(sable),
            nameof(blast_pssm),
            nameof(aaindex),
            nameof(geometry),
            nameof(iupred2a),
            nameof(stackdppred),
            nameof(r_peptides),
            nameof(r_protr),
            }.OrderBy(a=>a).ToArray();


        internal List<(string key, bool value)> key_value_list()
        {
            var data = new List<(string key, bool value)>()
            {
                (nameof(pse_aac), pse_aac),
                (nameof(sable), sable),
                //(nameof(mpsa_classification_data_subsequence),mpsa_classification_data_subsequence),
                (nameof(blast_pssm), blast_pssm),
                (nameof(aaindex), aaindex),
                (nameof(geometry), geometry),
                (nameof(iupred2a), iupred2a),
                (nameof(stackdppred), stackdppred),
                (nameof(r_peptides), r_peptides),
                (nameof(r_protr), r_protr)
            }.OrderBy(a=>a.key).ThenBy(a => a.value).ToList();

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
