using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;

namespace dimorphics_dataset
{
    public class protein_subsequence_info
    {
        internal int class_id;
        internal string class_name;
        
        internal string dimer_type;
        internal string parallelism;
        internal string symmetry_mode;
        
        internal string pdb_id;
        internal char chain_id = ' ';
        internal List<(int res_id, char i_code, char amino_acid)> res_ids = new List<(int res_id, char i_code, char amino_acid)>();
        internal string aa_subsequence;
        //internal string aa_before;
        //internal string aa_after;
        //internal string aa_protein;

        public static void save(string filename, List<protein_subsequence_info> psi_list)
        {
            var text = psi_list.Select(a => $@"{a.class_id},{a.class_name},{a.dimer_type},{a.parallelism},{a.symmetry_mode},{a.pdb_id},{a.chain_id},{string.Join(";", a.res_ids.Select(b => $@"{b.amino_acid}{b.i_code}{b.res_id}").ToList())},{a.aa_subsequence}").ToList();
            text.Insert(0, $@"{nameof(class_id)},{nameof(class_name)},{nameof(dimer_type)},{nameof(parallelism)},{nameof(symmetry_mode)},{nameof(pdb_id)},{nameof(chain_id)},{nameof(res_ids)},{nameof(aa_subsequence)}");
            io_proxy.WriteAllLines(filename, text, nameof(protein_subsequence_info), nameof(save));
        }

        public static List<protein_subsequence_info> load(string filename)
        {
            if (!File.Exists(filename) || new FileInfo(filename).Length == 0) return null;

            var text = io_proxy.ReadAllLines(filename,nameof(protein_subsequence_info), nameof(load));

            var data = text.Skip(1).Select(a =>
            {
                var b = a.Split(',');
                var j = 0;
                return new protein_subsequence_info()
                {
                    class_id = int.Parse(b[j++], NumberStyles.Integer, CultureInfo.InvariantCulture),
                    class_name = b[j++],
                    dimer_type = b[j++],
                    parallelism = b[j++],
                    symmetry_mode = b[j++],
                    pdb_id = b[j++],
                    chain_id = b[j++][0],
                    res_ids = b[j++].Split(';').Select(c => (res_id: int.Parse(c.Substring(2), NumberStyles.Integer, CultureInfo.InvariantCulture), i_code: c[1], amino_acid: c[0])).ToList(),
                    aa_subsequence = b[j++],
                    //aa_before = b[j++],
                    //aa_after = b[j++],
                    //aa_protein = b[j++],
                };
            }).ToList();

            return data;
        }
    }
}