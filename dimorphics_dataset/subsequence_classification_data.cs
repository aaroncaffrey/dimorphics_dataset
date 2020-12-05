using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;

namespace dimorphics_dataset
{
   
    internal class subsequence_classification_data
    {
        internal int class_id;
        internal string class_name;
        internal string dimer_type;
        internal string parallelism;
        internal string symmetry_mode;
        internal string pdb_id;
        internal char chain_id;

        internal subsequence_classification_data_region interface_region;
        internal subsequence_classification_data_region nh_flank_region;
        internal subsequence_classification_data_region nh_contact_region;
        internal subsequence_classification_data_region chain_region;
        
        internal subsequence_classification_data()
        {

        }

        internal (string region_name, subsequence_classification_data_region region)[] get_regions()
        {
            return new (string region_name, subsequence_classification_data_region region)[]
            {
                (nameof(interface_region), interface_region),
                (nameof(nh_flank_region), nh_flank_region),
                (nameof(nh_contact_region), nh_contact_region),
                (nameof(chain_region), chain_region),
            };
        }

        internal void init_nh_flanking(int neighbourhood_flanking_size = 6, bool load_ss_predictions = true, bool load_foldx_energy = true)
        {
            const string module_name = nameof(subsequence_classification_data);
            const string method_name = nameof(init_nh_flanking);

            if (neighbourhood_flanking_size % 2 != 0 || neighbourhood_flanking_size % 3 != 0)
            {
                throw new Exception($@"{module_name}.{method_name}: Number must be divisible by both 2 and 3");
            }
            

            var interface_atoms_array_indexes = interface_region.master_atoms.Select(a => chain_region.master_atoms.IndexOf(a)).Where(a => a != -1).ToList();
            var start_array_index = interface_atoms_array_indexes.Min(); // first res_id in interface or subseq
            var end_array_index = interface_atoms_array_indexes.Max(); // last res_id in interface or subseq

            var all_array_indexes = interface_region.master_atoms.Select(a => a.residue_index).ToList();
            var total_indexes_before_start = all_array_indexes.Count(a => a < start_array_index);
            var total_indexes_after_end = all_array_indexes.Count(a => a > end_array_index);

            var neighbourhood_size_before_start = neighbourhood_flanking_size;
            var neighbourhood_size_after_end = neighbourhood_flanking_size;

            if (total_indexes_before_start < neighbourhood_flanking_size && total_indexes_after_end > neighbourhood_flanking_size)
            {
                neighbourhood_size_after_end += (neighbourhood_flanking_size - total_indexes_before_start);
            }

            if (total_indexes_after_end < neighbourhood_flanking_size && total_indexes_before_start > neighbourhood_flanking_size)
            {
                neighbourhood_size_before_start += (neighbourhood_flanking_size - total_indexes_after_end);
            }

            var neighbourhood_flanking_master_atoms = chain_region.master_atoms.Where((a, i) => (i >= start_array_index - neighbourhood_size_before_start && i < start_array_index) || (i > end_array_index && i <= end_array_index + neighbourhood_size_after_end)).Except(interface_region.atoms).ToList();
            var neighbourhood_flanking_atoms = neighbourhood_flanking_master_atoms.SelectMany(a => a.amino_acid_atoms).Distinct().ToList();
            
            nh_flank_region = new subsequence_classification_data_region(this, neighbourhood_flanking_atoms,load_ss_predictions,load_foldx_energy);
        }

        internal void init_nh_contacts(double nh_max_dist = 5.0, bool should_include_interface = false, bool load_ss_predictions = true, bool load_foldx_energy = true)
        {
            const string module_name = nameof(subsequence_classification_data);
            const string method_name = nameof(init_nh_contacts);

            if (nh_max_dist > 8.0)
            {
                throw new Exception($"{module_name}.{method_name}: the specified maximum atomic distance is too long for contacts");
            }

            var contact_table = interface_region.atoms.First().intramolecular_contact_table;

            var contact_indexes = new List<int>();
            foreach (var atom in interface_region.atoms)
            {
                var distances = contact_table[atom.intramolecular_contact_table_index];

                for (var i = 0; i < distances.Length; i++)
                {
                    if (i == atom.intramolecular_contact_table_index)
                    {
                        continue;
                    }

                    if (distances[i] <= nh_max_dist)
                    {
                        contact_indexes.Add(i);
                    }
                }
            }

            // get list of contacts
            var neighbourhood_3d_contacts = chain_region.atoms.Where(a => contact_indexes.Contains(a.intramolecular_contact_table_index)).ToList();
            neighbourhood_3d_contacts = should_include_interface ? neighbourhood_3d_contacts.Union(interface_region.atoms).ToList() : neighbourhood_3d_contacts.Except(interface_region.atoms).ToList();
            neighbourhood_3d_contacts = neighbourhood_3d_contacts.Distinct().SelectMany(a => a.amino_acid_atoms).Distinct().OrderBy(a => a.chain_id).ThenBy(a => a.residue_index).ThenBy(a => a.i_code).ToList();
            
            nh_contact_region=new subsequence_classification_data_region(this, neighbourhood_3d_contacts, load_ss_predictions, load_foldx_energy);
        }

        internal static List<string> get_row_comments_headers(subsequence_classification_data instance_meta_data)
        {
            var comment_headers = new List<string>
            {
                $@"row_index",
                $@"{nameof(instance_meta_data.pdb_id)}",
                $@"{nameof(instance_meta_data.chain_id)}", 
                $@"{nameof(instance_meta_data.dimer_type)}", 
                $@"{nameof(instance_meta_data.class_id)}", 
                $@"{nameof(instance_meta_data.class_name)}", 
                $@"{nameof(instance_meta_data.parallelism)}", 
                $@"{nameof(instance_meta_data.symmetry_mode)}"
            };

            foreach (var region in instance_meta_data.get_regions())
            {
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.unique_id)}");

                comment_headers.Add($@"{region.region_name}_{nameof(region.region.aa_sequence.Length)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.aa_sequence)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.res_ids)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.dssp3_monomer)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.dssp3_multimer)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.dssp_monomer)}");
                comment_headers.Add($@"{region.region_name}_{nameof(region.region.dssp_multimer)}");
                comment_headers.AddRange(region.region?.ss_predictions?.Select(a => $@"{region.region_name}_{a.format}").ToList() ?? new List<string>());
            }

            return comment_headers;
        }

        internal static List<string> get_row_comments(int row_index, subsequence_classification_data instance_meta_data)
        {
            var row_comments = new List<string>();

            var regions = instance_meta_data.get_regions();

            var if_region = regions.First(a => a.region_name == nameof(interface_region));

            row_comments.Add(row_index.ToString(CultureInfo.InvariantCulture));

            row_comments.Add($@"{instance_meta_data.pdb_id}");
            row_comments.Add($@"{instance_meta_data.chain_id}");
            row_comments.Add($@"{instance_meta_data.dimer_type}");
            row_comments.Add($@"{instance_meta_data.class_id:+#;-#;+0}");
            row_comments.Add($@"{instance_meta_data.class_name}");
            row_comments.Add($@"{instance_meta_data.parallelism}");
            row_comments.Add($@"{instance_meta_data.symmetry_mode}");

            foreach (var region in regions)
            {
                
                row_comments.Add($@"{region.region.unique_id()}");

                row_comments.Add($@"{region.region.aa_sequence.Length}");
                row_comments.Add($@"{region.region.aa_sequence}");
                row_comments.Add($@"{string.Join($@" ", region.region.res_ids.Select(a => $@"{a.amino_acid}{a.residue_index}{(a.i_code != default && !char.IsWhiteSpace(a.i_code) ? $@"{(char.IsDigit(a.i_code) ? /*"i"*/$@"" : $@"")}{a.i_code}" : $@"")}").ToList())}");
                row_comments.Add($@"{region.region.dssp3_monomer}");
                row_comments.Add($@"{region.region.dssp3_multimer}");
                row_comments.Add($@"{region.region.dssp_monomer}");
                row_comments.Add($@"{region.region.dssp_multimer}");
                row_comments.AddRange(region.region?.ss_predictions?.Select(a => a.prediction).ToList() ?? new List<string>());
            }

            return row_comments;
        }
    }
}
