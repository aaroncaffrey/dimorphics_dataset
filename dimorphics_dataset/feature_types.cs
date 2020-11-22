using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types
    {
        internal feature_types_1d feature_types_interface_1d;
        internal feature_types_1d feature_types_neighbourhood_1d;
        internal feature_types_1d feature_types_chain_1d;

        internal feature_types_2d feature_types_interface_2d;
        internal feature_types_2d feature_types_neighbourhood_2d;
        internal feature_types_2d feature_types_chain_2d;

        internal feature_types_3d feature_types_interface_3d;
        internal feature_types_3d feature_types_neighbourhood_3d;
        internal feature_types_3d feature_types_chain_3d;

        //internal void AsArray2d()
        //{
        //    var x= new[]{
        //        (nameof(feature_types_interface_2d), feature_types_interface_2d),
        //        (nameof(feature_types_neighbourhood_2d), feature_types_neighbourhood_2d),
        //        (nameof(feature_types_chain_2d), feature_types_chain_2d),
        //        (nameof(feature_types_interface_3d), feature_types_interface_3d),
        //        (nameof(feature_types_neighbourhood_3d), feature_types_neighbourhood_3d),
        //        (nameof(feature_types_chain_3d), feature_types_chain_3d)
        //           };
        //
        //    return x;
        //}

        internal string get_output_file_tag()
        {
            var x = new List<string>();

            // 1d tags
            if (feature_types_interface_1d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1i);
            }

            if (feature_types_neighbourhood_1d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1n);
            }

            if (feature_types_chain_1d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area1p);
            }

            // 2d tags
            if (feature_types_interface_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2i);
            }

            if (feature_types_neighbourhood_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2n);
            }

            if (feature_types_chain_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area2p);
            }
             
            // 3d tags
            if (feature_types_interface_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3i);
            }

            if (feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3n);
            }

            if (feature_types_chain_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add(cmd_params.area3p);
            }

            return String.Join("", x);
        }

        internal static feature_types feature_types_params(string[] area)
        {
            var module_name = nameof(dimorphics_dataset.feature_types);
            var method_name = nameof(feature_types_params);

            if (area == null || area.Length == 0)
            {
                return null;
            }

            var do_1d_interface = area?.Any(a => string.Equals(a, cmd_params.area1i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_1d_nh = area?.Any(a => string.Equals(a, cmd_params.area1n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_1d_protein = area?.Any(a => string.Equals(a, cmd_params.area1p, StringComparison.InvariantCultureIgnoreCase)) ?? false;

            var do_2d_interface = area?.Any(a => string.Equals(a, cmd_params.area2i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_2d_nh = area?.Any(a => string.Equals(a, cmd_params.area2n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_2d_protein = area?.Any(a => string.Equals(a, cmd_params.area2p, StringComparison.InvariantCultureIgnoreCase)) ?? false;

            var do_3d_interface = area?.Any(a => string.Equals(a, cmd_params.area3i, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_3d_nh = area?.Any(a => string.Equals(a, cmd_params.area3n, StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_3d_protein = area?.Any(a => string.Equals(a, cmd_params.area3p, StringComparison.InvariantCultureIgnoreCase)) ?? false;


            var feature_types = new feature_types()
            {
                // 1d
                feature_types_interface_1d = !do_1d_interface
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_1d_interface,
                        sequence_geometry_classification_data = do_1d_interface,
                        //mpsa_classification_data_subsequence = do_1d_interface,
                        intrinsically_unordered_data = do_1d_interface,
                        aa_index_classification_data = do_1d_interface,
                        sable_classification_data = do_1d_interface,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_1d_interface,
                        r_peptides = do_1d_interface,
                        r_protr = do_1d_interface,
                    },
                feature_types_neighbourhood_1d = !do_1d_nh
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_1d_nh,
                        sequence_geometry_classification_data = do_1d_nh,
                        //mpsa_classification_data_subsequence = do_1d_nh,
                        intrinsically_unordered_data = do_1d_nh,
                        aa_index_classification_data = do_1d_nh,
                        sable_classification_data = do_1d_nh,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_1d_nh,
                        r_peptides = do_1d_nh,
                        r_protr = do_1d_nh,
                    },
                feature_types_chain_1d = !do_1d_protein
                    ? null
                    : new feature_types_1d()
                    {
                        pse_aac_sequence_classification_data = do_1d_protein,
                        sequence_geometry_classification_data = do_1d_protein,
                        //mpsa_classification_data_subsequence = do_1d_protein,
                        intrinsically_unordered_data = do_1d_protein,
                        aa_index_classification_data = do_1d_protein,
                        sable_classification_data = do_1d_protein,
                        dna_binding_prediction_data = do_1d_protein,
                        blast_pssm_subsequence_classification_data = do_1d_protein,
                        r_peptides = do_1d_protein,
                        r_protr = do_1d_protein,
                    },

                // 2d
                feature_types_interface_2d = !do_2d_interface
                    ? null
                    : new feature_types_2d()
                    {
                        mpsa_classification_data_subsequence = do_2d_interface,
                    },
                feature_types_neighbourhood_2d = !do_2d_nh
                    ? null
                    : new feature_types_2d()
                    {
                        
                        mpsa_classification_data_subsequence = do_2d_nh,
                        
                    },
                feature_types_chain_2d = !do_2d_protein
                    ? null
                    : new feature_types_2d()
                    {
                        mpsa_classification_data_subsequence = do_2d_protein,
                        
                    },


                // 3d

                feature_types_interface_3d = !do_3d_interface
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_interface,
                        intramolecular_classification_data = do_3d_interface,

                        foldx_classification_data = do_3d_interface,
                        tortuosity_classification_data = do_3d_interface,
                        ring_classification_data = do_3d_interface,
                        pse_ssc_dssp_classification_data = do_3d_interface,
                    },
                feature_types_neighbourhood_3d = !do_3d_nh
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_nh,
                        intramolecular_classification_data = do_3d_nh,

                        foldx_classification_data = do_3d_nh,
                        tortuosity_classification_data = do_3d_nh,
                        ring_classification_data = do_3d_nh,
                        pse_ssc_dssp_classification_data = do_3d_nh,
                    },
                feature_types_chain_3d = !do_3d_protein
                    ? null
                    : new feature_types_3d()
                    {
                        sasa_classification_data = do_3d_protein,
                        intramolecular_classification_data = do_3d_protein,

                        foldx_classification_data = do_3d_protein,
                        tortuosity_classification_data = do_3d_protein,
                        ring_classification_data = do_3d_protein,
                        pse_ssc_dssp_classification_data = do_3d_protein,
                    }
            };

            // 1d
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_interface_1d)} = {feature_types.feature_types_interface_1d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_1d)} = {feature_types.feature_types_neighbourhood_1d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_1d)} = {feature_types.feature_types_chain_1d}", module_name, method_name);

            //2d
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_interface_2d)} = {feature_types.feature_types_interface_2d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_2d)} = {feature_types.feature_types_neighbourhood_2d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_2d)} = {feature_types.feature_types_chain_2d}", module_name, method_name);

            //3d
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_interface_3d)} = {feature_types.feature_types_interface_3d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}", module_name, method_name);
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_3d)} = {feature_types.feature_types_chain_3d}", module_name, method_name);
            return feature_types;
        }
    }
}
