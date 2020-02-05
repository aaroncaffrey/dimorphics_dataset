using System;
using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    public class feature_types
    {
        internal feature_types_2d feature_types_interface_2d;
        internal feature_types_2d feature_types_neighbourhood_2d;
        internal feature_types_2d feature_types_chain_2d;
        internal feature_types_3d feature_types_interface_3d;
        internal feature_types_3d feature_types_neighbourhood_3d;
        internal feature_types_3d feature_types_chain_3d;

        //public void AsArray2d()
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

        public string get_output_file_tag()
        {
            var x = new List<string>();

            if (feature_types_interface_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("2i");
            }

            if (feature_types_neighbourhood_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("2n");
            }

            if (feature_types_chain_2d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("2p");
            }

            if (feature_types_interface_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("3i");
            }

            if (feature_types_neighbourhood_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("3n");
            }

            if (feature_types_chain_3d?.AsArray()?.Any(a => a.value) ?? false)
            {
                x.Add("3p");
            }

            return String.Join("", x);
        }

        public static feature_types feature_types_params(string[] area)
        {
            var do_2d_interface = area?.Any(a => string.Equals(a, "2i", StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_2d_nh = area?.Any(a => string.Equals(a, "2n", StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_2d_protein = area?.Any(a => string.Equals(a, "2p", StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_3d_interface = area?.Any(a => string.Equals(a, "3i", StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_3d_nh = area?.Any(a => string.Equals(a, "3n", StringComparison.InvariantCultureIgnoreCase)) ?? false;
            var do_3d_protein = area?.Any(a => string.Equals(a, "3p", StringComparison.InvariantCultureIgnoreCase)) ?? false;


            var feature_types = new feature_types()
            {
                feature_types_interface_2d = !do_2d_interface
                    ? null
                    : new feature_types_2d()
                    {
                        pse_aac_sequence_classification_data = do_2d_interface,
                        sequence_geometry_classification_data = do_2d_interface,
                        mpsa_classification_data_subsequence = do_2d_interface,
                        intrinsically_unordered_data = do_2d_interface,
                        aa_index_classification_data = do_2d_interface,
                        sable_classification_data = do_2d_interface,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_2d_interface,
                        r_peptides = do_2d_interface,
                        r_protr = do_2d_interface,
                    },
                feature_types_neighbourhood_2d = !do_2d_nh
                    ? null
                    : new feature_types_2d()
                    {
                        pse_aac_sequence_classification_data = do_2d_nh,
                        sequence_geometry_classification_data = do_2d_nh,
                        mpsa_classification_data_subsequence = do_2d_nh,
                        intrinsically_unordered_data = do_2d_nh,
                        aa_index_classification_data = do_2d_nh,
                        sable_classification_data = do_2d_nh,
                        dna_binding_prediction_data = false, //must be false - protein level only
                        blast_pssm_subsequence_classification_data = do_2d_nh,
                        r_peptides = do_2d_nh,
                        r_protr = do_2d_nh,
                    },
                feature_types_chain_2d = !do_2d_protein
                    ? null
                    : new feature_types_2d()
                    {
                        pse_aac_sequence_classification_data = do_2d_protein,
                        sequence_geometry_classification_data = do_2d_protein,
                        mpsa_classification_data_subsequence = do_2d_protein,
                        intrinsically_unordered_data = do_2d_protein,
                        aa_index_classification_data = do_2d_protein,
                        sable_classification_data = do_2d_protein,
                        dna_binding_prediction_data = do_2d_protein,
                        blast_pssm_subsequence_classification_data = do_2d_protein,
                        r_peptides = do_2d_protein,
                        r_protr = do_2d_protein,
                    },
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

            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_interface_2d)} = {feature_types.feature_types_interface_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_neighbourhood_2d)} = {feature_types.feature_types_neighbourhood_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_2d)} = {feature_types.feature_types_chain_2d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_interface_3d)} = {feature_types.feature_types_interface_3d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine(
                $@"{nameof(feature_types.feature_types_neighbourhood_3d)} = {feature_types.feature_types_neighbourhood_3d}",
                nameof(program), nameof(feature_types_params));
            io_proxy.WriteLine($@"{nameof(feature_types.feature_types_chain_3d)} = {feature_types.feature_types_chain_3d}",
                nameof(program), nameof(feature_types_params));
            return feature_types;
        }
    }
}
