using System.Collections.Generic;

namespace dimorphics_dataset
{
    //internal class foldx_energy_differences
    //{
    //    internal List<(bool is_fragment, bool is_repaired, double[] energy)> ala_scan;
    //    internal List<(bool is_fragment, bool is_repaired, double[] energy)> build_model_position_scan;
    //    internal List<(bool is_fragment, bool is_repaired, double[] energy)> build_model_subsequence_mutant;
    //    internal List<(bool is_fragment, bool is_repaired, double[] energy)> position_scan;
    //    //internal List<(bool is_fragment, bool is_repaired, double[] energy)> stability;
    //}

    internal class foldx_energy_differences
    {
        //internal (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) foldx_ala_scanning_result_protein;

        internal (string cmd_line, string wait_filename, List<foldx_ala_scanning_result> data) foldx_ala_scanning_result_subsequence;

        internal (string cmd_line, string wait_filename, List<foldx_position_scanning_result> data) foldx_position_scanning_result_subsequence;

        internal (string cmd_line, string wait_filename, List<foldx_energy_terms_ps> data) foldx_buildmodel_position_scan_result_subsequence;

        internal (string cmd_line, string wait_filename, List<foldx_energy_terms_sm> data) foldx_buildmodel_subsequence_mutant_result_subsequence;
    }
}
