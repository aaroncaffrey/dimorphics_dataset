namespace dimorphics_dataset
{

    public class info_blast_pssm_options
    {
        internal bool make_unsplit_sequence = true;
        internal bool make_split_sequence = true;
        
        internal bool make_standard_encoding = true;
        internal bool make_distance_transform = true; // false
        
        internal bool normalise_none = false;
        internal bool normalise_whole_pssm = true;
        internal bool normalise_subsequence = true;// false
        internal bool normalise_encoded_parts = true;// false
        internal bool normalise_encoded_vector = true;// false

        internal bool encode_standard_vector = true;
        internal bool encode_distance_vector = true; // false
        internal bool encode_interval_vector = true; // false

        // separate groups with the following statistics:
        internal bool encode_min = false;
        internal bool encode_max = false;
        internal bool encode_mean = true;
        internal bool encode_mean_sd = false;
        internal bool encode_median = false;
        internal bool encode_mode = false;
        internal bool encode_range = false;

        // encoding size (when using standard AA alphabet) - e.g. 20x20 = 400 could be 4x4 = 16 for reduced alphabet
        internal bool size_1 = true;
        internal bool size_20 = true;
        internal bool size_210 = true;
        internal bool size_400 = true;
        
        internal bool db_nr_local_1e_4 = false; // 0.0001
        internal bool db_nr_local_def = true; 
        internal bool db_nr_remote_1e_4 = false;
        
        internal bool db_sp_local_1e_4 = false;
        internal bool db_sp_local_def = true;
        internal bool db_sp_remote_1e_4 = false;
        
        internal bool db_ur90_local_def = true;


        public void set_all(bool value)
        {
            make_unsplit_sequence = value;
            make_split_sequence = value;

            make_standard_encoding = value;
            make_distance_transform = value;

            normalise_none = value;
            normalise_whole_pssm = value;
            normalise_subsequence = value;
            normalise_encoded_parts = value;
            normalise_encoded_vector = value;

            encode_standard_vector = value;
            encode_distance_vector = value;
            encode_interval_vector = value;


            encode_min = value;
            encode_max = value;
            encode_mean = value;
            encode_mean_sd = value;
            encode_median = value;
            encode_mode = value;
            encode_range = value;


            size_1 = value;
            size_20 = value;
            size_210 = value;
            size_400 = value;

            db_nr_local_1e_4 = value;
            db_nr_local_def = value;
            db_nr_remote_1e_4 = value;

            db_sp_local_1e_4 = value;
            db_sp_local_def = value;
            db_sp_remote_1e_4 = value;

            db_ur90_local_def = value;
        }
    }

}