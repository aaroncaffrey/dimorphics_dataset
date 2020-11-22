using System.Collections.Generic;
using System.Linq;

namespace dimorphics_dataset
{
    internal class feature_types_2d
    {
        
        //internal bool sable_classification_data = false;
        internal bool mpsa_classification_data_subsequence;
        

        internal feature_types_2d()
        {

        }
        internal feature_types_2d(bool enable)
        {
        
            //sable_classification_data = enable;
            mpsa_classification_data_subsequence = enable;
        
        }

        internal List<(string key, bool value)> AsArray()
        {
            var data = new List<(string key, bool value)>()
                {
        
                    //( nameof(sable_classification_data                    ),  sable_classification_data                         ) ,
                    ( nameof(mpsa_classification_data_subsequence         ),  mpsa_classification_data_subsequence              ) ,
        
                };

            return data;
        }

        internal string key_value_list_hr()
        {
            var data = AsArray();

            var ret = $@"({string.Join(", ", data.Select(a => $"{a.key} = {a.value}").ToList())})";

            return ret;
        }

        public override string ToString()
        {
            return key_value_list_hr();
            //throw new NotImplementedException();
        }
    }
}
