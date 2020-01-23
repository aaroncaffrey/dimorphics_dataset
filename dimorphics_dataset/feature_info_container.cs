using System;
using System.Collections.Generic;
using System.Text;
using Newtonsoft.Json;

namespace dimorphics_dataset
{
    public class feature_info_container
    {
        public List<feature_info> feautre_info_list;

        public static string serialise_json(feature_info_container feature_info_container)
        {
            var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All, };
            var serialised_serialise_json = JsonConvert.SerializeObject(feature_info_container, json_settings);
            return serialised_serialise_json;
        }

        public static feature_info_container deserialise(string serialized_json)
        {
            var json_settings = new JsonSerializerSettings() { PreserveReferencesHandling = PreserveReferencesHandling.All };
            feature_info_container feature_info_container;

            try
            {
                feature_info_container = JsonConvert.DeserializeObject<feature_info_container>(serialized_json, json_settings);
            }
            catch (Exception e)
            {
                feature_info_container = null;

                io_proxy.WriteLine(e.ToString());
            }

            return feature_info_container;
        }
    }
}
